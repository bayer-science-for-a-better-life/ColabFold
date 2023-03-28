"""
WIP: potential use case: create a manually designed MSA with paired and unpaired sequences.
"""
import logging
import logging.config
import numpy as np
import os
import random
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple, Sequence


def _to_a3m(sequences: Sequence[str]) -> str:
    """Converts sequences to an a3m file."""
    names = ['sequence %d' % i for i in range(1, len(sequences) + 1)]
    a3m = []
    for sequence, name in zip(sequences, names):
        a3m.append(u'>' + name + u'\n')
        a3m.append(sequence + u'\n')
    return ''.join(a3m)


class Kalign:
    """
    Python wrapper of the Kalign binary. Source.

    AlphaFold:
    https://huggingface.co/spaces/simonduerr/ProteinMPNN/blob/main/af_backprop/alphafold/data/tools/kalign.py
    """

    def __init__(self, *, binary_path: str, gap_open: int = 11, gap_extension=0.85, terminal_gap=0.45):
        """Initializes the Python Kalign wrapper.
        Args:
          binary_path: The path to the Kalign binary.
        Raises:
          RuntimeError: If Kalign binary not found within the path.
        """
        self.binary_path = binary_path
        self.gap_open = gap_open
        self.gap_extension = gap_extension
        self.terminal_gap = terminal_gap

    def align(self, sequences: Sequence[str]) -> str:
        """Aligns the sequences and returns the alignment in A3M string.
        Args:
          sequences: A list of query sequence strings. The sequences have to be at
            least 6 residues long (Kalign requires this). Note that the order in
            which you give the sequences might alter the output slightly as
            different alignment tree might get constructed.
        Returns:
          A string with the alignment in a3m format.
        Raises:
          RuntimeError: If Kalign fails.
          ValueError: If any of the sequences is less than 6 residues long.
        """
        logging.info('Aligning %d sequences', len(sequences))

        for s in sequences:
            if len(s) < 6:
                raise ValueError('Kalign requires all sequences to be at least 6 '
                                 'residues long. Got %s (%d residues).' % (s, len(s)))

        # with utils.tmpdir_manager(base_dir='/tmp') as query_tmp_dir:
        with tempfile.TemporaryDirectory() as query_tmp_dir:
            input_fasta_path = os.path.join(query_tmp_dir, 'input.fasta')
            output_a3m_path = os.path.join(query_tmp_dir, 'output.a3m')

            with open(input_fasta_path, 'w') as f:
                f.write(_to_a3m(sequences))

            cmd = [
                self.binary_path,
                '-q',
                '-i', input_fasta_path,
                '-o', output_a3m_path,
                '-c', 'input'
                      '-format', 'fasta',
                '-gapopen', str(self.gap_open),
                '-gapextension', str(self.gap_extension),
                '-tgpe', str(self.terminal_gap),

            ]

            logging.info('Launching subprocess "%s"', ' '.join(cmd))
            _ = subprocess.check_output(cmd)

            with open(output_a3m_path) as f:
                a3m = f.read()

            return a3m


def sample_sequence(seqs: List[List[str]], n: int, seed=0, use_top: bool = True, use_weights: bool = True):
    """
    Sample n sequences from a list of sequences.

    Either sample from the top (not really sampling) or truly sample from all the sequences with equal probability.

    Parameters:
        seqs: List of tuples (name, sequence)
        n: Number of sequences to sample
        seed: Random seed
        use_top: If True, sample from the top n sequences. If False, sample from all sequences with equal probability.
    """
    if use_top:
        return seqs[:n]
    else:
        rng = random.Random(seed)

        if use_weights:
            weights = [1 / (i + 1) for i in range(len(seqs))]

            return rng.choices(seqs, weights=weights, k=n)
        else:
            return rng.sample(seqs, n)


def pair_sequences(sampled_pairs_query1, sampled_pairs_query2, seed: int = 0) -> Tuple[List[str], List[str]]:
    """
    Take the results from the alignments / sampling runs and pair them randomly.

    This is not optimal since the OG AF paper only pairs sequences by organism. => experimental

    Parameters:
        sampled_pairs_query1: List of tuples (name, sequence) for the first query
        sampled_pairs_query2: List of tuples (name, sequence) for the second query
        seed: Random seed
    """
    seq1_copy = sampled_pairs_query1.copy()
    seq2_copy = sampled_pairs_query2.copy()

    rng_1 = random.Random(seed)
    rng_1.shuffle(seq1_copy)

    rng_2 = random.Random(seed)
    rng_2.shuffle(seq2_copy)

    return seq1_copy, seq2_copy


def concat_sequences(seqs1, seqs2, complex_aln_length):
    """
    Concatenate two lists of sequences.

    Parameters:
        seqs1: List of tuples (name, sequence) for the first query
        seqs2: List of tuples (name, sequence) for the second query
    """
    joined_seqs = []
    for seq1_tuple, seq2_tuple in zip(seqs1, seqs2):
        header_1 = seq1_tuple[0].split("\t")
        header_2 = seq2_tuple[0].split("\t")
        new_header = header_1[0] + "_" + header_2[0].replace(">", "")
        seq1 = seq1_tuple[1]
        seq2 = seq2_tuple[1]
        assert len(seq1[1]) + len(seq2[1]) == complex_aln_length
    seqs1, seqs2 = sampled_rnd_pairs_query1, sampled_rnd_pairs_query2
    return seqs1 + seqs2


def get_sequence_order(ref_seq, *args):
    return np.array([ref_seq.index(seq) for seq in args]).argsort()


def parse_a3m(a3m_file):
    with open(a3m_file, "r") as f:
        lines = f.readlines()
    seqs = []
    seq_tuple = []
    for line in lines:
        if line.startswith(">"):
            seq_tuple.append(line.strip())
        else:
            seq_tuple.append(line.strip())
            seqs.append(seq_tuple)
            seq_tuple = []
    return seqs


if __name__ == '__main__':
    # init logging
    logging.basicConfig(level=logging.INFO)
    # define the path to the a3m files & number of proteins in complex
    res_dir = Path("<insert_path_her>")
    n_proteins = 3

    max_paired = 500
    max_unpaired = 500
    max_unpaired_share = int(max_unpaired * 0.5)

    # sequence objects
    seq_data = {f"seq{i}": f"{i}.a3m" for i in range(n_proteins)}
    seq1_a3m = res_dir / seq_data["seq0"]
    seq2_a3m = res_dir / seq_data["seq1"]
    seq3_a3m = res_dir / seq_data["seq2"]

    seqs_msa1 = parse_a3m(seq1_a3m)[1:]
    seqs_msa2 = parse_a3m(seq2_a3m)[1:]
    seqs_msa3 = parse_a3m(seq3_a3m)[1:]

    que_seq1 = seqs_msa1[0][1]
    que_seq2 = seqs_msa2[0][1]
    ref_seq = seqs_msa3[0][1]

    # get order of sequences in the paired MSA
    seq_order = get_sequence_order(ref_seq, que_seq1, que_seq2)
    logging.info(f"Sequence order: {seq_order}")
    aln_msa1 = len(seqs_msa1)
    aln_msa2 = len(seqs_msa2)
    aln_msa3 = len(seqs_msa3)

    max_unpaired_share_q1 = min(max_unpaired_share, aln_msa1)
    max_unpaired_share_q2 = min(max_unpaired_share, aln_msa2)
    max_unpaired_share = min(max_paired, max_unpaired_share_q1, max_unpaired_share_q2)

    logging.info(f"Alignments: Paired: {aln_msa3:_}, Unpaired #1: {aln_msa1:_}, Unpaired #2: {aln_msa2:_}")
    logging.info(seqs_msa1[0])
    logging.info(seqs_msa2[0])
    logging.info(seqs_msa3[0])

    len_query_complex = len(seqs_msa3[0][1])
    len_query1 = len(seqs_msa1[0][1])
    len_query2 = len(seqs_msa2[0][1])
    logging.info(
        f"Lengths: Paired: {len_query_complex:_}, Unpaired #1: {len_query1:_}, Unpaired #2: {len_query2:_} (Sum: {len_query1 + len_query2:_})")

    complex_aln_length = len(seqs_msa3[0][1])
    # sample query 1 & 2 individually
    sampled_pairs_query1 = sample_sequence(seqs_msa1, max_unpaired_share_q1, seed=0, use_top=True)
    sampled_pairs_query2 = sample_sequence(seqs_msa2, max_unpaired_share_q2, seed=0, use_top=True)
    sampled_rnd_pairs_query1, sampled_rnd_pairs_query2 = \
        pair_sequences(sampled_pairs_query1, sampled_pairs_query2, seed=0)

    # write a function that randomly pairs the sequences from the two query sequences
    new_msa = [ref_seq]
    for i in range(len(sampled_rnd_pairs_query1)):
        new_msa.append(sampled_rnd_pairs_query1[i][1] + sampled_rnd_pairs_query2[i][1])

    # write a function that samples from a list of sequences but with decreasing probability from top to bottom
    sampled_alns_query1 = sample_sequence(seqs_msa1, max_unpaired_share_q1, seed=0, use_weights=True)
    sampled_alns_query2 = sample_sequence(seqs_msa2, max_unpaired_share_q2, seed=0, use_weights=True)

    # for i in range(len(sampled_alns_query1)):
    #     new_msa.append(sampled_alns_query1[i][1])
    #
    # for i in range(len(sampled_alns_query2)):
    #     new_msa.append(sampled_alns_query2[i][1])
    #
    # aln_lengths_query1 = [len(seq[1]) for seq in seqs_msa1]
    # aln_lengths_query2 = [len(seq[1]) for seq in seqs_msa2]
    # aln_lengths_query3 = [len(seq[1]) for seq in seqs_msa3]

    # logging.info(f"New MSA length: {len(new_msas):_}")
    # logging.info(f"MSA lengths: ")
    # logging.info(f"Paired: {len(sampled_pairs_complex):_}")
    # logging.info(f"Unpaired #1: {len(sampled_pairs_query1):_}")
    # logging.info(f"Unpaired #2: {len(sampled_pairs_query2):_}")

    aligner = Kalign(binary_path="kalign")
    a3m_res = aligner.align(new_msa[0:3])
    print("done")
    with open(res_dir / "kalign_test.a3m", "w") as f:
        f.write(a3m_res)

    print(a3m_res)

    # dummy_ref = que_seq1 + "W" * 10 + que_seq2
    # aligner = Align.PairwiseAligner()
    # aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    # aligner.open_gap_score = -10
    # aligner.mode = "local"
    # alignments = aligner.align(new_msa[0], sampled_alns_query1[50][1].replace("-", "").upper())[0]
    # print(alignments)
    # alignments = aligner.align(dummy_ref, sampled_alns_query1[50][1].replace("-", "").upper())[0]
    # print(alignments)
    # alignments = aligner.align(sampled_alns_query1[0][1], sampled_alns_query1[50][1].replace("-", "").upper())[0]
    # print(alignments)
    #
    # aligner.score(sampled_alns_query1[0][1].upper(), sampled_alns_query1[50][1].upper())
