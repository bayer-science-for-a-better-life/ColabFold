"""ColabPulldown. Simple script to concatenate sequences from two fasta files.

Provide the bait.fasta (your screening sequence) and the candidates.fasta (your library) to concatenate the sequences.

Note that the header is cut after the first space and special characters are replaced by "_".
"""
import argparse
import logging
import numpy as np
import time
from Bio import SeqIO


# Function to replace special characters in protein names
def replace_special_characters(sequence):
    """
    Replaces special characters in a protein sequence with underscores.

    Args:
        sequence (str): The protein sequence to process.

    Returns:
        str: The processed protein sequence with special characters replaced by underscores.
    """
    return sequence.replace(' ', '_').replace('*', '_').replace('|', '_').replace(':', '_')


if __name__ == '__main__':
    # Define the argument parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bait', type=str, help='Input fasta file for bait sequence')
    parser.add_argument('--candidates', type=str, help='Input fasta file for candidate sequences')
    parser.add_argument('--output', type=str, help='Output fasta file for concatenated sequences')
    parser.add_argument('--reciprocal', action='store_true', help='Combine sequences in both directions')
    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(funcName)s(): %(message)s",
        handlers=[
            logging.StreamHandler(),
        ],
    )

    bait = args.bait
    candidates = args.candidates
    output_file = args.output

    # Log the start time
    start_time = time.perf_counter()
    logging.info(f'Starting the bioinformatics analysis at {start_time}')

    # Parse the bait fasta file
    bait_records = list(SeqIO.parse(bait, 'fasta'))
    logging.info(f'Number of entries in bait fasta: {len(bait_records)}')

    # Parse the candidates fasta file
    candidate_records = list(SeqIO.parse(candidates, 'fasta'))
    logging.info(f'Number of entries in candidates fasta: {len(candidate_records)}')

    # Replace special characters in protein names
    for record in bait_records:
        record.id = replace_special_characters(record.id)

    for record in candidate_records:
        record.id = replace_special_characters(record.id)

    # Write a new fasta file through concatenation of the candidate sequences with the bait
    counter = 0
    with open(output_file, 'w') as output_handle:
        for candidate_record in candidate_records:
            concatenated_sequence = str(bait_records[0].seq) + ':' + str(candidate_record.seq)
            output_handle.write(f'>{bait_records[0].id}_and_{candidate_record.id}\n{concatenated_sequence}\n')
            counter += 1

            if args.reciprocal:
                concatenated_sequence = str(candidate_record.seq) + ':' + str(bait_records[0].seq)
                output_handle.write(f'>{candidate_record.id}_and_{bait_records[0].id}\n{concatenated_sequence}\n')
                counter += 1

    logging.info(f'Wrote {counter} concatenated sequences to {output_file}')

    # Log the end time and the output file written
    end_time = time.perf_counter()
    duration = np.round((end_time - start_time) / 60, 2)
    logging.info(f'Completed in {duration}. Output file written to {output_file}')
