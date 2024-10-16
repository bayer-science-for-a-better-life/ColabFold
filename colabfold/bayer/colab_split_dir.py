"""ColabSplit. Split a collection of alignment files to a defined number of directories.

This enables better usage of the GPU for predictions based on a single FASTA file.
"""
import argparse
import logging
import numpy as np
import shutil
from pathlib import Path


def batch(iterable, n=1):
    """
    Splits an iterable into batches of size n.

    Args:
        iterable (iterable): The iterable to split.
        n (int): The size of each batch.

    Yields:
        iterable: Batches of the original iterable.
    """
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]


def split_a3m_files(input_dir: str, output_dir: str, number_of_splits: int):
    """
    Splits .a3m files from the input directory into a specified number of subdirectories in the output directory.

    Args:
        input_dir (str): The directory containing the .a3m files to split.
        output_dir (str): The directory where the split files will be stored.
        number_of_splits (int): The number of subdirectories to create.
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    a3m_files = list(input_dir.glob('*.a3m'))
    files_per_split = int(np.ceil(len(a3m_files) / number_of_splits))
    logging.info(f"Input Dir: {input_dir}")
    logging.info(f"Output Dir: {output_dir}")
    logging.info(f"Number of files: {len(a3m_files)}")
    logging.info(f"Number of desired splits (jobs): {number_of_splits}")
    logging.info(f"Number of files per split: {files_per_split}")

    for job_count, a3m_files_to_move in enumerate(batch(a3m_files, files_per_split)):
        gpu_output_dir = Path(output_dir) / f'job_{str(job_count).zfill(2)}'

        for a3m_file in a3m_files_to_move:
            logging.info(f"Moving {a3m_file} to {gpu_output_dir}")
            output_file_path = gpu_output_dir
            output_file_path.mkdir(parents=True, exist_ok=True)
            shutil.copy(a3m_file, output_file_path)

    # verify that the number of files matches
    output_files = list(output_dir.rglob('*.a3m'))
    output_files = set([i.stem for i in output_files])
    a3m_files = set([i.stem for i in a3m_files])
    missing_files = a3m_files - output_files
    logging.info(f"Missing Files: {missing_files}")
    logging.info(f"Number of Input Files: {len(output_files)}")
    logging.info(f"Number of Output Files: {len(a3m_files)}")
    logging.info(f"Number of Missing Files: {len(missing_files)}")
    assert output_files == a3m_files, f"Files not copied correctly. {output_files} != {a3m_files}"


if __name__ == '__main__':
    # Define the argument parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input_dir', type=str, help='Input directory containing .a3m files')
    parser.add_argument('--output_dir', type=str, help='Output directory for split files')
    parser.add_argument('--n_splits', type=str, help='Number of splits (subdirectories) to create')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(funcName)s(): %(message)s",
        handlers=[
            logging.StreamHandler(),
        ],
    )

    input_dir = args.input_dir
    output_dir = args.output_dir
    n_splits = int(args.n_splits)
    split_a3m_files(input_dir, output_dir, n_splits)
    logging.info('Completed in splitting')
