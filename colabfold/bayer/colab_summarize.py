"""ColabSummarize.

Collect the iptm and ptm metrics from the structure prediction jobs."""
import argparse
import json
import pandas as pd
from pathlib import Path
from tqdm import tqdm


def extract_ptm_iptm(json_file):
    """Parse the json file from ColabFold."""
    with open(json_file, 'r') as file:
        data = json.load(file)

    if 'iptm+ptm' in data:
        ptm = data['iptm+ptm']
        iptm = data['iptm']
    else:
        ptm = data['ptm']
        iptm = data['iptm']
        max_pae = data['max_pae']

    return {'ptm': ptm, 'iptm': iptm, "max_pae": max_pae}


def process_json_files(source_directory, destination_directory, csv_name):
    """Process json files."""
    destination_path = Path(destination_directory)
    destination_path.mkdir(parents=True, exist_ok=True)

    pattern = '*scores*.json'
    matching_files = list(Path(source_directory).rglob(pattern))

    result_data = pd.DataFrame(columns=['Complex', 'pTM', 'ipTM', 'max_pae', 'Ranking_confidence'])

    for json_file_path in tqdm(matching_files):
        # get metrics
        ptm_iptm_values = extract_ptm_iptm(json_file_path)
        ptm_value = ptm_iptm_values['ptm']
        iptm_value = ptm_iptm_values['iptm']
        max_pae_value = ptm_iptm_values['max_pae']

        # get name from file
        # we expect this scheme:
        # query1_and_query2_scores_rank_004_alphafold2_multimer_v3_model_1_seed_000.json
        protein_names = Path(json_file_path).stem.split("_scores")[0]

        file_data = pd.DataFrame({
            'Complex': [protein_names],
            'pTM': [ptm_value],
            'ipTM': [iptm_value],
            'max_pae': [max_pae_value],
            'Ranking_confidence': [0.2 * ptm_value + 0.8 * iptm_value]
        })

        result_data = pd.concat([result_data, file_data], ignore_index=True)
    result_data = result_data.sort_values(by=["Ranking_confidence"], ascending=False)
    result_data["rank"] = range(1, len(result_data) + 1)
    result_data = result_data.reset_index(drop=True)
    csv_file_path = destination_path / f"{csv_name}.csv"
    result_data.to_csv(csv_file_path, index=False)


def main():
    parser = argparse.ArgumentParser(description='Process JSON files.')
    parser.add_argument('--source_directory', type=str, default='source_directory_path',
                        help='Source directory containing the JSON files')
    parser.add_argument('--destination_directory', type=str, default='destination_directory_path',
                        help='Destination directory for JSON files and analysis')
    parser.add_argument('--csv_name', type=str, default='Output',
                        help='Name of the output CSV file')
    args = parser.parse_args()

    process_json_files(args.source_directory, args.destination_directory, args.csv_name)


if __name__ == '__main__':
    main()
