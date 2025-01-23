import pandas as pd
import scipy.stats as stats
import numpy as np
import argparse 

def read_args():
    parser = argparse.ArgumentParser(description="Calculate median discrepancies and perform Wilcoxon test.")
    parser.add_argument("--core", required=True, help="Path to the core CSV file.")
    parser.add_argument("--accessory", required=True, help="Path to the accessory CSV file.")
    parser.add_argument("--kmers", required=True, help="Path to the kmers text file.")

    return parser.parse_args()

# Load data from core and accessory files
def load_csv(file):
    return pd.read_csv(file)

# Load kmers from kmers.txt
def load_kmers(file):
    with open(file, 'r') as f:
        return set(line.strip().lower() for line in f)

def calculate_discrepancy_and_test(core_file, accessory_file, kmers_file):
    # Load core and accessory data
    core = load_csv(core_file)
    accessory = load_csv(accessory_file)

    # Load k-mers
    kmers = load_kmers(kmers_file)

    # Merge core and accessory on "word"
    merged = pd.merge(core, accessory, on='word', suffixes=('_core', '_accessory'))

    # Calculate score discrepancy
    merged['score_discrepancy'] = merged['score_core'] - merged['score_accessory']

    # Separate into k-mers and other words
    kmers_data = merged[merged['word'].isin(kmers)]
    other_data = merged[~merged['word'].isin(kmers)]

    # Calculate median discrepancies
    median_discrepancy_kmers = kmers_data['score_discrepancy'].median() if not kmers_data.empty else None
    median_discrepancy_others = other_data['score_discrepancy'].median() if not other_data.empty else None

    # Perform Wilcoxon rank-sum test
    if not kmers_data.empty and not other_data.empty:
        wilcox_stat, p_value = stats.ranksums(kmers_data['score_discrepancy'], other_data['score_discrepancy'])
    else:
        wilcox_stat, p_value = None, None

    return {
        'median_discrepancy_kmers': median_discrepancy_kmers,
        'median_discrepancy_others': median_discrepancy_others,
        'wilcoxon_statistic': wilcox_stat,
        'p_value': p_value
    }

# Example usage
if __name__ == "__main__":
    args = read_args()
    results = calculate_discrepancy_and_test(args.core, args.accessory, args.kmers)
    filename_base = args.core.split('_')[0].rsplit("/", 1)[-1]

    results = [filename_base,results['median_discrepancy_kmers'], results['median_discrepancy_others'], results['wilcoxon_statistic'],results['p_value']]
    print(','.join([str(x) for x in results]))
