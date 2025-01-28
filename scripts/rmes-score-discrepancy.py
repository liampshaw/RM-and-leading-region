import pandas as pd
import scipy.stats as stats
import numpy as np
import argparse 

def read_args():
    parser = argparse.ArgumentParser(description="Calculate median discrepancies and perform Wilcoxon test.")
    parser.add_argument("--file1", required=True, help="Path to the first CSV file (e.g. core genes).")
    parser.add_argument("--file2", required=True, help="Path to the other CSV file (e.g. accessory genes).")
    parser.add_argument("--kmers", required=True, help="Path to the kmers text file.")

    return parser.parse_args()

# Load data from a csv file
def load_csv(file):
    return pd.read_csv(file)

# Load kmers from kmers.txt
def load_kmers(file):
    with open(file, 'r') as f:
        return set(line.strip().lower() for line in f)

def calculate_discrepancy_and_test(first_file, second_file, kmers_file):
    # Load the data from both files
    df1 = load_csv(first_file)
    df2 = load_csv(second_file)

    # Load k-mers
    kmers = load_kmers(kmers_file)

    # Merge on "word"
    merged = pd.merge(df1, df2, on='word', suffixes=('_file1', '_file2'))

    # Calculate score discrepancy
    merged['score_discrepancy'] = merged['score_file1'] - merged['score_file2']

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
        'n_special_kmers' : len(kmers_data),
        'median_discrepancy_kmers': median_discrepancy_kmers,
        'median_discrepancy_others': median_discrepancy_others,
        'wilcoxon_statistic': wilcox_stat,
        'p_value': p_value
    }

# Example usage
if __name__ == "__main__":
    args = read_args()
    results = calculate_discrepancy_and_test(args.file1, args.file2, args.kmers)
    filename_base = args.file1.split('_')[0].rsplit("/", 1)[-1] # Need to replace this with more general approach

    results = [filename_base,results['n_special_kmers'], results['median_discrepancy_kmers'], results['median_discrepancy_others'], results['wilcoxon_statistic'],results['p_value']]
    print(','.join([str(x) for x in results]))
