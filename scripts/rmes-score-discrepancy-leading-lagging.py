import pandas as pd
import scipy.stats as stats
import numpy as np
import argparse 
import re
import glob
import os

def read_args():
    parser = argparse.ArgumentParser(description="Merge results for .")
    parser.add_argument("--directory", required=True, help="Search path for output file.")
    parser.add_argument("--kmers", required=True, help="Path to the kmers text file.")
    parser.add_argument("--k", required=True, help="Value of k.")
    parser.add_argument("--output", required=True, help="Output file prefix.")

    return parser.parse_args()

# Load data from a csv file
def load_csv(file):
    return pd.read_csv(file)

# Load kmers from kmers.txt
def load_kmers(file):
    with open(file, 'r') as f:
        return set(line.strip().lower() for line in f)

def merge_score_files(first_file, second_file, kmer_file):
    # Load the data from both files
    df1 = load_csv(first_file)
    df2 = load_csv(second_file)
    kmers = load_kmers(kmer_file)

    # Merge on "word"
    merged = pd.merge(df1, df2, on='word', suffixes=('_file1', '_file2'))

    # Calculate score discrepancy
    merged['score_discrepancy'] = merged['score_file1'] - merged['score_file2']
    merged['count_discrepancy'] = merged['count_file1'] - merged['count_file2']

    merged = merged[merged['word'].isin(kmers)]
    kmers_used = merged['word']
    merged = merged[['word', 'score_discrepancy', 'count_discrepancy']]
    sorted_df = merged.sort_values(by='word')
    return [list(sorted_df['score_discrepancy']) , list(sorted_df['count_discrepancy'])]

# For each plasmid, compute the discrepancies and then return the results 
# save them in a dictionary with motifs as keys
def merge_comparisons(directory, kmers_file, k):
    # Find all leading and lagging files in the given directory
    leading_files = sorted(glob.glob(os.path.join(directory, "*leading*"+str(k)+".csv")))
    lagging_files = sorted(glob.glob(os.path.join(directory, "*lagging*"+str(k)+".csv")))
    kmers = load_kmers(kmers_file)
    all_scores = []
    all_counts = []
    for leading_file, lagging_file in zip(leading_files, lagging_files):
        results = merge_score_files(leading_file, lagging_file, kmers_file)
        all_scores.append(results[0])
        all_counts.append(results[1])
    return [all_scores, all_counts]



# Example usage
if __name__ == "__main__":
    args = read_args()
    kmers = sorted(load_kmers(args.kmers))
    results = merge_comparisons(args.directory, args.kmers, args.k)
    # Only use k-mers of right length
    kmers_reduced = [k for k in kmers if len(k)==int(args.k)]

    scores = results[0]
    counts = results[1]
    # Convert scores list of lists into a DataFrame
    scores_df = pd.DataFrame(scores).T  # Transpose to align scores with kmers
    # Assign column names
    scores_df.columns = [f'score{i+1}' for i in range(len(scores))]
    scores_df.insert(0, 'kmer', kmers_reduced)  # Insert kmers as the first column
    scores_df.to_csv(args.output+'_scores.csv', index=False)

    counts_df = pd.DataFrame(counts).T  # Transpose to align scores with kmers
    # Assign column names
    counts_df.columns = [f'count{i+1}' for i in range(len(scores))]
    counts_df.insert(0, 'kmer', kmers_reduced)  # Insert kmers as the first column
    counts_df.to_csv(args.output+'_counts.csv', index=False)


#    results.to_csv(args.output, index=False)
