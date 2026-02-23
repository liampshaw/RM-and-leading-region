import re
import glob
import os
import numpy as np
import pandas as pd
from Bio import SeqIO

BIN_SIZE = 5000


# ---------- Parsing ----------

def parse_einverted_to_df(file_path):
    """
    Parse einverted file into DataFrame:
    columns = [contig, start, end]
    """
    contigs = []
    starts = []
    ends = []

    with open(file_path) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]

        if ":" in line:
            contig = line.split(":")[0]

            seq1 = lines[i+1]
            seq2 = lines[i+3]

            m1 = re.search(r"\s+(\d+)\s+[acgtACGT]+\s+(\d+)", seq1)
            m2 = re.search(r"\s+(\d+)\s+[acgtACGT]+\s+(\d+)", seq2)

            if m1 and m2:
                s1, e1 = int(m1.group(1)), int(m1.group(2))
                s2, e2 = int(m2.group(1)), int(m2.group(2))

                start = min(s1, e1, s2, e2)
                end = max(s1, e1, s2, e2)

                contigs.append(contig)
                starts.append(start)
                ends.append(end)

            i += 4
        else:
            i += 1

    df = pd.DataFrame({
        "contig": contigs,
        "start": np.array(starts, dtype=np.int32),
        "end": np.array(ends, dtype=np.int32)
    })

    return df


def load_lengths(fasta_path):
    return {
        rec.id: len(rec.seq)
        for rec in SeqIO.parse(fasta_path, "fasta")
    }


# ---------- Core vectorised computation ----------

def compute_density_df(hit_df, lengths, genome_name):
    """
    Vectorised binning + counting
    """

    # midpoint (vectorised)
    midpoints = ((hit_df["start"].values + hit_df["end"].values) // 2).astype(np.int32)

    hit_df = hit_df.copy()
    hit_df["mid"] = midpoints
    hit_df["bin"] = midpoints // BIN_SIZE
    hit_df["genome"] = genome_name

    # filter contigs not in fasta
    hit_df = hit_df[hit_df["contig"].isin(lengths.keys())]

    # count hits per bin
    grouped = (
        hit_df
        .groupby(["genome", "contig", "bin"], sort=False)
        .size()
        .rename("count")
        .reset_index()
    )

    # ---------- Build full bin table ----------
    records = []

    for contig, length in lengths.items():
        n_bins = (length // BIN_SIZE) + 1
        bins = np.arange(n_bins, dtype=np.int32)

        starts = bins * BIN_SIZE
        ends = np.minimum((bins + 1) * BIN_SIZE, length)

        df = pd.DataFrame({
            "genome": genome_name,
            "contig": contig,
            "bin": bins,
            "start": starts,
            "end": ends
        })

        records.append(df)

    bins_df = pd.concat(records, ignore_index=True)

    # ---------- Merge counts ----------
    merged = bins_df.merge(grouped, on=["genome", "contig", "bin"], how="left")
    merged["count"] = merged["count"].fillna(0).astype(np.int32)

    # density per bp
    merged["bin_len"] = merged["end"] - merged["start"]
    merged["density_per_bp"] = merged["count"] / merged["bin_len"]

    return merged[["genome", "contig", "start", "end", "count", "density_per_bp"]]


# ---------- Batch processing ----------

def match_files(inv_pattern, fasta_pattern):
    inv_files = sorted(glob.glob(inv_pattern))
    fasta_files = sorted(glob.glob(fasta_pattern))

    # match by basename prefix
    fasta_dict = {
        os.path.splitext(os.path.basename(f))[0]: f
        for f in fasta_files
    }

    pairs = []
    for inv in inv_files:
        base = os.path.splitext(os.path.basename(inv))[0]
        if base in fasta_dict:
            pairs.append((inv, fasta_dict[base]))
        else:
            print(f"WARNING: No FASTA match for {inv}")

    return pairs


def run(inv_pattern, fasta_pattern, output_tsv):
    pairs = match_files(inv_pattern, fasta_pattern)

    all_results = []

    for inv_file, fasta_file in pairs:
        genome_name = os.path.basename(inv_file)

        hit_df = parse_einverted_to_df(inv_file)
        lengths = load_lengths(fasta_file)

        result_df = compute_density_df(hit_df, lengths, genome_name)
        all_results.append(result_df)

    final_df = pd.concat(all_results, ignore_index=True)

    final_df.to_csv(output_tsv, sep="\t", index=False)


# ---------- CLI ----------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--inv", required=True, help="Pattern e.g. '*.inv'")
    parser.add_argument("--fasta", required=True, help="Pattern e.g. '*.fa'")
    parser.add_argument("-o", required=True, help="Output TSV")

    args = parser.parse_args()

    run(args.inv, args.fasta, args.o)
