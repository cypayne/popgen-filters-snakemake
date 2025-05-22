#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(description="Calculate LD decay curve from pairwise R² data.")
    parser.add_argument("input", help="LD input file (.gz or .txt)")
    parser.add_argument("output", help="Output file (TSV with binned average R²)")
    parser.add_argument("--window", type=int, default=1000000, help="Window size in base pairs")
    parser.add_argument("--bins", type=int, default=100, help="Number of distance bins")
    args = parser.parse_args()

    print("Reading file...")
    df = pd.read_csv(args.input, sep="\t", compression='infer',
            usecols=["CHR_A", "BP_A", "BP_B", "R2"], dtype=str)

    # Remove repeated header lines
    df = df[df["CHR_A"] != "CHR_A"]

    # Cast columns to numeric
    df["BP_A"] = pd.to_numeric(df["BP_A"])
    df["BP_B"] = pd.to_numeric(df["BP_B"])
    df["R2"] = pd.to_numeric(df["R2"])

    print("Processing chromosomes...")
    results = []

    for chrom in df["CHR_A"].unique():
        sub = df[df["CHR_A"] == chrom].copy()
        sub["Dist"] = sub["BP_B"] - sub["BP_A"]
        sub = sub[sub["Dist"] > 0]  # Remove any non-positive distances

    # Bin into quantiles
    sub["bin"] = pd.qcut(sub["Dist"], q=args.bins, labels=False, duplicates='drop') + 1

    # Aggregate
    grouped = sub.groupby("bin").agg(
        num=("R2", "count"),
        avg_R2=("R2", "mean"),
        std=("R2", "std")
    ).reset_index()

    # Convert bin index to midpoint distances
    grouped["distance"] = grouped["bin"] * (args.window / args.bins)
    grouped.insert(0, "chr", chrom)
    results.append(grouped[["chr", "distance", "num", "avg_R2", "std"]])

    # Combine all chromosomes
    final_df = pd.concat(results, ignore_index=True)
    final_df.to_csv(args.output, sep="\t", index=False)
    print(f"Done! Output written to: {args.output}")

if __name__ == "__main__":
    main()

