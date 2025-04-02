import pandas as pd
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Average values in genomic data within specified windows")
parser.add_argument("window_size", type=int, help="Size of the window (in base pairs)")
parser.add_argument("input_file", help="Path to the input file")
parser.add_argument("output_file", help="Path to the output file")

args = parser.parse_args()

# Read in the data
df = pd.read_csv(args.input_file, sep="\s+", header=None, names=["chromosome", "position", "value1", "value2", "LR"])

# Sort by chromosome and position
df.sort_values(by=["chromosome", "position"], inplace=True)

# Calculate the window for each position based on the window size
df["window"] = (df["position"] // args.window_size) * args.window_size

# Generate all possible windows per chromosome (even empty ones)
all_windows = []

for chrom in df["chromosome"].unique():
  max_pos = df[df["chromosome"] == chrom]["position"].max()  # Max position in the chromosome
  window_range = range(0, max_pos + args.window_size, args.window_size)
  for window in window_range:
    all_windows.append((chrom, window))

# Create a DataFrame of all possible windows per chromosome
all_windows_df = pd.DataFrame(all_windows, columns=["chromosome", "window"])

# Merge the all_windows_df with the original DataFrame to include empty windows
merged_df = pd.merge(all_windows_df, df, on=["chromosome", "window"], how="left")

# Group by chromosome and window, calculate the average for columns 3-5, and count SNPs
averaged_df = merged_df.groupby(["chromosome", "window"]).agg({
  "value1": "mean",
  "value2": "mean",
  "LR": "mean",
  "position": "count"  # Count the number of SNPs (rows) in each window
}).reset_index()

# Rename the position count column to 'num_snps'
averaged_df.rename(columns={"position": "num_snps"}, inplace=True)

# Explicitly write out NAs
averaged_df = averaged_df.where(pd.notnull(averaged_df), 'NA')

# Output the averaged results to a file
averaged_df.to_csv(args.output_file, sep="\t", header=False, index=False)

print(f"Averaged data with SNP counts has been saved to {args.output_file}")

