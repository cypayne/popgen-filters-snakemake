#!/usr/bin/env python3

## thin snps so that no two SNPs are within N bp of eachother
## usage: python3 ~/scripts/thin_snps_by_distance_bp.py in.bcf out.bcf -N 1000 
## cyp+chatgpt 06-2025

import argparse
import pysam

def filter_by_distance(infile, outfile, min_dist):
  bcf_in = pysam.VariantFile(infile, 'r')
  bcf_out = pysam.VariantFile(outfile, 'w', header=bcf_in.header)
  last_pos = {}  # track last included position per contig
  for rec in bcf_in.fetch():
    if rec.alts is None:
      continue
    # Only SNPs: one-base alt, single nucleotide
    if len(rec.ref) == 1 and all(len(alt) == 1 for alt in rec.alts):
      chrom = rec.contig
      pos = rec.pos  # 1-based
      prev = last_pos.get(chrom, -10**12)
      if pos - prev >= min_dist:
        bcf_out.write(rec)
        last_pos[chrom] = pos
    else:
      # You can decide: include or skip indels/MNPs
      bcf_out.write(rec)

  bcf_in.close()
  bcf_out.close()

def main():
  p = argparse.ArgumentParser(description="Filter SNPs by minimum distance.")
  p.add_argument("input", help="Input BCF/VCF file")
  p.add_argument("output", help="Output BCF/VCF file")
  p.add_argument("-N", "--min-dist", type=int, default=1,
                 help="Minimum distance between SNPs (default: 1 bp)")
  args = p.parse_args()

  filter_by_distance(args.input, args.output, args.min_dist)

if __name__ == "__main__":
  main()
