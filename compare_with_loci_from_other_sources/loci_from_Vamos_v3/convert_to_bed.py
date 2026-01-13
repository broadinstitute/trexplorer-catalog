"""Converts Vamos TR catalog files to a standardized BED format.

This script reads Vamos catalog TSV files (v2.1 or v3) and converts them to a
BED-formatted file with simplified motifs. The output contains columns:
    chrom, start_0based, end_1based, motif, motif_size

Usage:
    python3 convert_to_bed.py --version 3
"""

import argparse
import os

import pandas as pd

from str_analysis.utils.find_repeat_unit import find_repeat_unit_without_allowing_interruptions


def simplify_motif(motif):
    """Simplifies a motif by finding its minimal repeat unit.

    Args:
        motif: The motif sequence to simplify (e.g., "CAGCAG" -> "CAG").

    Returns:
        The simplified motif string representing the minimal repeat unit.
    """
    simplified_motif, _, _ = find_repeat_unit_without_allowing_interruptions(
        motif, allow_partial_repeats=False)
    return simplified_motif


p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("-v", "--version", choices=["2.1", "3"], default="3",
               help="Vamos catalog version to process (default: 3)")
args = p.parse_args()

if args.version == "2.1":
    input_path = "vamos_catalog.ori.v2.1.tsv.gz"
    df = pd.read_table(input_path)
    output_path = "vamosGenomicTR_v2.1.bed"
elif args.version == "3":
    input_path = "vamosGenomicTR_v3.0_oriMotifs.tsv.gz"
    # The v3 TSV has these columns (many unused):
    # chrom, start_1based, end_1based, motifs, _, STR_or_VNTR, motif_size,
    # consensus_motif, and several other fields we don't need
    df = pd.read_table(input_path, names=[
        "chrom", "start_1based", "end_1based", "motifs", "unused1",
        "STR_or_VNTR", "motif_size", "consensus_motif",
        "unused2", "unused3", "unused4", "unused5", "in_segdup",
        "unused6", "unused7", "unused8"
    ])
    output_path = "vamosGenomicTR_v3.0.bed"
else:
    p.error(f"Invalid version: {args.version}")

print(f"Parsed {input_path}")
print(df.iloc[0])

# Extract the first motif from the comma-separated list and simplify it
df["motif1"] = df["motifs"].str.split(",").str[0]
df["motif1"] = df["motif1"].apply(simplify_motif)
df["motif_size"] = df["motif1"].str.len()

# Convert to 0-based coordinates and select output columns
df["start_0based"] = df["start_1based"] - 1
df = df[["chrom", "start_0based", "end_1based", "motif1", "motif_size"]].copy()
df.sort_values(["chrom", "start_0based", "end_1based"], inplace=True)

print(f"{sum(df['motif_size'] == 1):,d} homopolymers")
print(f"{sum(df['motif_size'] == 2):,d} dinucleotides")

# Write output
df.to_csv(output_path, sep="\t", index=False, header=False)
os.system(f"bgzip -f {output_path}")
os.system(f"tabix -f {output_path}.gz")
print(f"Wrote {len(df):,d} rows to {output_path}.gz")

# Print summary statistics
stats = {
    "STRs (3-6bp)": sum(df['motif_size'] <= 6),
    "VNTRs (7+bp)": sum(df['motif_size'] > 6),
}
for key, count in stats.items():
    print(f"      {count:,d} {key}")
