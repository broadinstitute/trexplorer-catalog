"""
Extract tandem repeat loci from Danzi et al. 2025

Paper: "Detailed tandem repeat allele profiling in 1,027 long-read genomes reveals genome-wide patterns of pathogenicity"
URL: https://www.biorxiv.org/content/10.1101/2025.01.06.631535v1.full

Data source:
- Supplementary Table 3 (media-3.xlsx): Contains 923,089 TR loci with constraint metrics
- longestPureSegmentQuantiles.txt.gz: Additional file from authors containing motif sequences

Selection criteria (based on Figure 5C):
- OE_len <= 0.33 (constrained loci with observed/expected length ratio <= 0.33)
- OE_len >= 2 (expanded loci with observed/expected length ratio >= 2)
- combinedLPSStdev_percentile >= 0.8 (80th percentile of LPS standard deviation)
- combinedLPSStdev >= 0.2 (minimum LPS standard deviation)

Coordinate system: hg38/GRCh38 (no liftover needed)

Number of loci extracted: ~134,390 outlier loci from 923,089 total loci

Processing notes:
- Motifs are taken from longestPureSegmentQuantiles.txt.gz
- Multiple output files are generated with different boundary definitions:
  - wide_boundaries: Original boundaries from Danzi et al.
  - narrow_adotto_boundaries: Boundaries from overlapping loci in Adotto catalog v1.2
- Filtering variants based on motif matching with Adotto catalog
"""

import os
import pandas as pd

os.chdir(os.path.dirname(__file__))

# Parse Danzi et al. 2025 Supp. Table 3 and filter to outliers
input_supp_table3_path = "media-3.xlsx"
print(f"Reading {input_supp_table3_path}")
df = pd.read_excel(input_supp_table3_path)
print(f"Read {len(df):,d} intervals from {input_supp_table3_path}")

# Extract chromosome, start, end from TRID column
df[["chrom", "start_0based", "end_1based"]] = df["TRID"].str.split("_", expand=True)
df["start_0based"] = df["start_0based"].astype(int)
df["end_1based"] = df["end_1based"].astype(int)

# Compute percentile of combinedLPSStdev
df["combinedLPSStdev_percentile"] = df["combinedLPSStdev"].rank(pct=True)

# Filter based on thresholds from Figure 5C
OE_len_threshold1 = 0.33
OE_len_threshold2 = 2
combinedLPSStdev_percentile_min_threshold = 0.8
combinedLPSStdev_min_threshold = 0.2

total = len(df)
df = df[(df["OE_len"] <= OE_len_threshold1)
    | (df["OE_len"] >= OE_len_threshold2)
    | (df["combinedLPSStdev_percentile"] >= combinedLPSStdev_percentile_min_threshold)
    | (df["combinedLPSStdev"] >= combinedLPSStdev_min_threshold)
].copy()

print(f"Kept {len(df):,d} out of {total:,d} ({(len(df) / total):.1%}) intervals after filtering by outlier criteria")

# Load motifs from additional table
df2_path = "longestPureSegmentQuantiles.txt.gz"
print(f"Reading {df2_path}")
df2 = pd.read_table(df2_path)
print(f"Read {len(df2):,d} rows from {df2_path}")

# Keep only the first motif for each TRID (sorted by N_motif descending)
df2.sort_values(by=["N_motif"], inplace=True, ascending=False)
df2 = df2.groupby("TRID").first().reset_index()

# Create lookup dictionary
trid_to_motif_lookup = {
    trid: motif for trid, motif in zip(df2["TRID"], df2["longestPureSegmentMotif"])
}

# Add motifs to dataframe
df["motif"] = df["TRID"].map(trid_to_motif_lookup)

# Check for missing motifs
if sum(df['motif'].isna()) > 0:
    print(f"WARNING: {sum(df['motif'].isna()):,d} out of {len(df):,d} loci are missing motif information")
    df = df[df['motif'].notna()].copy()

# Write output BED file
output_bed_path = "Danzi_2025_loci.bed"
with open(output_bed_path, "wt") as f:
    for _, row in sorted(df.iterrows(), key=lambda x: (x[1]['chrom'], x[1]['start_0based'], x[1]['end_1based'])):
        f.write(f"{row['chrom']}\t{row['start_0based']}\t{row['end_1based']}\t{row['motif']}\n")

os.system(f"bgzip -f {output_bed_path}")
os.system(f"tabix -f {output_bed_path}.gz")
print(f"Wrote {len(df):,d} loci to {output_bed_path}.gz")
