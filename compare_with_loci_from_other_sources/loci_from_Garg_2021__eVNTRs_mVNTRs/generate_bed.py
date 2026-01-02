#!/usr/bin/env python3
"""
Extract TR loci from Garg 2021 supplementary tables (eVNTRs and mVNTRs)

Source:
- mmc5.xlsx (Table S4): eVNTRs - VNTRs associated with gene expression
- mmc6.xlsx (Table S5): mVNTRs - VNTRs associated with DNA methylation

Paper: Garg et al. 2021. "Pervasive cis effects of variation in copy number of
large tandem repeats on local DNA methylation and gene expression"
AJHG. https://doi.org/10.1016/j.ajhg.2021.03.012
"""

from collections import Counter
import os
import pandas as pd
import re
import sys
sys.path.append('../')
from compare_loci_utils import does_locus_pass_filters


os.chdir(os.path.dirname(__file__))


def parse_coordinate(coord_str):
    """Parse coordinate string like 'chr1:100238124-100238280' to (chrom, start, end)"""
    if pd.isna(coord_str):
        return None
    match = re.match(r'(chr[^:]+):(\d+)-(\d+)', str(coord_str))
    if match:
        return match.group(1), int(match.group(2)), int(match.group(3))
    return None


def parse_motif(motif_str):
    """Parse motif string, handling multiple motifs separated by commas"""
    if pd.isna(motif_str):
        return None
    # Take the first motif if multiple are listed
    motif = str(motif_str).split(',')[0].strip()
    # Remove any non-ACGT characters
    motif = re.sub(r'[^ACGTacgt]', '', motif)
    return motif.upper() if motif else None


# Read eVNTRs from mmc5.xlsx
print("Reading eVNTRs from mmc5.xlsx...")
df_evntr = pd.read_excel("mmc5.xlsx", header=1)

# The actual column names are in the first row of data, so we need header=1
# But the columns are still unnamed, so let's get the first data row as headers
evntr_cols = df_evntr.iloc[0].tolist()
df_evntr = df_evntr.iloc[1:].reset_index(drop=True)
df_evntr.columns = evntr_cols

print(f"  Read {len(df_evntr)} eVNTR associations")

# Read mVNTRs from mmc6.xlsx
print("Reading mVNTRs from mmc6.xlsx...")
df_mvntr = pd.read_excel("mmc6.xlsx", header=1)

mvntr_cols = df_mvntr.iloc[0].tolist()
df_mvntr = df_mvntr.iloc[1:].reset_index(drop=True)
df_mvntr.columns = mvntr_cols

print(f"  Read {len(df_mvntr)} mVNTR associations")

# Extract unique VNTR loci from eVNTRs
vntr_loci = {}  # key: (chrom, start, end), value: motif

skipped_evntr = 0
for idx, row in df_evntr.iterrows():
    coord = parse_coordinate(row.get('VNTR coordinate (hg38)'))
    if coord is None:
        skipped_evntr += 1
        continue

    chrom, start, end = coord
    motif = parse_motif(row.get('VNTR motif sequence'))

    if motif is None or len(motif) == 0:
        skipped_evntr += 1
        continue

    if len(motif) > 1000:  # Skip extremely long motifs
        skipped_evntr += 1
        continue

    key = (chrom, start, end)
    if key not in vntr_loci:
        vntr_loci[key] = motif

print(f"  Extracted {len(vntr_loci)} unique eVNTR loci (skipped {skipped_evntr})")

# Extract unique VNTR loci from mVNTRs
skipped_mvntr = 0
mvntr_count = 0
for idx, row in df_mvntr.iterrows():
    coord = parse_coordinate(row.get('VNTR coordinate (hg38)'))
    if coord is None:
        skipped_mvntr += 1
        continue

    chrom, start, end = coord
    motif = parse_motif(row.get('Motif sequence'))

    if motif is None or len(motif) == 0:
        skipped_mvntr += 1
        continue

    if len(motif) > 1000:
        skipped_mvntr += 1
        continue

    key = (chrom, start, end)
    if key not in vntr_loci:
        vntr_loci[key] = motif
        mvntr_count += 1

print(f"  Added {mvntr_count} unique mVNTR loci not in eVNTRs (skipped {skipped_mvntr})")

# Sort by genomic coordinates
output_rows = []
filter_reason_counts = Counter()
for (chrom, start_0based, end_1based), motif in vntr_loci.items():
    passes_filters, adjusted_motif, filter_reason = does_locus_pass_filters(
        chrom, start_0based, end_1based, motif,
        min_repeats_in_reference=2,
        min_adjusted_motif_purity=0.2,
    )

    if passes_filters:
        if len(adjusted_motif) != len(motif):
            print(f"Simplified {chrom}:{start_0based}-{end_1based} motif: {motif} => {adjusted_motif}")

        output_rows.append((chrom, start_0based, end_1based, adjusted_motif))
    else:
        filter_reason_counts[filter_reason] += 1

print(f"Kept {len(output_rows):,d} out of {len(vntr_loci):,d} loci after applying filters")
if filter_reason_counts:
    print("Filter reasons:")
    for reason, count in filter_reason_counts.most_common():
        print(f"  {count:10,d}  {reason}")
output_rows.sort(key=lambda x: (x[0].replace('chr', '').zfill(2), x[1], x[2]))

# Write output BED file
output_bed_path = "Garg_2021_eVNTRs_mVNTRs.bed"
with open(output_bed_path, "wt") as output_bed:
    for chrom, start, end, motif in output_rows:
        output_bed.write(f"{chrom}\t{start}\t{end}\t{motif}\n")

# Compress and index
os.system(f"bgzip -f {output_bed_path}")
os.system(f"tabix -f {output_bed_path}.gz")

print("=" * 50)
print(f"Wrote {len(output_rows):,d} loci to {output_bed_path}.gz")
print("=" * 50)
