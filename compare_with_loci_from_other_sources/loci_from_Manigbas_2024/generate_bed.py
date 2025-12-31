#!/usr/bin/env python3
"""
Generate BED file for Manigbas 2024 phenotype-associated TRs.

This script reads Supplementary Data 5 from Manigbas et al. 2024 and creates a BED file
with the 1,634 unique TRs showing causal associations with phenotypes in UK Biobank.

Paper: Manigbas MT, Mori MA, Navarro FC, et al.
A phenome-wide association study of tandem repeat variation in 168,554 individuals from the UK Biobank.
Nat Commun. 2024;15:10509.

Data source: 41467_2024_54678_MOESM3_ESM.xlsx
Sheet: Data S5_Causal analysis

Author: Generated for tandem repeat catalog comparison
"""

import pandas as pd
import os
import sys
sys.path.append('../')
from compare_loci_utils import does_locus_pass_filters


os.chdir(os.path.dirname(__file__))

# Input file
INPUT_FILE = "41467_2024_54678_MOESM3_ESM.xlsx"
SHEET_NAME = "Data S5_Causal analysis"

def main():
    print(f"Reading {INPUT_FILE}, sheet '{SHEET_NAME}'...")

    # Read Data S5, skipping title rows
    df = pd.read_excel(INPUT_FILE, sheet_name=SHEET_NAME, skiprows=2)

    print(f"  Loaded {len(df):,} TR-trait associations from causal analysis")
    print(f"  Columns: {list(df.columns[:10])}")

    # Extract unique TRs (same TR can be associated with multiple traits)
    unique_trs = df[['Chr', 'TR start, hg38', 'TR end, hg38', 'TR motif']].drop_duplicates()

    print(f"\n  Found {len(unique_trs):,} unique TRs")

    # Extract BED format columns
    bed_rows = []

    for idx, row in unique_trs.iterrows():
        chrom = row['Chr']
        start_0based = int(row['TR start, hg38'])
        end_1based = int(row['TR end, hg38'])
        motif = row['TR motif']

        # Handle missing motifs
        if pd.isna(motif):
            continue

        # Ensure chromosome has 'chr' prefix (already does)
        if not str(chrom).startswith('chr'):
            chrom = f'chr{chrom}'

        passes_filters, adjusted_motif = does_locus_pass_filters(
            chrom, start_0based, end_1based, motif,
            min_repeats_in_reference=2,
            min_adjusted_motif_purity=0.2,
        )

        if passes_filters:
            if len(adjusted_motif) != len(motif):
                print(f"Simplified {chrom}:{start_0based}-{end_1based} motif: {motif} => {adjusted_motif}")
        
            bed_rows.append((chrom, start_0based, end_1based, adjusted_motif))

    print(f"\n  Processed {len(bed_rows):,} phenotype-associated TR loci")

    # Sort by chromosome and position
    chrom_order = {f'chr{i}': i for i in range(1, 23)}
    chrom_order.update({'chrX': 23, 'chrY': 24, 'chrM': 25})

    def sort_key(row):
        chrom, start, end, motif = row
        return (chrom_order.get(chrom, 99), start)

    bed_rows.sort(key=sort_key)

    # Write BED file
    output_file = 'Manigbas_2024_phenotype_associated_TRs.bed'
    print(f"\nWriting {output_file}...")

    with open(output_file, 'w') as f:
        for chrom, start, end, motif in bed_rows:
            f.write(f'{chrom}\t{start}\t{end}\t{motif}\n')

    # Compress and index
    print(f"Compressing with bgzip...")
    os.system(f'bgzip -f {output_file}')
    os.system(f'tabix -p bed {output_file}.gz')

    print(f"\nDone! Created {output_file}.gz with {len(bed_rows):,} loci")
    print(f"First few entries:")
    for row in bed_rows[:5]:
        print(f"  {row[0]}\t{row[1]}\t{row[2]}\t{row[3]}")

if __name__ == '__main__':
    main()
