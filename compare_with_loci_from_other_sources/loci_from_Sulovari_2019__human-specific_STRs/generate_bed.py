#!/usr/bin/env python3
"""
Generate BED file for Sulovari 2019 human-specific expanded TRs.

This script reads Dataset S9 from Sulovari et al. 2019 and creates a BED file
with the 1,584 human-specific tandem repeat expansions.

Paper: Sulovari A, Li R, Audano PA, et al.
Human-specific tandem repeat expansion and differential gene expression during primate evolution.
Proc Natl Acad Sci U S A. 2019;116(46):23243-23253.

Data source: pnas.1912175116.sd09.xlsx
Dataset S9: Full list of ab initio and HSE loci and their tandem repeat
annotation in human and NHP haplotypes.

Author: Generated for tandem repeat catalog comparison
"""

import pandas as pd
import os
import sys
sys.path.append('../')
from compare_loci_utils import does_locus_pass_filters


os.chdir(os.path.dirname(__file__))

# Input file
INPUT_FILE = "pnas.1912175116.sd09.xlsx"

def main():
    print(f"Reading {INPUT_FILE}...")

    # Read Dataset S9, skipping the title row
    df = pd.read_excel(INPUT_FILE, skiprows=1)

    print(f"  Loaded {len(df):,} human-specific expanded TR loci")
    print(f"  Columns: {list(df.columns[:15])}...")

    # Extract relevant columns
    # CHR, START, END are the genomic coordinates (GRCh38/hg38)
    # GRCh38_motif is the repeat motif
    bed_rows = []

    for idx, row in df.iterrows():
        chrom = row['CHR']
        start = int(row['START'])
        end = int(row['END'])
        motif = row['GRCh38_motif']

        # Handle missing motifs
        if pd.isna(motif):
            continue

        # Ensure chromosome has 'chr' prefix
        if not str(chrom).startswith('chr'):
            chrom = f'chr{chrom}'

        passes_filters, adjusted_motif = does_locus_pass_filters(
            chrom, start, end, motif,
            min_repeats_in_reference=2,
            min_adjusted_motif_purity=0.2,
        )

        if passes_filters:
            if len(adjusted_motif) != len(motif):
                print(f"Simplified {chrom}:{start}-{end} motif: {motif} => {adjusted_motif}")
            bed_rows.append((chrom, start, end, adjusted_motif))


    print(f"\n  {len(bed_rows):,} out of {len(df):,} loci passed filters")

    # Sort by chromosome and position
    chrom_order = {f'chr{i}': i for i in range(1, 23)}
    chrom_order.update({'chrX': 23, 'chrY': 24, 'chrM': 25})

    def sort_key(row):
        chrom, start, end, motif = row
        return (chrom_order.get(chrom, 99), start)

    bed_rows.sort(key=sort_key)

    # Write BED file
    output_file = 'Sulovari_2019_human_specific_STRs.bed'
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
