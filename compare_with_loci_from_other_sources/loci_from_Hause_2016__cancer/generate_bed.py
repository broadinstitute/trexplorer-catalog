#!/usr/bin/env python3
"""
Generate BED file from Hause 2016 cancer MSI loci.

This script processes Supplementary Table 10 from Hause et al. 2016,
which contains 2,685 cancer-specific microsatellite instability loci.

Paper: Hause RJ, Pritchard CC, Shendure J, et al.
Classification and characterization of microsatellite instability across 18 cancer types.
Nat Med. 2016;22(11):1342-1350.

Data source: 41591_2016_BFnm4191_MOESM30_ESM.xlsx
Sheet: msih_cancer_specific_instabilit

Coordinates are in hg19 and will be converted to hg38.

Author: Generated for tandem repeat catalog comparison
"""

from collections import Counter
from liftover import get_lifter
import os
import pandas as pd
import re
import sys
sys.path.append('../')
from compare_loci_utils import does_locus_pass_filters

os.chdir(os.path.dirname(__file__))

INPUT_FILE = "41591_2016_BFnm4191_MOESM30_ESM.xlsx"

def parse_locus(locus_str):
    """
    Parse locus string like '2:148683681-148683698' to extract chr, start, end.
    Returns (chrom, start_1based, end_1based)
    """
    match = re.match(r'(\d+|X|Y):(\d+)-(\d+)', locus_str)
    if not match:
        return None, None, None

    chrom, start, end = match.groups()
    return f'chr{chrom}', int(start), int(end)

def parse_repeat_sequence(seq_str):
    """
    Parse repeat notation like '(A)8' to extract the motif.
    Returns None for complex repeats (multiple motifs) or if parsing fails.
    """
    if pd.isna(seq_str):
        return None

    matches = re.findall(r'\(([ACGT]+)\)', str(seq_str))

    if len(matches) != 1:
        return None

    return matches[0]

def main():
    print(f"Reading {INPUT_FILE}...")

    df = pd.read_excel(INPUT_FILE, sheet_name=0)

    print(f"  Loaded {len(df):,} MSI loci")
    print(f"  Columns: {list(df.columns)}")

    # Initialize liftOver from hg19 to hg38
    print("\nInitializing liftOver from hg19 to hg38...")
    converter = get_lifter('hg19', 'hg38', one_based=False)

    output_rows = []
    skipped_locus_parse = 0
    skipped_complex_motif = 0
    skipped_liftover = 0

    for idx, row in df.iterrows():
        # Parse locus string
        chrom_hg19, start_hg19_1based, end_hg19_1based = parse_locus(row['locus'])

        if chrom_hg19 is None:
            skipped_locus_parse += 1
            continue

        # Parse repeat sequence to get motif
        motif = parse_repeat_sequence(row['repeat_dna_sequence'])

        if motif is None:
            skipped_complex_motif += 1
            continue

        # Convert to 0-based for liftOver
        start_hg19_0based = start_hg19_1based - 1
        end_hg19_0based = end_hg19_1based - 1

        # Perform liftOver for start position
        liftover_result_start = converter[chrom_hg19][start_hg19_0based]

        if len(liftover_result_start) != 1 or len(liftover_result_start[0]) != 3:
            skipped_liftover += 1
            continue

        chrom_hg38, start_hg38_0based, strand_hg38 = liftover_result_start[0]

        # LiftOver end position
        liftover_result_end = converter[chrom_hg19][end_hg19_0based]
        if len(liftover_result_end) != 1:
            skipped_liftover += 1
            continue

        _, end_hg38_0based, _ = liftover_result_end[0]
        end_hg38_1based = end_hg38_0based + 1

        # Handle cases where liftOver reverses the coordinates
        if end_hg38_1based < start_hg38_0based:
            start_hg38_0based, end_hg38_1based = end_hg38_0based, start_hg38_0based + 1

        passes_filters, adjusted_motif = does_locus_pass_filters(
            chrom_hg38, start_hg38_0based, end_hg38_1based, motif,
            min_repeats_in_reference=2,
            min_adjusted_motif_purity=0.2,
        )

        if passes_filters:
            if len(adjusted_motif) != len(motif):
                print(f"Simplified {chrom_hg38}:{start_hg38_0based}-{end_hg38_1based} motif: {motif} => {adjusted_motif}")

            output_rows.append((chrom_hg38, start_hg38_0based, end_hg38_1based, adjusted_motif))

    print(f"\nProcessed {len(output_rows):,} loci successfully")
    total_skipped = skipped_locus_parse + skipped_complex_motif + skipped_liftover
    print(f"Skipped {total_skipped:,} loci:")
    if skipped_locus_parse:
        print(f"  - {skipped_locus_parse:,} due to locus parsing errors")
    if skipped_complex_motif:
        print(f"  - {skipped_complex_motif:,} due to complex/compound motifs")
    if skipped_liftover:
        print(f"  - {skipped_liftover:,} due to liftOver failures")

    # Print motif size distribution
    motif_size_counts = Counter(len(row[3]) for row in output_rows)
    print(f"\nMotif size distribution:")
    for size in sorted(motif_size_counts.keys()):
        count = motif_size_counts[size]
        label = {1: "homopolymers", 2: "dinucleotides", 3: "trinucleotides",
                 4: "tetranucleotides", 5: "pentanucleotides", 6: "hexanucleotides"}.get(size, f"{size}bp")
        print(f"  {size}bp ({label}): {count:,} ({100*count/len(output_rows):.1f}%)")

    # Sort by chromosome and position
    chrom_order = {f'chr{i}': i for i in range(1, 23)}
    chrom_order.update({'chrX': 23, 'chrY': 24, 'chrM': 25})

    def sort_key(row):
        chrom, start, end, motif = row
        return (chrom_order.get(chrom, 99), start)

    output_rows.sort(key=sort_key)

    # Write output BED file
    output_file = "Hause_2016_cancer_MSI_loci.bed"
    print(f"\nWriting {output_file}...")

    with open(output_file, "w") as f:
        for chrom, start, end, motif in output_rows:
            f.write(f"{chrom}\t{start}\t{end}\t{motif}\n")

    # Compress and index
    print(f"Compressing with bgzip...")
    os.system(f"bgzip -f {output_file}")
    os.system(f"tabix -p bed {output_file}.gz")

    print(f"\nDone! Created {output_file}.gz with {len(output_rows):,} loci")
    print(f"First few entries:")
    for row in output_rows[:5]:
        print(f"  {row[0]}\t{row[1]}\t{row[2]}\t{row[3]}")

if __name__ == '__main__':
    main()
