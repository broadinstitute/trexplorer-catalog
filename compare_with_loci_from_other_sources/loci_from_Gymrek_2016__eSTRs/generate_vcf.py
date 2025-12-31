#!/usr/bin/env python3
"""
Generate VCF file for Gymrek 2016 eSTRs.

This script reads the Gymrek 2016 supplementary CSV (hg19 coordinates), performs
liftOver to hg38, and creates a VCF file where each STR is represented as an
insertion of the motif at the STR start position.

Data source: NIHMS784396-supplement-Supplemental_File_TR_eQTLs.csv

Paper: Gymrek M, Willems T, Guilmatre A, et al.
Abundant contribution of short tandem repeats to gene expression variation in humans.
Nat Genet. 2016;48(1):22-29.
"""

from liftover import get_lifter
import os
import pandas as pd
import pysam

os.chdir(os.path.dirname(__file__))

# Input files
CSV_FILE = "NIHMS784396-supplement-Supplemental_File_TR_eQTLs.csv"
REFERENCE_FASTA = os.path.expanduser("~/hg38.fa")

# Output file
OUTPUT_VCF = "Gymrek_2016_eSTRs.vcf"


def main():
    # Initialize liftOver from hg19 to hg38
    print("Initializing liftOver (hg19 -> hg38)...")
    converter = get_lifter('hg19', 'hg38', one_based=False)

    # Load reference genome
    print(f"Loading reference genome: {REFERENCE_FASTA}")
    fasta = pysam.FastaFile(REFERENCE_FASTA)

    # Read eSTR data
    print(f"Reading eSTR data from {CSV_FILE}...")
    df = pd.read_csv(CSV_FILE)
    print(f"  Loaded {len(df):,} eSTR associations")

    # Filter to significant eSTRs only (gene-level FDR 5%)
    df = df[df['signif.estr'] == True]
    print(f"  Filtered to {len(df):,} significant eSTR associations")

    # Get unique STR loci (same locus may be associated with multiple genes)
    unique_strs = df[['chrom', 'str.start', 'motif']].drop_duplicates()
    print(f"  Found {len(unique_strs):,} unique significant eSTR loci")

    # Generate VCF entries
    vcf_records = []
    skipped_liftover = 0
    skipped_ref = 0

    for i, (_, row) in enumerate(unique_strs.iterrows()):
        if i % 10000 == 0 and i > 0:
            print(f"  Processed {i:,}/{len(unique_strs):,} loci...")

        chrom_hg19 = row['chrom']
        pos_1based_hg19 = int(row['str.start'])
        motif = str(row['motif'])

        if pd.isna(motif) or motif == 'nan':
            skipped_ref += 1
            continue

        # Convert to 0-based position for liftOver
        pos_0based_hg19 = pos_1based_hg19 - 1

        # Perform liftOver
        liftover_result = converter[chrom_hg19][pos_0based_hg19]

        if len(liftover_result) != 1 or len(liftover_result[0]) != 3:
            skipped_liftover += 1
            continue

        chrom_hg38, pos_0based_hg38, strand_hg38 = liftover_result[0]

        try:
            # Get reference base at the position before the STR (for VCF anchor base)
            # The anchor base is at pos_0based - 1 (the base immediately before the STR)
            anchor_pos_0based = pos_0based_hg38 - 1
            if anchor_pos_0based < 0:
                skipped_ref += 1
                continue

            ref_base = fasta.fetch(chrom_hg38, anchor_pos_0based, anchor_pos_0based + 1).upper()

            if not ref_base or ref_base == 'N':
                skipped_ref += 1
                continue

            # Create ALT allele as ref_base + motif (representing an insertion)
            alt_allele = ref_base + motif

            # VCF uses 1-based positions - POS is the anchor base position
            pos_1based_hg38 = anchor_pos_0based + 1

            vcf_record = {
                'CHROM': chrom_hg38,
                'POS': pos_1based_hg38,
                'ID': '.',
                'REF': ref_base,
                'ALT': alt_allele,
                'QUAL': '.',
                'FILTER': 'PASS',
                'INFO': f'MOTIF={motif};SOURCE=Gymrek2016'
            }

            vcf_records.append(vcf_record)

        except Exception as e:
            skipped_ref += 1
            continue

    print(f"\nSkipped {skipped_liftover:,} loci (liftOver failed)")
    print(f"Skipped {skipped_ref:,} loci (no motif or invalid reference)")

    # Create DataFrame and sort
    vcf_df = pd.DataFrame(vcf_records)

    # Sort by chromosome and position
    chrom_order = {f'chr{i}': i for i in range(1, 23)}
    chrom_order.update({'chrX': 23, 'chrY': 24, 'chrM': 25})

    vcf_df['chrom_num'] = vcf_df['CHROM'].map(chrom_order)
    vcf_df = vcf_df.dropna(subset=['chrom_num'])
    vcf_df = vcf_df.sort_values(['chrom_num', 'POS']).drop('chrom_num', axis=1)

    # Write VCF file
    print(f"\nWriting VCF file: {OUTPUT_VCF}")
    with open(OUTPUT_VCF, 'w') as f:
        # Write VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##reference={REFERENCE_FASTA}\n")
        f.write("##INFO=<ID=MOTIF,Number=1,Type=String,Description=\"STR repeat motif\">\n")
        f.write("##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Data source\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for _, row in vcf_df.iterrows():
            f.write(f"{row['CHROM']}\t{int(row['POS'])}\t{row['ID']}\t{row['REF']}\t{row['ALT']}\t{row['QUAL']}\t{row['FILTER']}\t{row['INFO']}\n")

    # Compress and index
    print("Compressing with bgzip...")
    os.system(f"bgzip -f {OUTPUT_VCF}")
    os.system(f"tabix -p vcf {OUTPUT_VCF}.gz")

    print(f"\nDone! Created {OUTPUT_VCF}.gz with {len(vcf_df):,} variants")


if __name__ == '__main__':
    main()
