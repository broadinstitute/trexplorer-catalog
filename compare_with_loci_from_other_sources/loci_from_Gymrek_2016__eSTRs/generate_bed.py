#!/usr/bin/env python3
"""
Generate BED file from Gymrek 2016 eSTR data.

This script processes the supplementary CSV file from:
Gymrek, M., Willems, T., Guilmatre, A. et al.
Abundant contribution of short tandem repeats to gene expression variation in humans.
Nat Genet 48, 22–29 (2016).
https://doi.org/10.1038/ng.3461

The paper identified 2,060 significant expression STRs (eSTRs) at gene-level FDR of 5%.
Coordinates are in hg19 and will be converted to hg38.
"""

from liftover import get_lifter
import os
import pandas as pd
from pprint import pprint

os.chdir(os.path.dirname(__file__))

# TODO: Update this filename after downloading the supplementary CSV file
# The file should be named something like "ng.3461-S3.csv" or similar
# Download from: https://www.nature.com/articles/ng.3461 (Supplementary Information section)
INPUT_FILE = "REPLACE_WITH_ACTUAL_FILENAME.csv"

# Check if file exists
if not os.path.exists(INPUT_FILE):
    print(f"ERROR: Input file '{INPUT_FILE}' not found!")
    print("Please download the supplementary CSV file from:")
    print("https://www.nature.com/articles/ng.3461")
    print("Update the INPUT_FILE variable in this script with the actual filename.")
    exit(1)

# Initialize liftOver from hg19 to hg38
converter = get_lifter('hg19', 'hg38', one_based=False)

# Read the eSTR data
# TODO: Update column names after inspecting the actual file
# Expected columns: chromosome, position/start, end(?), motif/repeat_unit, gene, p-value, etc.
df = pd.read_csv(INPUT_FILE)
print(f"Read {len(df):,d} loci from {INPUT_FILE}")

# Print first row to understand the structure
if len(df) > 0:
    print("\nFirst row:")
    pprint(df.iloc[0].to_dict())
    print("\nColumns:", df.columns.tolist())

# TODO: Uncomment and adapt this code after inspecting the file structure
"""
output_rows = []
skipped_count = 0

for _, row in df.iterrows():
    # TODO: Update these column names based on actual CSV structure
    chrom_hg19 = row["chromosome"]  # or "chr" or "chrom"

    # Ensure chromosome has 'chr' prefix
    if not chrom_hg19.startswith('chr'):
        chrom_hg19 = f'chr{chrom_hg19}'

    # TODO: Update based on how position is specified in the file
    # It might be a single position, or start/end coordinates
    start_hg19_0based = row["start"] - 1  # Convert to 0-based if 1-based in file
    # end_hg19_1based = row["end"]  # if end column exists

    # Get the repeat motif
    motif = row["motif"]  # or "repeat_unit" or similar column name

    # Perform liftOver
    liftover_result = converter[chrom_hg19][start_hg19_0based]

    if len(liftover_result) != 1 or len(liftover_result[0]) != 3:
        print(f"WARNING: Liftover result for {chrom_hg19}:{start_hg19_0based} has {len(liftover_result)} elements: {liftover_result}. Skipping...")
        skipped_count += 1
        continue

    chrom_hg38, start_hg38_0based, strand_hg38 = liftover_result[0]

    # TODO: Calculate end coordinate based on the data
    # This might be: end = start + (motif_length * num_repeats)
    # or the end coordinate might be provided in the file
    end_hg38_1based = start_hg38_0based + len(motif) * 10  # PLACEHOLDER - update based on actual data

    output_rows.append((chrom_hg38, start_hg38_0based, end_hg38_1based, motif))

print(f"\nProcessed {len(output_rows):,d} loci successfully")
print(f"Skipped {skipped_count:,d} loci due to liftOver issues")

# Write output BED file
output_bed_path = "Gymrek_2016_eSTRs.bed"
with open(output_bed_path, "wt") as output_bed:
    for output_row in sorted(output_rows):
        chrom_hg38, start_hg38_0based, end_hg38_1based, motif = output_row
        output_bed.write(f"{chrom_hg38}\t{start_hg38_0based}\t{end_hg38_1based}\t{motif}\n")

# Compress and index
os.system(f"bgzip -f {output_bed_path}")
os.system(f"tabix -f {output_bed_path}.gz")
print(f"\nWrote {len(output_rows):,d} loci to {output_bed_path}.gz")
"""

print("\n" + "="*80)
print("SCRIPT TEMPLATE - NEEDS TO BE COMPLETED")
print("="*80)
print("1. Download the supplementary CSV file from the paper")
print("2. Update INPUT_FILE variable with the actual filename")
print("3. Inspect the file structure by running this script")
print("4. Uncomment and update the processing code based on actual column names")
print("5. Test the script to ensure correct BED file generation")
print("="*80)
