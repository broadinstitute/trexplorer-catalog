import os
import pandas as pd

os.chdir(os.path.dirname(__file__))

# Read supplementary table S1 from Tanudisastro et al. 2024
# Paper: "Polymorphic tandem repeats influence cell type-specific gene expression across the human immune landscape"
# URL: https://www.biorxiv.org/content/10.1101/2024.11.02.621562v2
#
# Table S1 contains 55,310 unique sc-eTR loci (single-cell expression quantitative trait TR loci)
# Coordinates are already in hg38, so no liftover needed

df = pd.read_table("./TableS1v0.1.tsv")
print(f"Read {len(df):,d} rows from ./TableS1v0.1.tsv")

# Extract unique loci (chr, start, end, motif)
# The table has multiple rows per locus (one per cell type and eGene association)
df_unique = df[['chromosome', 'start coordinate (hg38)', 'end coordinate (hg38)', 'motif']].drop_duplicates()
print(f"Found {len(df_unique):,d} unique loci")

output_rows = []
for _, row in df_unique.iterrows():
    chrom = row['chromosome']
    start_0based = row['start coordinate (hg38)']
    end_1based = row['end coordinate (hg38)']
    motif = row['motif']
    output_rows.append((chrom, start_0based, end_1based, motif))

output_bed_path = "Tanudisastro_2024_loci.bed"
with open(output_bed_path, "wt") as output_bed:
    for output_row in sorted(output_rows):
        chrom, start_0based, end_1based, motif = output_row
        output_bed.write(f"{chrom}\t{start_0based}\t{end_1based}\t{motif}\n")

os.system(f"bgzip -f {output_bed_path}")
os.system(f"tabix -f {output_bed_path}.gz")
print(f"Wrote {len(output_rows):,d} loci to {output_bed_path}.gz")
