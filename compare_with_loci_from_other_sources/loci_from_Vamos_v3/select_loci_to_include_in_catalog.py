"""Filters Vamos TR loci based on quality criteria for catalog inclusion.

This script reads the BED file produced by convert_to_bed.py and applies
filtering criteria to select high-quality TR loci suitable for inclusion
in the tandem repeat catalog.

Filtering criteria:
    - Minimum motif size: 7bp (VNTRs only)
    - Maximum locus width: 10,000bp
    - Minimum repeats in reference: 2
    - Minimum adjusted motif purity: 0.33
    - Only ACGT bases in motifs (no ambiguous bases)

Output:
    vamosGenomicTR_v3.0.loci_to_include_in_catalog.bed.gz
"""

import os
import sys

import pandas as pd
from tqdm import tqdm

tqdm.pandas()

sys.path.append('../')
from compare_loci_utils import select_loci

os.chdir(os.path.dirname(__file__))

table_path = "vamosGenomicTR_v3.0.bed.gz"
output_path = "vamosGenomicTR_v3.0.loci_to_include_in_catalog.bed"

print(f"Processing {table_path}")
df = pd.read_table(table_path, names=["chrom", "start_0based", "end_1based", "motif", "motif_size"])
total = len(df)
print(f"Loaded {total:,d} loci")

df = select_loci(
    df,
    min_motif_size=7,
    max_locus_width=10_000,
    min_repeats_in_reference=2,
    #min_adjusted_motif_purity=0.33,
    min_adjusted_motif_purity=0.85,    
    adjust_motifs_to_maximize_purity=True,
    drop_duplicates=True,
    keep_only_motifs_with_ACGT_bases=True,
)

df[["chrom", "start_0based", "end_1based", "adjusted_motif"]].to_csv(
    output_path, index=False, header=False, sep="\t")

os.system(f"bgzip -f {output_path}")
os.system(f"tabix -f {output_path}.gz")

print("=" * 50)
print(f"Wrote {len(df):,d} loci to {output_path}.gz ({len(df)/total*100:.1f}% of input)")
print("=" * 50)
