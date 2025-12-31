#%%
import os
import pandas as pd

import sys
sys.path.append('../')
from compare_loci_utils import select_loci
#from compare_loci_with_catalog import OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_ABOVE_0_66

os.chdir(os.path.dirname(__file__))

table_path = "functional_vntrs.bed"

print("-" * 100)
print(f"Processing {table_path}")
df = pd.read_table(table_path, names=["chrom", "start_0based", "end_1based", "motif"])
total = len(df)

# Filter to loci that don't match well in catalog (overlap_score < 5, equivalent to <= 4)
df = select_loci(df,
    #max_overlap_score=OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_ABOVE_0_66 - 1,
    min_repeats_in_reference=2,
    min_adjusted_motif_purity=0.2,
    adjust_motifs_to_maximize_purity=True,
    drop_duplicates=True)

output_filename_prefix = "known_functional_VNTRs"
output_tsv_path = f"{output_filename_prefix}.loci_to_include_in_catalog.tsv.gz"
output_bed_path = f"{output_filename_prefix}.loci_to_include_in_catalog.bed"

df.to_csv(output_tsv_path, sep="\t", index=False)
print(f"Wrote {len(df):,d} out of {total:,d} loci to {output_tsv_path}")

df_bed = df[["chrom", "start_0based", "end_1based", "adjusted_motif"]]
df_bed.sort_values(by=["chrom", "start_0based", "end_1based"], inplace=True)
df_bed.to_csv(output_bed_path, sep="\t", index=False, header=False)
os.system(f"bgzip -f {output_bed_path}")
os.system(f"tabix -f {output_bed_path}.gz")
print(f"Wrote {len(df):,d} out of {total:,d} loci to {output_bed_path}.gz")

print(f"Selected {len(df):,d} out of {total:,d} loci with the following motif distribution:")
df["motif_size"] = df["motif"].str.len()
print(df["motif_size"].value_counts())




"""

set(zip(df.overlap_score, df.overlap))

{(0, 'absent from TRExplorer_v2 catalog'),   # match found? is no
 (1, 'Jaccard similarity ≤ 0.2'),            # match found? is no
 (2, '0.2 < Jaccard similarity ≤ 0.33'),
 (3, '0.33 < Jaccard similarity ≤ 0.5'),
 (4, '0.5 < Jaccard similarity ≤ 0.66'),
 (5, 'Jaccard similarity > 0.66'),
 (6, 'diff ≤ 2 repeats'),
 (7, 'exact match')}

 set(zip(df.motif_match_score, df.motif_match))

{(0, 'absent from TRExplorer_v2 catalog'),
 (1, 'different motif length'),
 (2, 'same motif length'),
 (3, 'same motif')}

Columns:
    chrom, 
    overlap_score,
    motif_match_score, 
    overlap, 
    motif_match, 
    match_found?,
    TRExplorer_v2_start_0based, 
    TRExplorer_v2_end_1based,
    TRExplorer_v2_canonical_motif, 
    TRExplorer_v2_reference_region,
    TRExplorer_v2_reference_repeat_count,
    TRExplorer_v2_reference_region_size
"""

#%%
