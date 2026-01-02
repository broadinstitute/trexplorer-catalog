#%%
import os
import pandas as pd
import sys
from tqdm import tqdm
tqdm.pandas()

sys.path.append('../')
from compare_loci_utils import select_loci
#from compare_loci_with_catalog import OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BELOW_0_2

os.chdir(os.path.dirname(__file__))

table_path = "hg38.hipstr_reference.catalog.bed.gz"
output_path = "hg38.hipstr_reference.loci_to_include_in_catalog.bed"

print(f"Processing {table_path}")
df = pd.read_table(table_path, names=["chrom", "start_0based", "end_1based", "motif", "motif_size"])
total = len(df)

df = select_loci(df,
    #max_overlap_score=OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BELOW_0_2,
    max_locus_width=3_000,
    min_repeats_in_reference=2,
    min_adjusted_motif_purity=0.2,
    adjust_motifs_to_maximize_purity=True,
    drop_duplicates=True)

df[["chrom", "start_0based", "end_1based", "adjusted_motif"]].to_csv(
    output_path, index=False, header=False, sep="\t")

os.system(f"bgzip -f {output_path}")
os.system(f"tabix -f {output_path}.gz")

print("=" * 50)
print(f"Wrote {len(df):,d} loci to {output_path}.gz")
print("=" * 50)


"""
"overlap" column:

exact match                          930434
0.2 < Jaccard similarity ≤ 0.66      367844
diff ≤ 2 repeats                     153648
absent from TRExplorer_v2 catalog     86605  **
Jaccard similarity > 0.66             84080
Jaccard similarity ≤ 0.2              16332  **


Table columns:

['chrom',
 'start_0based',
 'end_1based',
 'motif',
 'simplified_motif',
 'canonical_motif',
 'reference_region',
 'reference_repeat_count',
 'reference_region_size',
 'HipSTR_start_0based',
 'HipSTR_end_1based',
 'HipSTR_motif',
 'HipSTR_simplified_motif',
 'HipSTR_canonical_motif',
 'HipSTR_reference_region',
 'HipSTR_reference_repeat_count',
 'HipSTR_reference_region_size',
 'HipSTR_purity_of_motif',
 'HipSTR_purity_of_motif_length',
 'optimal_motif',
 'optimal_canonical_motif',
 'optimal_motif_length',
 'optimal_motif_purity',
 'optimal_motif_quality_score',
 'optimal_motif_match',
 'overlap_score',
 'motif_match_score',
 'overlap',
 'motif_match',
 'match_found?',
 'TRExplorer_v2_start_0based',
 'TRExplorer_v2_end_1based',
 'TRExplorer_v2_motif',
 'TRExplorer_v2_canonical_motif',
 'TRExplorer_v2_reference_region',
 'TRExplorer_v2_reference_repeat_count',
 'TRExplorer_v2_reference_region_size',
 'TRExplorer_v2_purity_of_motif',
 'TRExplorer_v2_purity_of_motif_length']
"""
