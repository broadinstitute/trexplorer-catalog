#%%
import os
import pandas as pd
import pyfaidx
import sys
sys.path.append('../')
from tqdm import tqdm

tqdm.pandas()


from compare_loci_with_catalog import OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BELOW_0_2

from str_analysis.utils.find_motif_utils import adjust_motif_to_maximize_purity_in_interval


os.chdir(os.path.dirname(__file__))

table_path = "hg38.hipstr_reference.catalog.overlap_with_TRExplorer_v2.tsv.gz"

output_path = "hg38.hipstr_reference.loci_to_include_in_catalog.bed"


pyfaidx_reference_fasta_obj = pyfaidx.Fasta(os.path.expanduser("~/hg38.fa"), one_based_attributes=False, as_raw=True)

print(f"Processing {table_path}")
df = pd.read_table(table_path)

df = df[df["overlap_score"] <= OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BELOW_0_2].copy()
df = df[(df["end_1based"] - df["start_0based"]) < 3_000]

print(f"Adjusting boundaries and motifs for {len(df):,d} loci...")
def adjust_locus(row):
    try:
        adjusted_motif, purity = adjust_motif_to_maximize_purity_in_interval(
            pyfaidx_reference_fasta_obj,
            f"chr{str(row.chrom).replace('chr', '')}",
            row.HipSTR_start_0based,
            row.HipSTR_end_1based,
            row.HipSTR_motif,
        )
        was_adjusted = adjusted_motif != row.HipSTR_motif
        return pd.Series({
            'adjusted_motif': adjusted_motif,
            'was_adjusted': was_adjusted,
            'purity': purity
        })

    except ValueError as e:
        print(f"WARNING: {e}")
        return pd.Series({
            'adjusted_motif': row.HipSTR_motif,
            'was_adjusted': False,
            'purity': float('nan'),
        })

df[['adjusted_motif', 'was_adjusted', 'purity']] = df.progress_apply(
    adjust_locus, axis=1)

print(f"Adjusted {df['was_adjusted'].sum():,d} out of {len(df):,d} loci")

df[["chrom", "HipSTR_start_0based", "HipSTR_end_1based", "adjusted_motif"]].to_csv(
    output_path, index=False, header=False, sep="\t")

os.system(f"bgzip -f {output_path}")
os.system(f"tabix -f {output_path}.gz")

print(f"Wrote {len(df):,d} loci to {output_path}.gz")


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
