#%%
import os
import pandas as pd
import pyfaidx

from str_analysis.utils.find_motif_utils import find_highest_purity_motif_length_for_interval


os.chdir(os.path.dirname(__file__))

table_path = "hg38.hipstr_reference.catalog.overlap_with_TRExplorer_v2.tsv.gz"

output_path = "hg38.hipstr_reference.loci_to_include_in_catalog.bed"


pyfaidx_reference_fasta_obj = pyfaidx.Fasta(os.path.expanduser("~/hg38.fa"), one_based_attributes=False, as_raw=True)

print(f"Processing {table_path}")
df = pd.read_table(table_path)

from pprint import pprint
pprint(df["overlap"].value_counts())
df = df[df["overlap"].isin({
    "Jaccard similarity ≤ 0.2",
    "absent from TRExplorer_v2 catalog",
})].copy()

def compute_highest_purity_motif(row):
    optimal_motif, _, _ = find_highest_purity_motif_length_for_interval(
        pyfaidx_reference_fasta_obj, str(row.chrom), row.HipSTR_start_0based, row.HipSTR_end_1based,
        min_motif_length=len(row.HipSTR_motif), max_motif_length=len(row.HipSTR_motif))

    return optimal_motif


df["HighestPurityMotifOfSameLengthAsOriginal"] = df.apply(compute_highest_purity_motif, axis=1)

print(f"Will change the motif for", sum(df["HighestPurityMotifOfSameLengthAsOriginal"] != df["HipSTR_motif"]), "out of",
      len(df), "loci")

df[["chrom", "HipSTR_start_0based",  "HipSTR_end_1based", "HighestPurityMotifOfSameLengthAsOriginal"]].to_csv(
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
