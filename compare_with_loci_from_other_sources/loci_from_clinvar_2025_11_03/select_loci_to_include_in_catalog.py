#%%
import collections
import os
import pandas as pd

os.chdir(os.path.dirname(__file__))

table_path = "clinvar_2025_11_03.tandem_repeats.overlap_with_TRExplorer_v2.tsv.gz"


print(f"Processing {table_path}")
df = pd.read_table(table_path)
total = len(df)

df["motif_size"] = df["motif"].str.len()

df = df[df["overlap_score"] < 2].copy()

# remove loci that overlap with self (ie. with other loci in this table after filtering)
df.sort_values(by=["chrom", "start_0based", "end_1based"], inplace=True)
previous_row = None
rows_to_keep = []
zero_width_counter = 0
non_unique_counter = 0
for index, row in df.iterrows():
    if row["start_0based"] >= row["end_1based"]:
        zero_width_counter += 1
        continue
    is_unique=True
    if previous_row is not None and previous_row["chrom"] == row["chrom"]:
        intersection_size = min(row["end_1based"], previous_row["end_1based"]) - max(row["start_0based"], previous_row["start_0based"])
        union_size = max(row["end_1based"], previous_row["end_1based"]) - min(row["start_0based"], previous_row["start_0based"])
        jaccard_similarity = intersection_size / union_size
        if jaccard_similarity >= 0.2 or (jaccard_similarity > 0 and row["motif"] == previous_row["motif"]):
            is_unique = False
            non_unique_counter += 1

    if is_unique:
        previous_row = row
        rows_to_keep.append(row)

print(f"Keeping {len(rows_to_keep):,d} out of {total:,d} loci after removing {zero_width_counter:,d} zero-width loci and {non_unique_counter:,d} non-unique loci")
df = pd.DataFrame(rows_to_keep)

print(f"Selected {len(df):,d} out of {total:,d} loci with the following motif distribution:") 
print(df["motif_size"].value_counts())

output_filename_prefix = table_path.replace(".overlap_with_TRExplorer_v2.tsv.gz", "") 
output_tsv_path = f"{output_filename_prefix}.loci_to_include_in_catalog.tsv.gz"
output_bed_path = f"{output_filename_prefix}.loci_to_include_in_catalog.bed"

df.to_csv(output_tsv_path, sep="\t", index=False)
print(f"Wrote {len(df):,d} out of {total:,d} loci to {output_tsv_path}")

df_bed = df[["chrom", "start_0based", "end_1based", "motif"]]
df_bed.to_csv(output_bed_path, sep="\t", index=False, header=False)
os.system(f"bgzip -f {output_bed_path}")
os.system(f"tabix -f {output_bed_path}.gz")
print(f"Wrote {len(df):,d} out of {total:,d} loci to {output_bed_path}.gz")




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