#%%
import os
import pandas as pd

os.chdir(os.path.dirname(__file__))

table_paths = [
    "TableS1v0.1.overlap_with_TRExplorer_v2.tsv.gz",
    #"TableS3v0.1.overlap_with_TRExplorer_v2.tsv.gz",
    #"fm_etrs_coordinates.overlap_with_TRExplorer_v2.tsv.gz",
]

for table_path in table_paths:
    print("-" * 100)
    print(f"Processing {table_path}")
    df = pd.read_table(table_path)
    total = len(df)
    df["motif_size"] = df["motif"].str.len()

    df = df.drop_duplicates(subset=["chrom", "start_0based", "end_1based", "motif"])
    if len(df) != total:
        print(f"Dropped {total - len(df):,d} duplicate loci")

    # filter to loci that don't match anything in the catalog and have motif size > 2
    df = df[df["match_found?"].str.startswith("no")].copy()
    df = df[df["motif_size"] > 2].copy()


    assert set(df.overlap_score) == {0, 1}, f"Expected overlap_score to be 0 or 1. Found scores: {set(df.overlap_score)}"
    print(f"Selected {len(df):,d} out of {total:,d} loci with the following motif distribution:") 

    # if motif size > 7, then set motif size to 7
    df.loc[df["motif_size"] > 7, "motif_size"] = 7
    print(df["motif_size"].value_counts())

    output_filename_prefix = table_path.replace(".overlap_with_TRExplorer_v2.tsv.gz", "") 
    output_tsv_path = f"{output_filename_prefix}.loci_to_include_in_catalog.tsv.gz"
    output_bed_path = f"{output_filename_prefix}.loci_to_include_in_catalog.bed"
    
    df.to_csv(output_tsv_path, sep="\t", index=False)
    print(f"Wrote {len(df):,d} out of {total:,d} loci to {output_tsv_path}")

    df_bed = df[["chrom", "start_0based", "end_1based", "motif"]].copy()
    df_bed.sort_values(by=["chrom", "start_0based", "end_1based"], inplace=True)
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
    sc-eTRs_from_TableS1_start_0based,
    sc-eTRs_from_TableS1_end_1based, 
    sc-eTRs_from_TableS1_motif,
    sc-eTRs_from_TableS1_simplified_motif,
    sc-eTRs_from_TableS1_canonical_motif,
    sc-eTRs_from_TableS1_reference_region,
    sc-eTRs_from_TableS1_reference_repeat_count,
    sc-eTRs_from_TableS1_reference_region_size, 
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