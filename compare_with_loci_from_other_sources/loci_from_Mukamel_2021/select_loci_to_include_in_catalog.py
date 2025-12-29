#%%
import os
import pandas as pd
import pyfaidx

import sys
sys.path.append('../')
from compare_loci_with_catalog import OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_ABOVE_0_66
from str_analysis.utils.find_motif_utils import adjust_motif_and_boundaries_to_maximize_purity

os.chdir(os.path.dirname(__file__))

table_paths = [
    "vntrs_in_ST1.overlap_with_TRExplorer_v2.tsv.gz",
]

for table_path in table_paths:
    print("-" * 100)
    print(f"Processing {table_path}")
    df = pd.read_table(table_path)
    total = len(df)
    df["motif_size"] = df["motif"].str.len()

    # filter to loci that don't match anything in the catalog and have motif size > 2
    df = df[df["overlap_score"] < OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_ABOVE_0_66].copy()

    print(f"Selected {len(df):,d} out of {total:,d} loci with the following motif distribution:")
    print(df["motif_size"].value_counts())

    # Load reference genome and adjust boundaries/motifs
    print("Loading reference genome...")
    pyfaidx_reference_fasta_obj = pyfaidx.Fasta(os.path.expanduser("~/hg38.fa"), one_based_attributes=False, as_raw=True)

    print(f"Adjusting boundaries and motifs for {len(df):,d} loci...")
    def adjust_locus(row):
        adjusted_start, adjusted_end, adjusted_motif, was_adjusted, purity = adjust_motif_and_boundaries_to_maximize_purity(
            pyfaidx_reference_fasta_obj, str(row.chrom), row.start_0based, row.end_1based, row.motif
        )
        return pd.Series({
            'adjusted_start_0based': adjusted_start,
            'adjusted_end_1based': adjusted_end,
            'adjusted_motif': adjusted_motif,
            'was_adjusted': was_adjusted,
            'purity': purity
        })

    df[['adjusted_start_0based', 'adjusted_end_1based', 'adjusted_motif', 'was_adjusted', 'purity']] = df.apply(adjust_locus, axis=1)

    print(f"Adjusted {df['was_adjusted'].sum():,d} out of {len(df):,d} loci")

    # Update the main columns with adjusted values
    df['start_0based'] = df['adjusted_start_0based']
    df['end_1based'] = df['adjusted_end_1based']
    df['motif'] = df['adjusted_motif']

    output_filename_prefix = table_path.replace(".overlap_with_TRExplorer_v2.tsv.gz", "") 
    output_tsv_path = f"{output_filename_prefix}.loci_to_include_in_catalog.tsv.gz"
    output_bed_path = f"{output_filename_prefix}.loci_to_include_in_catalog.bed"
    
    df.to_csv(output_tsv_path, sep="\t", index=False)
    print(f"Wrote {len(df):,d} out of {total:,d} loci to {output_tsv_path}")

    df_bed = df[["chrom", "start_0based", "end_1based", "motif"]]
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