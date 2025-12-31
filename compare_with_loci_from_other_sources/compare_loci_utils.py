import os
import pandas as pd
import pyfaidx

from str_analysis.utils.find_motif_utils import adjust_motif_to_maximize_purity_in_interval


def compute_adjusted_motif_columns(df, reference_genome_path="~/hg38.fa"):
    print(f"Computing adjusted motifs for {len(df):,d} loci")
    pyfaidx_reference_fasta_obj = pyfaidx.Fasta(os.path.expanduser(reference_genome_path), one_based_attributes=False, as_raw=True)

    def adjust_locus(row):
        try:
            adjusted_motif, purity = adjust_motif_to_maximize_purity_in_interval(
                pyfaidx_reference_fasta_obj,
                f"chr{str(row.chrom).replace('chr', '')}",
                row["start_0based"],
                row["end_1based"],
                row["motif"],
            )

            return pd.Series({
                "adjusted_motif": adjusted_motif,
                "was_adjusted": adjusted_motif != row["motif"],
                "adjusted_motif_purity": purity
            })
        except ValueError as e:
            print(f"WARNING: {e}")
            return pd.Series({
                "adjusted_motif": row["motif"],
                "was_adjusted": False,
                "adjusted_motif_purity": pd.NA,
            })

    df[["adjusted_motif", "was_adjusted", "adjusted_motif_purity"]] = df.apply(adjust_locus, axis=1)
    print(f"Adjusted {df['was_adjusted'].sum():,d} out of {len(df):,d} loci")

    return df


def select_loci(
    df,
    adjust_motifs_to_maximize_purity=True,
    max_overlap_score=None,
    max_locus_width=None,
    min_repeats_in_reference=None,
    min_adjusted_motif_purity=None,
    drop_duplicates=True,
    reference_genome_path="~/hg38.fa",
):

    if max_overlap_score is not None:
        before_filter = len(df)
        df = df[df["overlap_score"] <= max_overlap_score].copy()
        print(f"After filtering to overlap_score <= {max_overlap_score}: kept {len(df):,d} out of {before_filter:,d} ({len(df)/before_filter:.1%}) loci")

    if max_locus_width is not None:
        before_filter = len(df)
        df = df[(df["end_1based"] - df["start_0based"]) < max_locus_width].copy()
        print(f"After filtering to locus width < {max_locus_width}bp: kept {len(df):,d} out of {before_filter:,d} ({len(df)/before_filter:.1%}) loci")

    if min_repeats_in_reference is not None:
        before_filter = len(df)
        df["repeat_count"] = (df["end_1based"] - df["start_0based"]) / df["motif"].str.len()
        df = df[df["repeat_count"] >= min_repeats_in_reference].copy()
        print(f"After filtering to >= {min_repeats_in_reference} repeats: kept {len(df):,d} out of {before_filter:,d} ({len(df)/before_filter:.1%}) loci")

    if adjust_motifs_to_maximize_purity:
        df = compute_adjusted_motif_columns(df, reference_genome_path=reference_genome_path)
        df["motif"] = df["adjusted_motif"]

    if min_adjusted_motif_purity is not None:
        nan_purity_count = df["adjusted_motif_purity"].isna().sum()
        if nan_purity_count > 0:
            print(f"Warning: {nan_purity_count:,d} loci have NaN purity and will be excluded")

        before_filter = len(df)
        df = df[df["adjusted_motif_purity"] >= min_adjusted_motif_purity].copy()
        print(f"After filtering for purity > {min_adjusted_motif_purity}: kept {len(df):,d} out of {before_filter:,d} ({len(df)/before_filter:.1%}) loci")

    if drop_duplicates:
        before_filter = len(df)
        df = df.drop_duplicates(subset=["chrom", "start_0based", "end_1based", "motif"])
        if len(df) != before_filter:
            print(f"Dropped {before_filter - len(df):,d} duplicate loci")

    if any(df["motif"].str.upper().str.contains("N")):
        raise ValueError("Some motifs contain 'N'")

    return df
