import os
import pandas as pd
import pyfaidx

from str_analysis.utils.find_motif_utils import adjust_motif_to_maximize_purity_in_interval

DEFAULT_REFERENCE_GENOME_PATH = "~/hg38.fa"



def does_locus_pass_filters(chrom, start_0based, end_1based, motif, reference_genome_path=DEFAULT_REFERENCE_GENOME_PATH,
    max_locus_width=None,
    min_repeats_in_reference=None,
    min_adjusted_motif_purity=None,
):
    """Check if a locus passes all specified filters.

    Args:
        chrom: Chromosome name
        start_0based: 0-based start coordinate
        end_1based: 1-based end coordinate
        motif: The repeat motif
        reference_genome_path: Path to reference genome FASTA
        max_locus_width: Maximum allowed locus width in bp (exclusive)
        min_repeats_in_reference: Minimum number of repeats required
        min_adjusted_motif_purity: Minimum purity after motif adjustment

    Returns:
        bool: True if the locus passes all filters, False otherwise
        str: Adjusted motif
    """
    locus_width = end_1based - start_0based

    if max_locus_width is not None:
        if locus_width >= max_locus_width:
            return False, None

    if min_repeats_in_reference is not None:
        if not motif:
            return False, None
        repeat_count = locus_width / len(motif)
        if repeat_count < min_repeats_in_reference:
            return False, None

    pyfaidx_reference_fasta_obj = None
    try:
        pyfaidx_reference_fasta_obj = pyfaidx.Fasta(
            os.path.expanduser(reference_genome_path), one_based_attributes=False, as_raw=True
        )

        adjusted_motif, purity = adjust_motif_to_maximize_purity_in_interval(
            pyfaidx_reference_fasta_obj,
            f"chr{str(chrom).replace('chr', '')}",
            start_0based,
            end_1based,
            motif,
        )

    except ValueError as e:
        print(f"WARNING: error while checking filters for locus: {chrom}:{start_0based}-{end_1based} {motif}: {e}. "
              f"Skipping locus..")
        return False, None
    finally:
        if pyfaidx_reference_fasta_obj is not None:
            pyfaidx_reference_fasta_obj.close()

    if min_adjusted_motif_purity is not None:
        if purity < min_adjusted_motif_purity:
            return False, None

    if adjusted_motif.upper().strip("ACGT") != "":
        raise ValueError(f"adjusted motif contains unexpected chars: {chrom}:{start_0based}-{end_1based} {motif} => {adjusted_motif}")

    return True, adjusted_motif


def compute_adjusted_motif_columns(df, reference_genome_path=DEFAULT_REFERENCE_GENOME_PATH):
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
    reference_genome_path=DEFAULT_REFERENCE_GENOME_PATH,
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
        # Filter out empty motifs to avoid division by zero
        df = df[df["motif"].str.len() > 0].copy()
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
        print(f"After filtering for purity >= {min_adjusted_motif_purity}: kept {len(df):,d} out of {before_filter:,d} ({len(df)/before_filter:.1%}) loci")

    if drop_duplicates:
        before_filter = len(df)
        df = df.drop_duplicates(subset=["chrom", "start_0based", "end_1based", "motif"])
        if len(df) != before_filter:
            print(f"Dropped {before_filter - len(df):,d} duplicate loci")

    if any(df["motif"].str.upper().str.contains("N")):
        raise ValueError("Some motifs contain 'N'")

    return df

