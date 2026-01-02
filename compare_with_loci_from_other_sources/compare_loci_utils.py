"""Utilities for filtering and selecting tandem repeat loci from external sources.

This module provides functions to filter loci based on various criteria (width,
repeat count, purity) and to adjust motifs to maximize purity against a reference
genome. These utilities are used by the select_loci_to_include_in_catalog.py
scripts in each loci_from_* subdirectory.

Typical usage:
    from compare_loci_utils import select_loci

    df = pd.read_table("loci.bed", names=["chrom", "start_0based", "end_1based", "motif"])
    df = select_loci(df, min_repeats_in_reference=2, min_adjusted_motif_purity=0.2)
"""

import os
import pandas as pd
import pyfaidx

from str_analysis.utils.find_motif_utils import adjust_motif_to_maximize_purity_in_interval

DEFAULT_REFERENCE_GENOME_PATH = "~/hg38.fa"

FILTER_LOCUS_ON_SUPERCONTIG = "Locus on supercontig"
FILTER_LOCUS_WIDTH_TOO_LARGE = "Locus width too large"
FILTER_MOTIF_EMPTY = "Motif is empty"
FILTER_REPEAT_COUNT_TOO_LOW = "Repeat count too low"
FILTER_PURITY_CALCULATION_ERROR = "Error calculating purity"
FILTER_PURITY_TOO_LOW = "Purity too low"


def does_locus_pass_filters(chrom, start_0based, end_1based, motif, reference_genome_path=DEFAULT_REFERENCE_GENOME_PATH,
    max_locus_width=None,
    min_repeats_in_reference=None,
    min_adjusted_motif_purity=None,
    keep_only_primary_chromosomes=True,
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
        keep_only_primary_chromosomes: Filter out loci on super-contigs

    Returns:
        tuple: A 3-tuple containing:
            - passes_filters (bool): True if the locus passes all filters
            - adjusted_motif (str or None): The adjusted motif if passes, None otherwise
            - filter_reason (str or None): The reason for filtering if fails, None otherwise
    """
    locus_width = end_1based - start_0based

    if keep_only_primary_chromosomes and "_" in chrom:
        return False, None, FILTER_LOCUS_ON_SUPERCONTIG

    if max_locus_width is not None:
        if locus_width >= max_locus_width:
            return False, None, FILTER_LOCUS_WIDTH_TOO_LARGE

    if min_repeats_in_reference is not None:
        if not motif:
            return False, None, FILTER_MOTIF_EMPTY
        repeat_count = locus_width / len(motif)
        if repeat_count < min_repeats_in_reference:
            return False, None, FILTER_REPEAT_COUNT_TOO_LOW

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
        return False, None, FILTER_PURITY_CALCULATION_ERROR
    finally:
        if pyfaidx_reference_fasta_obj is not None:
            pyfaidx_reference_fasta_obj.close()

    if min_adjusted_motif_purity is not None:
        if purity < min_adjusted_motif_purity:
            return False, None, FILTER_PURITY_TOO_LOW

    if adjusted_motif.upper().strip("ACGT") != "":
        raise ValueError(f"adjusted motif contains unexpected chars: {chrom}:{start_0based}-{end_1based} {motif} => {adjusted_motif}")

    return True, adjusted_motif, None


def compute_adjusted_motif_columns(df, reference_genome_path=DEFAULT_REFERENCE_GENOME_PATH):
    """Add columns for adjusted motif and purity to a DataFrame of loci.

    For each locus, attempts to find a rotation of the motif that maximizes
    the purity (fraction of bases matching the repeat pattern) in the reference
    genome sequence at that locus.

    Args:
        df: DataFrame with columns 'chrom', 'start_0based', 'end_1based', 'motif'.
        reference_genome_path: Path to reference genome FASTA file.

    Returns:
        The input DataFrame with three new columns added:
            - adjusted_motif: The motif rotation with highest purity
            - was_adjusted: True if adjusted_motif differs from original motif
            - adjusted_motif_purity: Fraction of bases matching the repeat pattern
    """
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
    keep_only_primary_chromosomes=True,
    reference_genome_path=DEFAULT_REFERENCE_GENOME_PATH,
):
    """Filter and process a DataFrame of tandem repeat loci.

    Applies a series of filters to select loci that meet quality criteria,
    optionally adjusts motifs to maximize purity, and removes duplicates.
    Progress is printed to stdout as each filter is applied.

    Args:
        df: DataFrame with columns 'chrom', 'start_0based', 'end_1based', 'motif'.
            May also contain 'overlap_score' if max_overlap_score is used.
        adjust_motifs_to_maximize_purity: If True, rotate motifs to maximize
            purity against the reference genome sequence.
        max_overlap_score: If set, exclude loci with overlap_score above this
            value. Requires 'overlap_score' column in df.
        max_locus_width: If set, exclude loci wider than this value (in bp).
        min_repeats_in_reference: If set, exclude loci with fewer than this
            many repeats (locus_width / motif_length).
        min_adjusted_motif_purity: If set, exclude loci with adjusted purity
            below this value (0.0 to 1.0).
        drop_duplicates: If True, remove duplicate loci based on coordinates
            and motif.
        keep_only_primary_chromosomes: If True, exclude loci on supercontigs
            (chromosomes containing '_' in their name).
        reference_genome_path: Path to reference genome FASTA file.

    Returns:
        Filtered DataFrame with the same columns as input, plus additional
        columns if adjust_motifs_to_maximize_purity is True:
            - adjusted_motif, was_adjusted, adjusted_motif_purity

    Raises:
        ValueError: If any remaining motifs contain 'N'.
    """
    if keep_only_primary_chromosomes:
        before_filter = len(df)
        df = df[~df["chrom"].astype(str).str.contains("_")].copy()
        if len(df) != before_filter:
            print(f"After filtering to primary chromosomes: kept {len(df):,d} out of {before_filter:,d} ({len(df)/before_filter:.1%}) loci")

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

