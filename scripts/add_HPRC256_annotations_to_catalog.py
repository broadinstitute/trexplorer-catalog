"""Add HPRC256 population statistics annotations to a TR catalog JSON file.

HPRC256 is a large long-read sequencing cohort (Human Pangenome Reference Consortium)
used to compute allele frequency statistics at tandem repeat loci. This script annotates
catalog records with these population-level statistics.

This script owns every HPRC256_* field: each run clears them first and then writes only the
ones that apply, so running it again over an already-annotated catalog replaces the annotations
rather than layering onto them.
"""

import argparse
import gzip
import ijson
import math
import os
import pandas as pd
import simplejson as json
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.file_utils import download_local_copy
from catalog_annotation_utils import (clear_previous_annotations, find_unlisted_annotation_fields,
                                      print_annotation_replacement_summary)

DEFAULT_TSV_PATH = "gs://tandem-repeat-catalog/v2.0/hprc_lps.2025_12.per_locus_and_motif.256_samples.tsv.gz"

# Every field this script writes. A locus that drops out of the TSV, or whose row no longer has a
# value for one of the optional fields, must lose whatever a previous run gave it, so these are
# cleared per record before the new ones are written.
HPRC256_FIELDS = (
    "HPRC256_AlleleHistogram",
    "HPRC256_BiallelicHistogram",
    "HPRC256_MinAllele",
    "HPRC256_ModeAllele",
    "HPRC256_MaxAllele",
    "HPRC256_UniqueAlleleLengths",
    "HPRC256_NumCalledAlleles",
    "HPRC256_Stdev",
    "HPRC256_Median",
    "HPRC256_99thPercentile",
    "HPRC256_StdevRankByMotif",
    "HPRC256_StdevRankTotalNumberByMotif",
)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Add HPRC256 population statistics annotations to a TR catalog JSON file."
    )
    parser.add_argument("--tsv-path", default=DEFAULT_TSV_PATH,
                        help="Path to the HPRC256 TSV file (can be a GCS path)")
    parser.add_argument("--show-progress-bar", action="store_true",
                        help="Show a progress bar")
    parser.add_argument("-o", "--output-catalog-json-path",
                        help="Path of the output catalog JSON file with HPRC256 annotations")
    parser.add_argument("catalog_json_path", help="Path of the JSON catalog to annotate")
    args = parser.parse_args()

    # Download TSV if it's a remote path
    tsv_path = download_local_copy(args.tsv_path)
    if not os.path.isfile(tsv_path):
        parser.error(f"{args.tsv_path} file not found")

    if not os.path.isfile(args.catalog_json_path):
        parser.error(f"{args.catalog_json_path} file not found")

    if not args.output_catalog_json_path:
        args.output_catalog_json_path = args.catalog_json_path.replace(".json.gz", ".with_HPRC256_annotations.json.gz").replace(".json", ".with_HPRC256_annotations.json.gz")

    # Load and process the TSV
    print(f"Parsing {args.tsv_path}")
    df = pd.read_table(tsv_path)

    expected_columns = {"locus_id", "motif", "allele_size_histogram", "min_allele", "mode_allele",
                        "stdev", "median", "99th_percentile", "max_allele", "unique_allele_lengths",
                        "num_called_alleles"}
    missing_columns = expected_columns - set(df.columns)
    if missing_columns:
        parser.error(f"{args.tsv_path} is missing expected columns: {missing_columns}")

    print(f"Loaded {len(df):,d} rows from TSV")

    # The annotation lookup below is keyed by locus_id alone, so a duplicate locus_id would
    # silently overwrite an earlier row's annotations rather than raising. This is the failure mode
    # of https://github.com/PacificBiosciences/trgt-lps/issues/5, just one step downstream of it.
    duplicate_locus_ids = df["locus_id"][df["locus_id"].duplicated()].unique()
    if len(duplicate_locus_ids) > 0:
        parser.error(f"{args.tsv_path} has {len(duplicate_locus_ids):,d} duplicate locus_id values, "
                     f"e.g. {list(duplicate_locus_ids[:5])}")

    # Compute stdev rank by motif
    print("Computing stdev ranks by canonical motif")
    df["canonical_motif"] = df["motif"].apply(compute_canonical_motif)
    df_grouped_by_motif = df.groupby("canonical_motif")
    df["stdev_rank_by_motif"] = df_grouped_by_motif["stdev"].rank(ascending=False)
    df["stdev_rank_total_number_by_motif"] = df_grouped_by_motif["locus_id"].transform("count")

    # Build annotation lookup
    annotation_lookup = {}
    row_iterator = df.iterrows()
    if args.show_progress_bar:
        row_iterator = tqdm.tqdm(row_iterator, total=len(df), unit=" records", unit_scale=True,
                                 desc="Building annotation lookup")

    for _, row in row_iterator:
        locus_id = row["locus_id"]
        annotations = {}

        # String fields
        annotations["HPRC256_AlleleHistogram"] = row["allele_size_histogram"]
        biallelic_histogram = row.get("biallelic_histogram")
        if pd.notna(biallelic_histogram):
            annotations["HPRC256_BiallelicHistogram"] = biallelic_histogram

        # Integer fields
        if pd.notna(row["min_allele"]):
            annotations["HPRC256_MinAllele"] = int(row["min_allele"])
        if pd.notna(row["mode_allele"]):
            annotations["HPRC256_ModeAllele"] = int(row["mode_allele"])
        if pd.notna(row["max_allele"]):
            annotations["HPRC256_MaxAllele"] = int(row["max_allele"])
        if pd.notna(row["unique_allele_lengths"]):
            annotations["HPRC256_UniqueAlleleLengths"] = int(row["unique_allele_lengths"])
        if pd.notna(row["num_called_alleles"]):
            annotations["HPRC256_NumCalledAlleles"] = int(row["num_called_alleles"])

        # Float fields
        if pd.notna(row["stdev"]):
            annotations["HPRC256_Stdev"] = row["stdev"]
        if pd.notna(row["median"]):
            annotations["HPRC256_Median"] = row["median"]
        if pd.notna(row["99th_percentile"]):
            annotations["HPRC256_99thPercentile"] = row["99th_percentile"]

        # Computed rank fields
        if pd.notna(row["stdev_rank_by_motif"]):
            annotations["HPRC256_StdevRankByMotif"] = int(row["stdev_rank_by_motif"])
        if pd.notna(row["stdev_rank_total_number_by_motif"]):
            annotations["HPRC256_StdevRankTotalNumberByMotif"] = int(row["stdev_rank_total_number_by_motif"])

        annotation_lookup[locus_id] = annotations

    print(f"Built annotation lookup with {len(annotation_lookup):,d} entries")

    unlisted_fields = find_unlisted_annotation_fields(annotation_lookup, HPRC256_FIELDS)
    if unlisted_fields:
        parser.error(f"HPRC256_FIELDS doesn't list {sorted(unlisted_fields)}, so re-running this script "
                     f"over an already-annotated catalog would leave those fields stale")

    # Annotate the catalog
    input_locus_counter = annotated_locus_counter = 0
    locus_with_previous_annotation_counter = locus_that_lost_annotation_counter = 0
    print(f"Adding HPRC256 annotations to {args.catalog_json_path}")
    fopen = gzip.open if args.catalog_json_path.endswith("gz") else open
    with fopen(args.catalog_json_path, "rt") as f:
        f2open = gzip.open if args.output_catalog_json_path.endswith("gz") else open
        with f2open(args.output_catalog_json_path, "wt") as f2:
            iterator = ijson.items(f, "item")
            if args.show_progress_bar:
                iterator = tqdm.tqdm(iterator, unit=" records", unit_scale=True,
                                     desc="Annotating catalog")
            f2.write("[")
            for i, record in enumerate(iterator):
                locus_id = record["LocusId"]
                input_locus_counter += 1

                had_previous_annotation = clear_previous_annotations(record, HPRC256_FIELDS)
                if had_previous_annotation:
                    locus_with_previous_annotation_counter += 1

                if locus_id in annotation_lookup:
                    record.update(annotation_lookup[locus_id])
                    annotated_locus_counter += 1
                elif had_previous_annotation:
                    locus_that_lost_annotation_counter += 1
                if i > 0:
                    f2.write(", ")
                f2.write(json.dumps(record, use_decimal=True, indent=4))
            f2.write("]")

    print(f"Annotated {annotated_locus_counter:,d} out of {input_locus_counter:,d} loci "
          f"({annotated_locus_counter/max(1, input_locus_counter):.1%})")
    print_annotation_replacement_summary("HPRC256", f"are no longer in {args.tsv_path}",
                                         locus_with_previous_annotation_counter,
                                         locus_that_lost_annotation_counter)
    print(f"Wrote annotated catalog to {args.output_catalog_json_path}")


if __name__ == "__main__":
    main()
