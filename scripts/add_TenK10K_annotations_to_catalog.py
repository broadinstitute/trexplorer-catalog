"""Add TenK10K population statistics annotations to a TR catalog JSON file.

TenK10K is a large short-read sequencing cohort used to compute allele frequency statistics
at tandem repeat loci. This script annotates catalog records with these population-level statistics.
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

DEFAULT_TSV_PATH = "gs://tandem-repeat-catalog/v2.0/tenk10k_str_mt_rows.reformatted.tsv.gz"


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Add TenK10K population statistics annotations to a TR catalog JSON file."
    )
    parser.add_argument("--tsv-path", default=DEFAULT_TSV_PATH,
                        help="Path to the TenK10K TSV file (can be a GCS path)")
    parser.add_argument("--show-progress-bar", action="store_true",
                        help="Show a progress bar")
    parser.add_argument("-o", "--output-catalog-json-path",
                        help="Path of the output catalog JSON file with TenK10K annotations")
    parser.add_argument("catalog_json_path", help="Path of the JSON catalog to annotate")
    args = parser.parse_args()

    # Download TSV if it's a remote path
    tsv_path = download_local_copy(args.tsv_path)
    if not os.path.isfile(tsv_path):
        parser.error(f"{args.tsv_path} file not found")

    if not os.path.isfile(args.catalog_json_path):
        parser.error(f"{args.catalog_json_path} file not found")

    if not args.output_catalog_json_path:
        args.output_catalog_json_path = args.catalog_json_path.replace(".json.gz", ".with_TenK10K_annotations.json.gz").replace(".json", ".with_TenK10K_annotations.json.gz")

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
        annotations["TenK10K_AlleleHistogram"] = row["allele_size_histogram"]
        biallelic_histogram = row.get("biallelic_histogram")
        if pd.notna(biallelic_histogram):
            annotations["TenK10K_BiallelicHistogram"] = biallelic_histogram

        # Integer fields
        if pd.notna(row["min_allele"]):
            annotations["TenK10K_MinAllele"] = int(row["min_allele"])
        if pd.notna(row["mode_allele"]):
            annotations["TenK10K_ModeAllele"] = int(row["mode_allele"])
        if pd.notna(row["max_allele"]):
            annotations["TenK10K_MaxAllele"] = int(row["max_allele"])
        if pd.notna(row["unique_allele_lengths"]):
            annotations["TenK10K_UniqueAlleleLengths"] = int(row["unique_allele_lengths"])
        if pd.notna(row["num_called_alleles"]):
            annotations["TenK10K_NumCalledAlleles"] = int(row["num_called_alleles"])

        # Float fields
        if pd.notna(row["stdev"]):
            annotations["TenK10K_Stdev"] = row["stdev"]
        if pd.notna(row["median"]):
            annotations["TenK10K_Median"] = row["median"]
        if pd.notna(row["99th_percentile"]):
            annotations["TenK10K_99thPercentile"] = row["99th_percentile"]

        # Computed rank fields
        if pd.notna(row["stdev_rank_by_motif"]):
            annotations["TenK10K_StdevRankByMotif"] = int(row["stdev_rank_by_motif"])
        if pd.notna(row["stdev_rank_total_number_by_motif"]):
            annotations["TenK10K_StdevRankTotalNumberByMotif"] = int(row["stdev_rank_total_number_by_motif"])

        annotation_lookup[locus_id] = annotations

    print(f"Built annotation lookup with {len(annotation_lookup):,d} entries")

    # Annotate the catalog
    input_locus_counter = annotated_locus_counter = 0
    print(f"Adding TenK10K annotations to {args.catalog_json_path}")
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
                if locus_id in annotation_lookup:
                    record.update(annotation_lookup[locus_id])
                    annotated_locus_counter += 1
                if i > 0:
                    f2.write(", ")
                f2.write(json.dumps(record, use_decimal=True, indent=4))
            f2.write("]")

    print(f"Annotated {annotated_locus_counter:,d} out of {input_locus_counter:,d} loci "
          f"({annotated_locus_counter/max(1, input_locus_counter):.1%})")
    print(f"Wrote annotated catalog to {args.output_catalog_json_path}")


if __name__ == "__main__":
    main()
