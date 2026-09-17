"""Add AoU1027 population statistics annotations to a TR catalog JSON file.

AoU1027 is a cohort of 1027 long-read sequencing samples from the All of Us Research Program
used to compute allele frequency statistics at tandem repeat loci. This script annotates
catalog records with these population-level statistics.

This script owns every AoU1027_* field: each run clears them first and then writes only the
ones that apply, so running it again over an already-annotated catalog replaces the annotations
rather than layering onto them.
"""

import argparse
import gzip
import ijson
import os
import simplejson as json
import tqdm

from str_analysis.utils.file_utils import download_local_copy
from catalog_annotation_utils import (clear_previous_annotations, find_unlisted_annotation_fields,
                                      print_annotation_replacement_summary)

DEFAULT_TSV_PATH = "gs://tandem-repeat-catalog/v2.0/AoULR_phase1_TRGT_Weisburd_v1_combined.txt.gz"

# Every field this script writes. A locus that drops out of the TSV, or whose row no longer has a
# value for one of the optional fields, must lose whatever a previous run gave it, so these are
# cleared per record before the new ones are written.
AOU1027_FIELDS = (
    "AoU1027_MinAllele",
    "AoU1027_ModeAllele",
    "AoU1027_MaxAllele",
    "AoU1027_Stdev",
    "AoU1027_Median",
    "AoU1027_99thPercentile",
    "AoU1027_NumCalledAlleles",
    "AoU1027_StdevRankByMotif",
    "AoU1027_StdevRankTotalNumberByMotif",
    "AoU1027_OE_Length",
    "AoU1027_OE_LengthPercentile",
)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Add AoU1027 population statistics annotations to a TR catalog JSON file."
    )
    parser.add_argument("--tsv-path", default=DEFAULT_TSV_PATH,
                        help="Path to the AoU1027 TSV file (can be a GCS path)")
    parser.add_argument("--show-progress-bar", action="store_true",
                        help="Show a progress bar")
    parser.add_argument("-o", "--output-catalog-json-path",
                        help="Path of the output catalog JSON file with AoU1027 annotations")
    parser.add_argument("catalog_json_path", help="Path of the JSON catalog to annotate")
    args = parser.parse_args()

    # Download TSV if it's a remote path
    tsv_path = download_local_copy(args.tsv_path)
    if not os.path.isfile(tsv_path):
        parser.error(f"{args.tsv_path} file not found")

    if not os.path.isfile(args.catalog_json_path):
        parser.error(f"{args.catalog_json_path} file not found")

    if not args.output_catalog_json_path:
        args.output_catalog_json_path = args.catalog_json_path.replace(".json.gz", ".with_AoU1027_annotations.json.gz").replace(".json", ".with_AoU1027_annotations.json.gz")

    # Parse the AoU1027 TSV file
    print(f"Parsing {args.tsv_path}")
    annotation_lookup = {}

    fopen = gzip.open if tsv_path.endswith("gz") else open
    with fopen(tsv_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        col_indices = {col: i for i, col in enumerate(header)}

        expected_columns = {"TRID2", "longestPureSegmentMotif", "0thPercentile", "Mode", "Stdev",
                           "50thPercentile", "99thPercentile", "100thPercentile", "numCalledAlleles",
                           "StdevRankByMotif", "StdevRankTotalNumberByMotif", "OE_len", "OE_len_percentile"}
        missing_columns = expected_columns - set(col_indices.keys())
        if missing_columns:
            parser.error(f"{args.tsv_path} is missing expected columns: {missing_columns}")

        line_iterator = f
        if args.show_progress_bar:
            line_iterator = tqdm.tqdm(f, unit=" lines", unit_scale=True,
                                      desc="Building annotation lookup")

        for line in line_iterator:
            fields = line.rstrip("\n").split("\t")
            locus_id = fields[col_indices["TRID2"]]

            # annotation_lookup below is keyed by locus_id alone, so a duplicate would silently
            # overwrite an earlier row's annotations rather than raising. This is the failure mode
            # of https://github.com/PacificBiosciences/trgt-lps/issues/5, just one step downstream.
            if locus_id in annotation_lookup:
                parser.error(f"{args.tsv_path} has duplicate TRID2 value {locus_id}")

            motif_size = len(fields[col_indices["longestPureSegmentMotif"]])

            annotations = {}

            # Integer fields (converted from bp to repeat units)
            annotations["AoU1027_MinAllele"] = int(float(fields[col_indices["0thPercentile"]])) // motif_size
            annotations["AoU1027_ModeAllele"] = int(float(fields[col_indices["Mode"]])) // motif_size
            annotations["AoU1027_MaxAllele"] = int(float(fields[col_indices["100thPercentile"]])) // motif_size

            # Float fields (converted from bp to repeat units)
            annotations["AoU1027_Stdev"] = float(fields[col_indices["Stdev"]]) / motif_size
            annotations["AoU1027_Median"] = float(fields[col_indices["50thPercentile"]]) // motif_size
            annotations["AoU1027_99thPercentile"] = float(fields[col_indices["99thPercentile"]]) // motif_size

            # Integer fields (no unit conversion)
            annotations["AoU1027_NumCalledAlleles"] = int(fields[col_indices["numCalledAlleles"]])
            annotations["AoU1027_StdevRankByMotif"] = int(fields[col_indices["StdevRankByMotif"]])
            annotations["AoU1027_StdevRankTotalNumberByMotif"] = int(fields[col_indices["StdevRankTotalNumberByMotif"]])

            # Float fields (can be empty string)
            oe_len_value = fields[col_indices["OE_len"]]
            if oe_len_value != "":
                annotations["AoU1027_OE_Length"] = float(oe_len_value)

            oe_len_percentile_value = fields[col_indices["OE_len_percentile"]]
            if oe_len_percentile_value != "":
                annotations["AoU1027_OE_LengthPercentile"] = float(oe_len_percentile_value)

            annotation_lookup[locus_id] = annotations

    print(f"Built annotation lookup with {len(annotation_lookup):,d} entries")

    unlisted_fields = find_unlisted_annotation_fields(annotation_lookup, AOU1027_FIELDS)
    if unlisted_fields:
        parser.error(f"AOU1027_FIELDS doesn't list {sorted(unlisted_fields)}, so re-running this script "
                     f"over an already-annotated catalog would leave those fields stale")

    # Annotate the catalog
    input_locus_counter = annotated_locus_counter = 0
    locus_with_previous_annotation_counter = locus_that_lost_annotation_counter = 0
    print(f"Adding AoU1027 annotations to {args.catalog_json_path}")
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

                had_previous_annotation = clear_previous_annotations(record, AOU1027_FIELDS)
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
    print_annotation_replacement_summary("AoU1027", f"are no longer in {args.tsv_path}",
                                         locus_with_previous_annotation_counter,
                                         locus_that_lost_annotation_counter)
    print(f"Wrote annotated catalog to {args.output_catalog_json_path}")


if __name__ == "__main__":
    main()
