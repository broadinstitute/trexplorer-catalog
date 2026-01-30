"""Add AoU1027 population statistics annotations to a TR catalog JSON file.

AoU1027 is a cohort of 1027 long-read sequencing samples from the All of Us Research Program
used to compute allele frequency statistics at tandem repeat loci. This script annotates
catalog records with these population-level statistics.
"""

import argparse
import gzip
import ijson
import os
import simplejson as json
import tqdm

from str_analysis.utils.file_utils import download_local_copy

DEFAULT_TSV_PATH = "gs://tandem-repeat-catalog/v2.0/AoULR_phase1_TRGT_Weisburd_v1_combined.txt.gz"


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

    # Annotate the catalog
    input_locus_counter = annotated_locus_counter = 0
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
