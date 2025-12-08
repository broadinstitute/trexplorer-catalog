"""This script takes a BED file of variation clusters and isolated TRs, as well as a JSON file of all tandem repeats
and writes out a TRGT catalog with all input variation clusters as well as any tandem repeats from the input catalog
that aren't embedded in variation clusters (ie. are isolated repeats).
"""


import argparse
import collections
import gzip
import os
import simplejson as json
import re
import tqdm

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator
from str_analysis.convert_expansion_hunter_catalog_to_trgt_catalog import convert_expansion_hunter_record_to_trgt_rows

MINIMUM_CHANGE_TO_BOUNDARIES_THRESHOLD = 6

def run(cmd):
    print(cmd)
    os.system(cmd)


def parse_info_field(info_field):
    """Parse a TRGT catalog info field into a python dictionary"""
    result = {}
    for key_value in info_field.split(";"):
        key_value = key_value.split("=")
        if len(key_value) != 2:
            raise ValueError(f"Invalid key-value pair '{key_value}' in line {info_field}")
        key, value = key_value
        result[key] = value
    return result


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("-o", "--output-bed-path", help="Path of output BED file.")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    parser.add_argument("input_isolated_repeats_and_variation_clusters_bed_path", 
                        help="Path of the input variation clusters BED file")
    parser.add_argument("input_repeat_catalog", help="Catalog of all tandem repeats in JSON or BED format")
    args = parser.parse_args()

    if ".bed" not in args.input_isolated_repeats_and_variation_clusters_bed_path:
        parser.error("--input-variation-clusters-bed-path must have a '.bed' suffix")

    if not args.output_bed_path:
        args.output_bed_path = re.sub(".bed(.b?gz)?$", "", args.input_isolated_repeats_and_variation_clusters_bed_path)
        args.output_bed_path += ".processed.bed"
    elif args.output_bed_path.endswith(".bed.gz"):
        args.output_bed_path = re.sub(".gz$", "", args.output_bed_path)
    elif not args.output_bed_path.endswith(".bed"):
        parser.error("--output-bed-path must have a '.bed' suffix")

    # copy all lines from input_isolated_repeats_and_variation_clusters_bed_path to the output file while
    # parsing the LocusIds of all TRs
    locus_ids_in_variation_clusters = set()
    locus_ids_isolated_repeats = set()
    output_bed_file = open(args.output_bed_path, "wt")
    output_row_counter = 0
    fopen = gzip.open if args.input_isolated_repeats_and_variation_clusters_bed_path.endswith("gz") else open
    with fopen(args.input_isolated_repeats_and_variation_clusters_bed_path, "rt") as f:
        if args.show_progress_bar:
            f = tqdm.tqdm(f, unit=" records", unit_scale=True)

        for line in f:
            fields = line.strip("\n").split("\t")

            info_field_dict = parse_info_field(fields[3])

            # Parse LocusIds from the 'ID' field
            for locus_id in info_field_dict["ID"].split(","):
                if locus_id.count("-") != 3:
                    raise ValueError(f"Unexpected locus_id '{locus_id}'")
                if locus_id in locus_ids_in_variation_clusters:
                    raise ValueError(f"locus_id '{locus_id}' occurs more than once")

                if "VC" in info_field_dict["STRUC"]:
                    locus_ids_in_variation_clusters.add(locus_id)
                elif "TR" in info_field_dict["STRUC"]:
                    locus_ids_isolated_repeats.add(locus_id)
                else:
                    raise ValueError(f"STRUC does not contain either 'VC' or 'TR' in line {fields}")

            output_bed_file.write(line)
            output_row_counter += 1

    all_locus_ids_in_variation_clusters_file = locus_ids_in_variation_clusters | locus_ids_isolated_repeats
    print(f"Parsed {len(locus_ids_in_variation_clusters):,d} variation clusters "
          f"and {len(locus_ids_isolated_repeats):,d} isolated repeats from "
          f"{args.input_isolated_repeats_and_variation_clusters_bed_path}")

    locus_ids_missing_from_variation_clusters_file = set()
    all_locus_ids_in_catalog = set()

    filtered_tr_counter = 0
    for record_i, record in enumerate(get_variant_catalog_iterator(
            args.input_repeat_catalog, show_progress_bar=args.show_progress_bar)):
        all_locus_ids_in_catalog.add(record["LocusId"])
        if record["LocusId"] not in all_locus_ids_in_variation_clusters_file:
            locus_ids_missing_from_variation_clusters_file.add(record["LocusId"])
            for output_row in convert_expansion_hunter_record_to_trgt_rows(record_i, record):
                filtered_tr_counter += 1
                info_field_dict = parse_info_field(output_row[3])
                info_field_dict["STRUC"] = f"<TR:FILTERED{filtered_tr_counter}>"
                output_row[3] = ";".join(f"{key}={value}" for key, value in info_field_dict.items())
                output_bed_file.write("\t".join(map(str, output_row)) + "\n")
                output_row_counter += 1


    print(f"Restored {len(locus_ids_missing_from_variation_clusters_file):,d} TRs that were not in the variation "
          f"clusters file")
    output_bed_file.close()

    unexpected_locus_ids_in_variation_cluster_catalog = all_locus_ids_in_variation_clusters_file - all_locus_ids_in_catalog
    if unexpected_locus_ids_in_variation_cluster_catalog:
        raise ValueError(f"{len(unexpected_locus_ids_in_variation_cluster_catalog)} locus IDs in the variation cluster "
            f"catalog were not found in the input repeat catalog: {unexpected_locus_ids_in_variation_cluster_catalog}")

    run(f"bedtools sort -i {args.output_bed_path} | bgzip > {args.output_bed_path}.sorted")
    run(f"mv {args.output_bed_path}.sorted {args.output_bed_path}.gz")
    os.remove(args.output_bed_path)

    print(f"Wrote {output_row_counter:,d} rows and "
          f"{len(all_locus_ids_in_variation_clusters_file | locus_ids_missing_from_variation_clusters_file):,d} "
          f"unique LocusIds to {args.output_bed_path}.gz")


if __name__ == "__main__":
    main()

