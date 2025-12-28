"""This script takes a TRGT catalog BED file and outputs a new BED file in LongTR format.
It expects the TRGT locus ID to have the format "chr1-123456-123467-CAG" for isolated repeats, and
"chr1-123456-123467-CAG,chr1-123467-123489-CTG" for compound definitions that contain multiple adjacent TRs.
"""

import argparse
import collections
import gzip
import os
import re
import tqdm


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


def compute_dominant_motif(info_field_dict):
    """Compute the dominant motif for a TR locus. For compound definitions (those that span multiple
    adjacent TRs), this is the motif of the constituent TR that spans the largest interval (for all other loci).

    Args:
        info_field_dict (dict): A dictionary of info fields for a TR locus

    Return:
        str: The dominant motif for the TR locus.
    """
    motifs = []
    for locus_id in info_field_dict["ID"].split(","):
        if locus_id.count("-") != 3:
            raise ValueError(f"Unexpected locus_id '{locus_id}'")

        chrom, start_0based, end, motif = locus_id.split("-")
        motifs.append((int(end) - int(start_0based), motif))
        
    _, dominant_motif = max(motifs)
    return dominant_motif


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("-o", "--output-bed-path", help="Path of output BED file.")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    parser.add_argument("input_trgt_catalog_bed_path", help="Path of the input TRGT catalog BED file")
    args = parser.parse_args()

    if ".bed" not in args.input_trgt_catalog_bed_path:
        parser.error(f"Input path {args.input_trgt_catalog_bed_path} must have a '.bed' suffix")

    if not args.output_bed_path:
        args.output_bed_path = re.sub(".bed(.b?gz)?$", "", args.input_trgt_catalog_bed_path)
        args.output_bed_path = args.output_bed_path.replace(".TRGT", "")
        args.output_bed_path += ".LongTR.bed"
    elif not args.output_bed_path.endswith(".bed"):
        parser.error("--output-bed-path must have a '.bed' suffix")

    counter = collections.Counter()
    output_bed_file = open(args.output_bed_path, "wt")
    fopen = gzip.open if args.input_trgt_catalog_bed_path.endswith("gz") else open
    with fopen(args.input_trgt_catalog_bed_path, "rt") as f:
        if args.show_progress_bar:
            f = tqdm.tqdm(f, unit=" records", unit_scale=True)

        for line in f:
            counter["total"] += 1
            fields = line.strip("\n").split("\t")
            chrom = fields[0]
            start_0based = int(fields[1])
            end_1based = int(fields[2])

            if start_0based + 1 >= end_1based:
                # avoid "Region has a STOP <= START" error
                counter["skipped"] += 1
                if args.verbose:
                    print(f"WARNING: Skipping record #{counter['skipped']} because the interval has width {end_1based - start_0based}bp")
                continue

            info_field_dict = parse_info_field(fields[3])
            motif = compute_dominant_motif(info_field_dict)

            counter["output"] += 1
            output_bed_file.write("\t".join(map(str, [
                chrom,
                start_0based + 1,  # LongTR BED files use 1-based coords
                end_1based,
                motif,
                info_field_dict["ID"],
            ])) + "\n")

    os.system(f"bgzip -f {args.output_bed_path}")

    print(f"Wrote {counter['output']:,d} out of {counter['total']:,d} rows to {args.output_bed_path}.gz")


if __name__ == "__main__":
    main()

