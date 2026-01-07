"""This script takes a TSV file of variation clusters and a JSON file of all tandem repeats
and writes out a TRGT catalog with:
1. Variation clusters (loci with non-zero offsets, grouped by vc_region)
2. Isolated repeats (loci with zero offsets, using original_region coordinates)
3. Any tandem repeats from the catalog that aren't in the TSV (restored from catalog)

Loci filtered due to DEPTH or EXTENSION are excluded entirely.
"""


import argparse
import collections
import gzip
import os
import re
import tqdm

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator
from str_analysis.convert_expansion_hunter_catalog_to_trgt_catalog import convert_expansion_hunter_record_to_trgt_rows


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
    parser.add_argument("input_variation_clusters_tsv_path",
                        help="Path of the input variation clusters TSV file")
    parser.add_argument("input_repeat_catalog", help="Catalog of all tandem repeats in JSON or BED format")
    args = parser.parse_args()

    if not args.output_bed_path:
        args.output_bed_path = re.sub(r"\.tsv(\.gz)?$", "", args.input_variation_clusters_tsv_path)
        args.output_bed_path += ".TRGT.bed"
    elif args.output_bed_path.endswith(".bed.gz"):
        args.output_bed_path = re.sub(r"\.gz$", "", args.output_bed_path)
    elif not args.output_bed_path.endswith(".bed"):
        parser.error("--output-bed-path must have a '.bed' suffix")

    # Data structures to collect variation clusters (grouped by vc_region)
    # and isolated repeats
    vc_region_to_loci = collections.defaultdict(list)  # vc_region -> list of (locus_id, motifs)
    isolated_repeats = []  # list of (chrom, start, end, locus_id, motifs)

    locus_ids_in_variation_clusters = set()
    locus_ids_isolated_repeats = set()
    locus_ids_filtered = set()

    fopen = gzip.open if args.input_variation_clusters_tsv_path.endswith("gz") else open
    with fopen(args.input_variation_clusters_tsv_path, "rt") as f:
        if args.show_progress_bar:
            f = tqdm.tqdm(f, unit=" records", unit_scale=True)

        header = None
        for line in f:
            fields = line.strip("\n").split("\t")

            # Parse header
            if header is None:
                header = fields
                continue

            # Parse TSV columns
            region_info = fields[0]
            original_region = fields[1]
            vc_start_offset = fields[2]
            vc_end_offset = fields[3]
            vc_region = fields[4] if len(fields) > 4 else ""

            # Extract info from region_info
            info_dict = parse_info_field(region_info)
            locus_id = info_dict["ID"]
            motifs = info_dict.get("MOTIFS", "")

            # Skip filtered loci
            if vc_end_offset in ("DEPTH", "EXTENSION"):
                locus_ids_filtered.add(locus_id)
                continue

            # Parse offsets
            try:
                start_offset = float(vc_start_offset) if vc_start_offset else 0.0
                end_offset = float(vc_end_offset) if vc_end_offset else 0.0
            except ValueError:
                print(f"WARNING: Could not parse offsets for locus {locus_id}: start='{vc_start_offset}', end='{vc_end_offset}'")
                continue

            # Determine if this is a VC or isolated TR
            if start_offset != 0 or end_offset != 0:
                # Variation cluster - group by vc_region
                if not vc_region:
                    print(f"WARNING: Non-zero offsets but empty vc_region for locus {locus_id}")
                    continue

                vc_region_to_loci[vc_region].append((locus_id, motifs))
                locus_ids_in_variation_clusters.add(locus_id)
            else:
                # Isolated repeat - store for later output
                chrom, start, end = parse_interval(original_region)
                isolated_repeats.append((chrom, start, end, locus_id, motifs))
                locus_ids_isolated_repeats.add(locus_id)

    print(f"Parsed TSV: {len(vc_region_to_loci):,d} variation clusters containing {len(locus_ids_in_variation_clusters):,d} loci, "
          f"{len(locus_ids_isolated_repeats):,d} isolated repeats, "
          f"{len(locus_ids_filtered):,d} filtered (excluded)")

    # Write output BED file
    output_bed_file = open(args.output_bed_path, "wt")
    output_row_counter = 0

    # Write variation clusters (one row per VC, with combined IDs and MOTIFs)
    vc_counter = 0
    for vc_region, loci in vc_region_to_loci.items():
        vc_counter += 1
        chrom, start, end = parse_interval(vc_region)

        # Combine IDs and MOTIFs from all loci in this VC
        combined_ids = ",".join(locus_id for locus_id, _ in loci)
        # Get unique motifs while preserving order
        seen_motifs = set()
        unique_motifs = []
        for _, motifs in loci:
            for m in motifs.split(","):
                if m and m not in seen_motifs:
                    seen_motifs.add(m)
                    unique_motifs.append(m)
        combined_motifs = ",".join(unique_motifs)

        # STRUC uses VC coordinates (without chr prefix)
        chrom_for_struc = chrom.replace("chr", "")
        struc = f"<VC:{chrom_for_struc}-{start}-{end}>"

        output_info = f"ID={combined_ids};MOTIFS={combined_motifs};STRUC={struc}"
        output_bed_file.write(f"{chrom}\t{start}\t{end}\t{output_info}\n")
        output_row_counter += 1

    # Write isolated repeats (one row per TR)
    for chrom, start, end, locus_id, motifs in isolated_repeats:
        struc = f"<TR:{locus_id}>"
        output_info = f"ID={locus_id};MOTIFS={motifs};STRUC={struc}"
        output_bed_file.write(f"{chrom}\t{start}\t{end}\t{output_info}\n")
        output_row_counter += 1

    all_locus_ids_in_tsv = locus_ids_in_variation_clusters | locus_ids_isolated_repeats

    # Restore any TRs from the catalog that weren't in the TSV
    locus_ids_missing_from_tsv = set()
    all_locus_ids_in_catalog = set()

    filtered_tr_counter = 0
    for record_i, record in enumerate(get_variant_catalog_iterator(
            args.input_repeat_catalog, show_progress_bar=args.show_progress_bar)):
        all_locus_ids_in_catalog.add(record["LocusId"])
        if record["LocusId"] not in all_locus_ids_in_tsv and record["LocusId"] not in locus_ids_filtered:
            locus_ids_missing_from_tsv.add(record["LocusId"])
            for output_row in convert_expansion_hunter_record_to_trgt_rows(record_i, record):
                filtered_tr_counter += 1
                info_field_dict = parse_info_field(output_row[3])
                info_field_dict["STRUC"] = f"<TR:FILTERED{filtered_tr_counter}>"
                output_row[3] = ";".join(f"{key}={value}" for key, value in info_field_dict.items())
                output_bed_file.write("\t".join(map(str, output_row)) + "\n")
                output_row_counter += 1

    print(f"Restored {len(locus_ids_missing_from_tsv):,d} TRs that were not in the variation clusters TSV")
    output_bed_file.close()

    # Check for unexpected locus IDs
    unexpected_locus_ids = (all_locus_ids_in_tsv | locus_ids_filtered) - all_locus_ids_in_catalog
    if unexpected_locus_ids:
        print(f"WARNING: {len(unexpected_locus_ids):,d} locus IDs in the TSV were not found in the input repeat catalog. "
              f"Examples: {list(unexpected_locus_ids)[:10]}")

    # Sort and compress output
    run(f"bedtools sort -i {args.output_bed_path} | bgzip > {args.output_bed_path}.sorted")
    run(f"mv {args.output_bed_path}.sorted {args.output_bed_path}.gz")
    os.remove(args.output_bed_path)

    total_unique_locus_ids = len(all_locus_ids_in_tsv | locus_ids_missing_from_tsv)
    print(f"Wrote {output_row_counter:,d} rows ({vc_counter:,d} VCs + {len(isolated_repeats):,d} isolated TRs + "
          f"{filtered_tr_counter:,d} restored) covering {total_unique_locus_ids:,d} unique LocusIds to {args.output_bed_path}.gz")


if __name__ == "__main__":
    main()
