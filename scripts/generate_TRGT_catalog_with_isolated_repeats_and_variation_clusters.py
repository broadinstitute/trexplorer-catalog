"""This script takes a TSV file of variation clusters and a JSON file of all tandem repeats
and writes out a TRGT catalog with:
1. Variation clusters: loci grouped by vc_region, for every vc_region that either is shared by
   more than one locus, or was genuinely extended (non-zero offset) beyond its one locus's own
   original_region. A cluster's widest member spans the whole vc_region, so it can't be told apart
   from a genuinely solo locus by its own offset alone when the vc_region has other members too;
   that case is decided by how many loci land on the same vc_region instead
   (https://github.com/PacificBiosciences/trgt-lps/issues/5).
2. Isolated repeats: every other locus (a singleton vc_region with zero offset, i.e. identical to
   the locus's own original_region), using its own original_region coordinates.
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
    """Run a shell command, raising if it fails.

    os.system reports only the exit status of the last command in a pipeline, so `a | b` looks
    successful whenever b succeeds, however badly a failed. bgzip writes a valid empty file when its
    input is empty, so a missing bedtools would otherwise produce an empty catalog that the rest of
    the pipeline happily copies into the release. Running under pipefail makes the pipeline's status
    that of the first command to fail.
    """
    print(cmd)
    if os.system(f"set -o pipefail; {cmd}") != 0:
        raise RuntimeError(f"Command failed: {cmd}")


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
    vc_regions_with_nonzero_offset = set()
    locus_id_to_original_region_and_motifs = {}  # locus_id -> (original_region, motifs), for every non-filtered locus

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

            try:
                start_offset = float(vc_start_offset) if vc_start_offset else 0.0
                end_offset = float(vc_end_offset) if vc_end_offset else 0.0
            except ValueError:
                print(f"WARNING: Could not parse offsets for locus {locus_id}: "
                      f"start='{vc_start_offset}', end='{vc_end_offset}'")
                continue

            if not vc_region:
                print(f"WARNING: Empty vc_region for locus {locus_id}")
                continue

            vc_region_to_loci[vc_region].append((locus_id, motifs))
            locus_id_to_original_region_and_motifs[locus_id] = (original_region, motifs)
            if start_offset != 0 or end_offset != 0:
                vc_regions_with_nonzero_offset.add(vc_region)

    # Keep a vc_region unless it has exactly one member whose own offset is 0/0. A vc_region shared
    # by more than one locus is a real cluster even if its widest member's own offset is 0/0 (that
    # member's interval spans the whole vc_region, so it looks solo by offset alone). A vc_region
    # with only one locus is real too if that locus's own offset is non-zero: the region was
    # genuinely extended beyond the locus's original_region, even though nothing else was merged
    # in. Only a singleton vc_region with 0/0 offset is just that locus's own region and nothing
    # more.
    vc_region_to_loci = {
        region: loci for region, loci in vc_region_to_loci.items()
        if len(loci) > 1 or region in vc_regions_with_nonzero_offset
    }
    locus_ids_in_variation_clusters = {locus_id for loci in vc_region_to_loci.values() for locus_id, _ in loci}

    # Everything else that had a row in the TSV is an isolated repeat.
    isolated_repeats = []
    for locus_id, (original_region, motifs) in locus_id_to_original_region_and_motifs.items():
        if locus_id in locus_ids_in_variation_clusters:
            continue
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

        # A variation cluster gets an ID derived from its own coordinates (without the chr prefix)
        # rather than from the IDs of the repeats it contains. Previously a VC containing a single
        # repeat was given that repeat's ID, which made the two indistinguishable in TRGT output
        # (https://github.com/PacificBiosciences/trgt-lps/issues/5). The list of repeats that the VC
        # contains moves into STRUC, which TRGT copies through to the VCF unchanged.
        output_info = (f"ID=VC:{chrom.replace('chr', '')}:{start}-{end};"
                       f"MOTIFS={combined_motifs};STRUC=<VC:{combined_ids}>")
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

    restored_tr_counter = 0
    for record_i, record in enumerate(get_variant_catalog_iterator(
            args.input_repeat_catalog, show_progress_bar=args.show_progress_bar)):
        all_locus_ids_in_catalog.add(record["LocusId"])
        if record["LocusId"] not in all_locus_ids_in_tsv and record["LocusId"] not in locus_ids_filtered:
            locus_ids_missing_from_tsv.add(record["LocusId"])
            for output_row in convert_expansion_hunter_record_to_trgt_rows(record_i, record):
                restored_tr_counter += 1
                info_field_dict = parse_info_field(output_row[3])
                # The converter writes a TRGT structure expression such as "(AG)n", but the rest of
                # this catalog labels each row with <TR:{locus id}> or <VC:{...}>, which is what
                # trgt-lps and the downstream tooling read to tell an isolated repeat from a
                # variation cluster. Use the same label here so every restored row names its locus.
                info_field_dict["STRUC"] = f"<TR:{info_field_dict['ID']}>"
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
          f"{restored_tr_counter:,d} restored) covering {total_unique_locus_ids:,d} unique LocusIds to {args.output_bed_path}.gz")


if __name__ == "__main__":
    main()
