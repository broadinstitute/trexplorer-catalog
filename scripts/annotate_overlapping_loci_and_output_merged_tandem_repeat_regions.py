"""Annotate overlapping loci and output merged tandem repeat regions.

This script reads a JSON catalog of tandem repeat loci and:
1. Groups overlapping loci using group_overlapping_loci
2. Adds two annotation fields to each locus:
   - NonOverlappingPurestLocus (0 or 1): marks the locus with highest purity in each overlap group
   - NonOverlappingLongestLocus (0 or 1): marks the locus with longest interval in each overlap group
3. Outputs the annotated catalog as JSON
4. Outputs merged intervals as TRGT format BED and simple BED files
"""

import argparse
import collections
import gzip
import ijson
import os
import re
import simplejson as json
import tqdm

from str_analysis.utils.eh_catalog_utils import group_overlapping_loci, parse_motifs_from_locus_structure
from str_analysis.utils.misc_utils import parse_interval


def get_motif_from_record(record):
    """Extract the primary motif from a catalog record.

    Uses the Motif field if present, otherwise parses from LocusStructure.
    """
    if "Motif" in record:
        return record["Motif"]

    motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
    return motifs[0] if motifs else ""


def get_interval_length(record):
    """Get the length of the reference region interval."""
    chrom, start, end = parse_interval(record["ReferenceRegion"])
    return end - start


def get_sort_key(record):
    """Get a sort key for deterministic ordering: (chrom, start, end, motif)."""
    chrom, start, end = parse_interval(record["ReferenceRegion"])
    motif = get_motif_from_record(record)
    # Convert chrom to sortable value (extract number if possible)
    chrom_num = chrom.replace("chr", "")
    try:
        chrom_sort = (0, int(chrom_num))
    except ValueError:
        # Handle chrX, chrY, chrM, etc.
        chrom_sort = (1, chrom_num)
    return (chrom_sort, start, end, motif)


def select_purest_locus(loci):
    """Select the locus with highest purity, using tie-breakers.

    Selection criteria (in order):
    1. Highest ReferenceRepeatPurity
    2. Longest interval (end - start)
    3. Smallest motif size
    4. Smallest genomic coordinates (chrom, start, end, motif)
    """
    def sort_key(record):
        purity = record.get("ReferenceRepeatPurity", 0)
        length = get_interval_length(record)
        motif = get_motif_from_record(record)
        motif_size = len(motif)
        coords = get_sort_key(record)
        # Negate purity and length for descending sort, keep motif_size and coords ascending
        return (-purity, -length, motif_size, coords)

    return min(loci, key=sort_key)


def select_longest_locus(loci):
    """Select the locus with longest interval, using tie-breakers.

    Selection criteria (in order):
    1. Longest interval (end - start)
    2. Highest ReferenceRepeatPurity
    3. Smallest motif size
    4. Smallest genomic coordinates (chrom, start, end, motif)
    """
    def sort_key(record):
        length = get_interval_length(record)
        purity = record.get("ReferenceRepeatPurity", 0)
        motif = get_motif_from_record(record)
        motif_size = len(motif)
        coords = get_sort_key(record)
        # Negate length and purity for descending sort, keep motif_size and coords ascending
        return (-length, -purity, motif_size, coords)

    return min(loci, key=sort_key)


def compute_merged_interval(loci):
    """Compute the union interval (min start, max end) for a group of loci."""
    min_start = float('inf')
    max_end = 0
    chrom = None

    for locus in loci:
        c, start, end = parse_interval(locus["ReferenceRegion"])
        if chrom is None:
            chrom = c
        min_start = min(min_start, start)
        max_end = max(max_end, end)

    return chrom, min_start, max_end


def get_sorted_locus_ids(loci):
    """Get locus IDs sorted by (chrom, start, end, motif)."""
    sorted_loci = sorted(loci, key=get_sort_key)
    return [locus["LocusId"] for locus in sorted_loci]


def get_unique_motifs(loci):
    """Get unique motifs from all loci, preserving order of first occurrence."""
    seen = set()
    motifs = []
    for locus in sorted(loci, key=get_sort_key):
        motif = get_motif_from_record(locus)
        if motif and motif not in seen:
            seen.add(motif)
            motifs.append(motif)
    return motifs


def main():
    parser = argparse.ArgumentParser(
        description="Annotate overlapping loci and output merged tandem repeat regions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("catalog_json_path", help="Path of the input JSON catalog (gzipped)")
    parser.add_argument("-o", "--output-catalog-json-path",
                        help="Path of the output annotated catalog JSON file")
    parser.add_argument("--output-trgt-bed-path",
                        help="Path of the output TRGT format BED file with merged intervals")
    parser.add_argument("--output-simple-bed-path",
                        help="Path of the output simple BED file with merged intervals")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    args = parser.parse_args()

    if not os.path.isfile(args.catalog_json_path):
        parser.error(f"File not found: {args.catalog_json_path}")

    # Set default output paths
    if not args.output_catalog_json_path:
        args.output_catalog_json_path = args.catalog_json_path.replace(
            ".json.gz", ".with_non_overlapping_annotations.json.gz"
        ).replace(".json", ".with_non_overlapping_annotations.json")
        if not args.output_catalog_json_path.endswith(".gz"):
            args.output_catalog_json_path += ".gz"

    if not args.output_trgt_bed_path:
        args.output_trgt_bed_path = args.output_catalog_json_path.replace(
            ".json.gz", ".merged_regions.trgt.bed"
        ).replace(".json", ".merged_regions.trgt.bed")

    if not args.output_simple_bed_path:
        args.output_simple_bed_path = args.output_catalog_json_path.replace(
            ".json.gz", ".merged_regions.bed"
        ).replace(".json", ".merged_regions.bed")

    # First pass: group overlapping loci and determine winners
    print(f"First pass: grouping overlapping loci from {args.catalog_json_path}")

    # Store which locus IDs are winners for each annotation
    purest_locus_ids = set()
    longest_locus_ids = set()

    # Store merged interval info for BED output
    # Each entry: (chrom, start, end, locus_ids, motifs)
    merged_intervals = []

    # Counters
    total_loci = 0
    single_locus_groups = 0
    multi_locus_groups = 0
    max_group_size = 0
    group_size_histogram = collections.Counter()

    fopen = gzip.open if args.catalog_json_path.endswith("gz") else open
    with fopen(args.catalog_json_path, "rt") as f:
        iterator = ijson.items(f, "item")
        if args.show_progress_bar:
            iterator = tqdm.tqdm(iterator, unit=" records", unit_scale=True, desc="Pass 1")

        for group in group_overlapping_loci(iterator, only_group_loci_with_similar_motifs=False):
            group_size = len(group)
            total_loci += group_size
            group_size_histogram[group_size] += 1
            max_group_size = max(max_group_size, group_size)

            if group_size == 1:
                # Single locus: mark as both purest and longest
                single_locus_groups += 1
                locus = group[0]
                purest_locus_ids.add(locus["LocusId"])
                longest_locus_ids.add(locus["LocusId"])

                # Store interval info
                chrom, start, end = parse_interval(locus["ReferenceRegion"])
                motif = get_motif_from_record(locus)
                merged_intervals.append((chrom, start, end, [locus["LocusId"]], [motif] if motif else []))
            else:
                # Multiple overlapping loci
                multi_locus_groups += 1

                # Select winners
                purest = select_purest_locus(group)
                longest = select_longest_locus(group)
                purest_locus_ids.add(purest["LocusId"])
                longest_locus_ids.add(longest["LocusId"])

                # Compute merged interval
                chrom, start, end = compute_merged_interval(group)
                locus_ids = get_sorted_locus_ids(group)
                motifs = get_unique_motifs(group)
                merged_intervals.append((chrom, start, end, locus_ids, motifs))

    print(f"Found {total_loci:,d} total loci")
    print(f"  - {single_locus_groups:,d} isolated loci (no overlaps)")
    print(f"  - {multi_locus_groups:,d} groups of overlapping loci")
    print(f"  - Largest overlap group: {max_group_size} loci")

    if args.verbose:
        print("Group size histogram:")
        for size in sorted(group_size_histogram.keys()):
            count = group_size_histogram[size]
            print(f"  {size} loci: {count:,d} groups")

    # Second pass: annotate catalog and write output
    print(f"Second pass: annotating catalog and writing to {args.output_catalog_json_path}")

    fopen = gzip.open if args.catalog_json_path.endswith("gz") else open
    with fopen(args.catalog_json_path, "rt") as f:
        f2open = gzip.open if args.output_catalog_json_path.endswith("gz") else open
        with f2open(args.output_catalog_json_path, "wt") as f2:
            iterator = ijson.items(f, "item")
            if args.show_progress_bar:
                iterator = tqdm.tqdm(iterator, unit=" records", unit_scale=True, desc="Pass 2")

            annotated_purest_count = 0
            annotated_longest_count = 0

            f2.write("[")
            for i, record in enumerate(iterator):
                locus_id = record["LocusId"]

                # Add annotations
                is_purest = 1 if locus_id in purest_locus_ids else 0
                is_longest = 1 if locus_id in longest_locus_ids else 0
                record["NonOverlappingPurestLocus"] = is_purest
                record["NonOverlappingLongestLocus"] = is_longest

                if is_purest:
                    annotated_purest_count += 1
                if is_longest:
                    annotated_longest_count += 1

                if i > 0:
                    f2.write(", ")
                f2.write(json.dumps(record, use_decimal=True, indent=4))

            f2.write("]")

    print(f"Annotated {total_loci:,d} loci")
    print(f"  - {annotated_purest_count:,d} loci marked as NonOverlappingPurestLocus=1")
    print(f"  - {annotated_longest_count:,d} loci marked as NonOverlappingLongestLocus=1")
    print(f"Wrote annotated catalog to {args.output_catalog_json_path}")

    # Write TRGT format BED file
    print(f"Writing TRGT format BED to {args.output_trgt_bed_path}")
    with open(args.output_trgt_bed_path, "w") as f:
        for chrom, start, end, locus_ids, motifs in merged_intervals:
            locus_id_str = ",".join(locus_ids)
            motifs_str = ",".join(motifs) if motifs else ""
            # Strip chr prefix for STRUC field (matching convention in other scripts)
            chrom_for_struc = chrom.replace("chr", "")
            struc = f"<MR:{chrom_for_struc}-{start}-{end}>"
            f.write(f"{chrom}\t{start}\t{end}\tID={locus_id_str};MOTIFS={motifs_str};STRUC={struc}\n")

    print(f"Wrote {len(merged_intervals):,d} merged intervals to {args.output_trgt_bed_path}")

    # Write simple BED file
    print(f"Writing simple BED to {args.output_simple_bed_path}")
    with open(args.output_simple_bed_path, "w") as f:
        for chrom, start, end, locus_ids, motifs in merged_intervals:
            locus_id_str = ",".join(locus_ids)
            f.write(f"{chrom}\t{start}\t{end}\t{locus_id_str}\n")

    print(f"Wrote {len(merged_intervals):,d} merged intervals to {args.output_simple_bed_path}")

    # Sort and compress BED files
    print("Sorting and compressing BED files...")
    for bed_path in [args.output_trgt_bed_path, args.output_simple_bed_path]:
        os.system(f"bedtools sort -i {bed_path} > {bed_path}.sorted && mv {bed_path}.sorted {bed_path}")
        os.system(f"bgzip -f {bed_path}")
        print(f"  Sorted and compressed: {bed_path}.gz")

    args.output_trgt_bed_path += ".gz"
    args.output_simple_bed_path += ".gz"

    # Final statistics
    print("\n=== Summary ===")
    print(f"Input catalog: {args.catalog_json_path}")
    print(f"Total loci: {total_loci:,d}")
    print(f"Overlap groups: {single_locus_groups + multi_locus_groups:,d}")
    print(f"  - Isolated loci (no overlaps): {single_locus_groups:,d}")
    print(f"  - Groups with overlaps: {multi_locus_groups:,d}")
    print(f"  - Largest group size: {max_group_size}")
    print(f"Output files:")
    print(f"  - Annotated catalog: {args.output_catalog_json_path}")
    print(f"  - TRGT format BED: {args.output_trgt_bed_path}")
    print(f"  - Simple BED: {args.output_simple_bed_path}")


if __name__ == "__main__":
    main()
