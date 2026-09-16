"""This script takes a TSV file of variation clusters and a JSON file of all tandem repeats
and writes out a TRGT catalog with:
1. One row for every repeat in the input catalog, carrying the same TRGT structure expression the
   repeat catalog itself uses. Repeats the TSV marked DEPTH or EXTENSION are flagged with a
   VC_FILTER key naming the reason rather than being dropped.
2. One additional row for every variation cluster, listing the repeats it contains in its STRUC field.

A repeat that belongs to a cluster therefore appears twice, once on its own and once as a member of
the cluster. The two rows usually cover different spans and are meant to coexist, so nothing here
chooses between them. A minority of repeat rows turn out to span exactly what some cluster row spans;
the counter below reports how many, matching on span alone, so such a pair is not necessarily a
cluster and one of its own members. Those rows are reported rather than dropped, since the catalog's
contract is one row per repeat plus one row per cluster. They can differ in MOTIFS as well as in ID
and STRUC, because a cluster's MOTIFS is the union over its members while a repeat's is its own.

Some repeats named by a cluster have no row of their own because they are not in the repeat catalog.
That is expected: the cluster table was computed against a different set of ids.
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

    vc_region_to_loci = collections.defaultdict(list)  # vc_region -> list of (locus_id, motifs)
    vc_regions_with_nonzero_offset = set()
    filter_reason_by_locus_id = {}  # locus_id -> "DEPTH" or "EXTENSION"

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

            region_info, _, vc_start_offset, vc_end_offset = fields[:4]
            vc_region = fields[4] if len(fields) > 4 else ""

            info_dict = parse_info_field(region_info)
            locus_id = info_dict["ID"]

            # Loci the cluster table could not process. They keep their row in the output and carry
            # the reason, so a caller can filter on it instead of finding the locus simply absent.
            if vc_end_offset in ("DEPTH", "EXTENSION"):
                filter_reason_by_locus_id[locus_id] = vc_end_offset
                continue

            try:
                start_offset = float(vc_start_offset) if vc_start_offset else 0.0
                end_offset = float(vc_end_offset) if vc_end_offset else 0.0
            except ValueError:
                print(f"WARNING: Could not parse offsets for locus {locus_id}: "
                      f"start='{vc_start_offset}', end='{vc_end_offset}'")
                continue

            if not vc_region:
                continue

            vc_region_to_loci[vc_region].append((locus_id, info_dict.get("MOTIFS", "")))
            if start_offset != 0 or end_offset != 0:
                vc_regions_with_nonzero_offset.add(vc_region)

    # Keep a vc_region unless it has exactly one member whose own offset is 0/0. A vc_region shared
    # by more than one locus is a real cluster even if its widest member's own offset is 0/0 (that
    # member's interval spans the whole vc_region, so it looks solo by offset alone). A vc_region
    # with only one locus is real too if that locus's own offset is non-zero: the region was
    # genuinely extended beyond the locus's original_region, even though nothing else was merged in.
    # Only a singleton vc_region with 0/0 offset is just that locus's own region and nothing more.
    vc_region_to_loci = {
        region: loci for region, loci in vc_region_to_loci.items()
        if len(loci) > 1 or region in vc_regions_with_nonzero_offset
    }

    print(f"Parsed TSV: {len(vc_region_to_loci):,d} variation clusters, "
          f"{len(filter_reason_by_locus_id):,d} loci flagged DEPTH or EXTENSION")

    output_bed_file = open(args.output_bed_path, "wt")
    counters = collections.Counter()

    # Write one row per repeat in the catalog
    for record_i, record in enumerate(get_variant_catalog_iterator(
            args.input_repeat_catalog, show_progress_bar=args.show_progress_bar)):
        # split_adjacent_repeats matches how the repeat catalog's own TRGT file is built (step 41 of
        # run_all_steps_to_generate_the_catalog.py), so a compound record becomes one row per adjacent
        # repeat here too rather than a single merged row. The converter's reference fasta argument is
        # only used to normalize chromosome names, which this catalog already writes with a chr prefix.
        for output_row in convert_expansion_hunter_record_to_trgt_rows(
                record_i, record, split_adjacent_repeats=True):
            # The converter already writes ID, MOTIFS and the "(CAG)n" structure expression, which is
            # what the repeat catalog's own TRGT file carries. Leave all three alone so the repeat
            # rows here stay comparable to that file, and append the filter reason if there is one.
            filter_reason = filter_reason_by_locus_id.get(parse_info_field(output_row[3])["ID"])
            if filter_reason:
                output_row[3] += f";VC_FILTER={filter_reason}"
                counters[f"repeat rows flagged VC_FILTER={filter_reason}"] += 1
            counters["repeat rows"] += 1
            if f"{output_row[0]}:{output_row[1]}-{output_row[2]}" in vc_region_to_loci:
                counters["repeat rows whose span equals a variation cluster's"] += 1
            output_bed_file.write("\t".join(map(str, output_row)) + "\n")

    # Write one additional row per variation cluster
    for vc_region, loci in vc_region_to_loci.items():
        chrom, start, end = parse_interval(vc_region)

        # Get unique motifs while preserving order
        seen_motifs = set()
        unique_motifs = []
        for _, motifs in loci:
            for motif in motifs.split(","):
                if motif and motif not in seen_motifs:
                    seen_motifs.add(motif)
                    unique_motifs.append(motif)

        # A variation cluster gets an ID derived from its own coordinates (without the chr prefix)
        # rather than from the IDs of the repeats it contains. Previously a VC containing a single
        # repeat was given that repeat's ID, which made the two indistinguishable in TRGT output
        # (https://github.com/PacificBiosciences/trgt-lps/issues/5). The list of repeats that the VC
        # contains moves into STRUC, which TRGT copies through to the VCF unchanged.
        output_info = (f"ID=VC:{chrom.replace('chr', '')}:{start}-{end};"
                       f"MOTIFS={','.join(unique_motifs)};"
                       f"STRUC=<VC:{','.join(locus_id for locus_id, _ in loci)}>")
        output_bed_file.write(f"{chrom}\t{start}\t{end}\t{output_info}\n")
        counters["variation cluster rows"] += 1

    output_bed_file.close()

    # Sort and compress output
    run(f"bedtools sort -i {args.output_bed_path} | bgzip > {args.output_bed_path}.sorted")
    run(f"mv {args.output_bed_path}.sorted {args.output_bed_path}.gz")
    os.remove(args.output_bed_path)

    for label, count in sorted(counters.items()):
        print(f"  {count:>12,d}  {label}")
    print(f"Wrote {counters['repeat rows'] + counters['variation cluster rows']:,d} rows "
          f"({counters['repeat rows']:,d} repeats + {counters['variation cluster rows']:,d} clusters) "
          f"to {args.output_bed_path}.gz")


if __name__ == "__main__":
    main()
