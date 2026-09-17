"""Add variation cluster annotations to catalog.

This script reads a variation clusters TSV file and adds the following annotations to the catalog.

Variation cluster membership is taken directly from the TSV: rows that share the same ``vc_region``
are members of the same variation cluster. Each VC's identity is the comma-joined list of its
member TRIDs. A TRID may be listed as a member of more than one VC; when that happens, the widest
VC wins. (Today's TSV lists each TRID under exactly one ``vc_region``, so each TRID has a
length-1 candidate list; the widest-wins logic is forward-compatible with future TSV formats.)

Annotations added to each catalog record:
- VariationCluster: The genomic interval of the chosen VC
- VariationClusterId: Comma-separated locus IDs of all members of the chosen VC
- VariationClusterMotifs: Comma-separated unique motifs across those members
- VariationClusterSizeDiff: (chosen vc_region size) - (locus original_region size)
- VariationClusterFilterReason: "DEPTH" or "EXTENSION" if the locus was filtered from variation clusters

This script owns those five fields: each run clears them first and then writes only the ones that
apply, so running it again over an already-annotated catalog replaces the annotations rather than
layering onto them.
"""

import argparse
import collections
import gzip
import ijson
import os
import simplejson as json
import tqdm

from str_analysis.utils.misc_utils import parse_interval
from catalog_annotation_utils import clear_previous_annotations, print_annotation_replacement_summary


# Every field this script writes. A locus that no longer belongs to a variation cluster, or that
# stops being filtered, must lose the values a previous run gave it, so these are cleared per record
# before the new ones are written.
VARIATION_CLUSTER_FIELDS = (
    "VariationCluster",
    "VariationClusterId",
    "VariationClusterMotifs",
    "VariationClusterSizeDiff",
    "VariationClusterFilterReason",
)


def parse_info_field(info_field):
    """Parse a TRGT catalog info field into a python dictionary"""
    result = {}
    for key_value in info_field.split(";"):
        key_value = key_value.split("=")
        if len(key_value) != 2:
            raise ValueError(f"Invalid key-value pair '{key_value}' in info field: {info_field}")
        key, value = key_value
        result[key] = value
    return result


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    parser.add_argument("--output-catalog-json-path", help="Path of the output catalog JSON file that includes variation cluster annotations")
    parser.add_argument("--generate-plot", action="store_true", help="Generate a plot of the size differences between variation clusters and simple repeats")
    parser.add_argument("variation_clusters_tsv_path", help="Path of the variation clusters TSV file")
    parser.add_argument("catalog_json_path", help="Path of the JSON catalog to annotate")
    args = parser.parse_args()

    for path in args.variation_clusters_tsv_path, args.catalog_json_path:
        if not os.path.isfile(path):
            parser.error(f"File not found: {path}")

    if not args.output_catalog_json_path:
        args.output_catalog_json_path = args.catalog_json_path.replace(".json", ".with_variation_clusters.json")

    # Per-locus data from the TSV.
    locus_id_to_original_length = {}     # locus_id -> end - start of original_region (for size_diff)
    locus_id_to_filter_reason = {}       # locus_id -> "DEPTH" or "EXTENSION"
    # locus_id -> list of vc_regions whose TSV-derived member set contains this locus_id.
    # Today each TRID has at most one entry (each TSV row pairs one TRID with one vc_region);
    # this stays a list to support future TSV formats where a TRID may appear in multiple rows.
    locus_id_to_candidate_vc_regions = collections.defaultdict(list)
    # vc_region (string "chr:start-end") -> ordered list of (locus_id, motifs) for its members.
    # Order is TSV-encounter order so the comma-joined VariationClusterId is deterministic.
    vc_region_to_members = collections.defaultdict(list)
    vc_region_to_length = {}             # vc_region -> end - start, for widest-VC pick
    vc_regions_with_nonzero_offset = set()

    all_tsv_locus_ids = set()
    input_locus_counter = 0
    loci_with_filter_reason = 0
    loci_with_nonzero_offsets = 0
    loci_with_zero_offset = 0

    if args.verbose:
        print(f"Parsing {args.variation_clusters_tsv_path}")

    fopen = gzip.open if args.variation_clusters_tsv_path.endswith("gz") else open
    with fopen(args.variation_clusters_tsv_path, "rt") as f:
        if args.show_progress_bar:
            f = tqdm.tqdm(f, unit=" records", unit_scale=True)

        header = None
        for line in f:
            fields = line.strip("\n").split("\t")

            # Parse header
            if header is None:
                header = fields
                continue

            input_locus_counter += 1

            # Parse TSV columns
            region_info = fields[0]
            original_region = fields[1]
            vc_start_offset = fields[2]
            vc_end_offset = fields[3]
            vc_region = fields[4] if len(fields) > 4 else ""

            # Extract locus ID from region_info
            info_dict = parse_info_field(region_info)
            locus_id = info_dict["ID"]
            motifs = info_dict.get("MOTIFS", "")
            all_tsv_locus_ids.add(locus_id)

            # Check if this locus was filtered
            if vc_end_offset in ("DEPTH", "EXTENSION"):
                locus_id_to_filter_reason[locus_id] = vc_end_offset
                loci_with_filter_reason += 1
                continue

            # Parse offsets as floats. A locus's own offset can't be used to decide cluster
            # membership on its own: the widest member of a real multi-locus cluster spans the
            # whole vc_region and so also has offset 0/0, same as a genuinely solo locus. But it
            # does matter for a singleton vc_region (see the filter after the loop): a lone locus
            # whose vc_region was genuinely extended past its own original_region is still real.
            try:
                start_offset = float(vc_start_offset) if vc_start_offset else 0.0
                end_offset = float(vc_end_offset) if vc_end_offset else 0.0
            except ValueError:
                print(f"WARNING: Could not parse offsets for locus {locus_id}: start='{vc_start_offset}', end='{vc_end_offset}'")
                continue

            if start_offset == 0 and end_offset == 0:
                loci_with_zero_offset += 1
            else:
                loci_with_nonzero_offsets += 1
                vc_regions_with_nonzero_offset.add(vc_region)

            # Record (locus_id, motifs) as a member of vc_region, and remember vc_region as a
            # candidate VC for this locus. Both directions are needed: vc_region->members to build
            # VariationClusterId/Motifs, locus_id->candidates to pick the widest VC per locus.
            vc_region_to_members[vc_region].append((locus_id, motifs))
            locus_id_to_candidate_vc_regions[locus_id].append(vc_region)

            if vc_region not in vc_region_to_length:
                _, vc_start, vc_end = parse_interval(vc_region)
                vc_region_to_length[vc_region] = vc_end - vc_start

            orig_chrom, orig_start, orig_end = parse_interval(original_region)
            locus_id_to_original_length[locus_id] = orig_end - orig_start

    # Keep a vc_region unless it has exactly one member whose own offset is 0/0. A vc_region shared
    # by more than one locus is a real cluster even if its widest member's own offset is 0/0 (that
    # member's interval spans the whole vc_region, so it looks solo by offset alone). A vc_region
    # with only one locus is real too if that locus's own offset is non-zero: the region was
    # genuinely extended beyond the locus's original_region, even though nothing else was merged
    # in. Only a singleton vc_region with 0/0 offset is just that locus's own region and nothing
    # more.
    vc_region_to_members = {
        region: members for region, members in vc_region_to_members.items()
        if len(members) > 1 or region in vc_regions_with_nonzero_offset
    }
    locus_id_to_candidate_vc_regions = {
        locus_id: [vc_region for vc_region in candidates if vc_region in vc_region_to_members]
        for locus_id, candidates in locus_id_to_candidate_vc_regions.items()
    }
    locus_id_to_candidate_vc_regions = {
        locus_id: candidates for locus_id, candidates in locus_id_to_candidate_vc_regions.items() if candidates
    }

    # For each locus, pick the widest vc_region from its candidate list. Today every list has
    # length 1, so "widest" is trivially the only entry; if future TSVs list a TRID under
    # multiple vc_regions, this picks the widest.
    locus_id_to_variation_cluster_interval = {}
    locus_id_to_variation_cluster_size_diff = {}
    size_diff_histogram = collections.Counter()

    for locus_id, candidate_vc_regions in locus_id_to_candidate_vc_regions.items():
        # Tie-break by the vc_region string so the choice is deterministic when multiple
        # candidates are equally wide.
        widest_vc_region = max(candidate_vc_regions, key=lambda vc: (vc_region_to_length[vc], vc))
        size_diff = vc_region_to_length[widest_vc_region] - locus_id_to_original_length[locus_id]

        locus_id_to_variation_cluster_interval[locus_id] = widest_vc_region
        locus_id_to_variation_cluster_size_diff[locus_id] = size_diff
        size_diff_histogram[size_diff] += 1

    loci_with_variation_cluster = len(locus_id_to_variation_cluster_interval)

    # Build locus_id -> VariationClusterId (comma-joined member TRIDs of the chosen VC) and
    # VariationClusterMotifs (comma-joined unique motifs across those members). The cluster's
    # full membership comes from the TSV grouping in vc_region_to_members; the (locus_id, motifs)
    # tuples are in TSV-encounter order so the comma-joined output is deterministic.
    locus_id_to_variation_cluster_id = {}
    locus_id_to_variation_cluster_motifs = {}
    for vc_region, members in vc_region_to_members.items():
        # If no locus chose this vc_region as its widest, no record will reference these strings;
        # skipping avoids spending memory on unused entries.
        if not any(locus_id_to_variation_cluster_interval.get(locus_id) == vc_region for locus_id, _ in members):
            continue
        variation_cluster_id = ",".join(locus_id for locus_id, _ in members)
        seen_motifs = set()
        unique_motifs = []
        for _, motifs in members:
            for m in motifs.split(","):
                if m and m not in seen_motifs:
                    seen_motifs.add(m)
                    unique_motifs.append(m)
        variation_cluster_motifs = ",".join(unique_motifs)
        for locus_id, _ in members:
            if locus_id_to_variation_cluster_interval.get(locus_id) == vc_region:
                locus_id_to_variation_cluster_id[locus_id] = variation_cluster_id
                locus_id_to_variation_cluster_motifs[locus_id] = variation_cluster_motifs

    if args.verbose:
        print(f"Parsed {input_locus_counter:,d} loci from {args.variation_clusters_tsv_path}")
        print(f"  - {loci_with_nonzero_offsets:,d} ({loci_with_nonzero_offsets/input_locus_counter:.1%}) have non-zero offsets in the TSV")
        print(f"  - {loci_with_filter_reason:,d} ({loci_with_filter_reason/input_locus_counter:.1%}) were filtered (DEPTH or EXTENSION)")
        print(f"  - {loci_with_zero_offset:,d} ({loci_with_zero_offset/input_locus_counter:.1%}) have zero offsets in the TSV")
        print(f"  - {loci_with_variation_cluster:,d} ({loci_with_variation_cluster/input_locus_counter:.1%}) "
              f"will get VariationCluster annotation (widest containing VC)")

    print(f"Annotating {args.catalog_json_path} with variation cluster annotations")
    fopen = gzip.open if args.catalog_json_path.endswith("gz") else open
    with fopen(args.catalog_json_path, "rt") as f:
        f2open = gzip.open if args.output_catalog_json_path.endswith("gz") else open
        with f2open(args.output_catalog_json_path, "wt") as f2:
            input_locus_counter = 0
            locus_with_vc_annotation_counter = 0
            locus_with_filter_annotation_counter = 0
            locus_without_annotation_counter = 0
            locus_with_previous_annotation_counter = 0
            locus_that_lost_annotation_counter = 0
            catalog_locus_ids = set()

            iterator = ijson.items(f, "item")
            if args.show_progress_bar:
                iterator = tqdm.tqdm(iterator, unit=" records", unit_scale=True)

            f2.write("[")
            for i, record in enumerate(iterator):
                locus_id = record["LocusId"]
                input_locus_counter += 1
                catalog_locus_ids.add(locus_id)

                had_previous_annotation = clear_previous_annotations(record, VARIATION_CLUSTER_FIELDS)
                if had_previous_annotation:
                    locus_with_previous_annotation_counter += 1

                if locus_id in locus_id_to_variation_cluster_interval:
                    record["VariationCluster"] = locus_id_to_variation_cluster_interval[locus_id]
                    record["VariationClusterId"] = locus_id_to_variation_cluster_id[locus_id]
                    record["VariationClusterMotifs"] = locus_id_to_variation_cluster_motifs[locus_id]
                    record["VariationClusterSizeDiff"] = locus_id_to_variation_cluster_size_diff[locus_id]
                    locus_with_vc_annotation_counter += 1
                elif locus_id in locus_id_to_filter_reason:
                    record["VariationClusterFilterReason"] = locus_id_to_filter_reason[locus_id]
                    locus_with_filter_annotation_counter += 1
                else:
                    locus_without_annotation_counter += 1
                    if had_previous_annotation:
                        locus_that_lost_annotation_counter += 1

                if i > 0:
                    f2.write(", ")
                f2.write(json.dumps(record, use_decimal=True, indent=4))

            f2.write("]")

    print(f"Annotated {input_locus_counter:,d} loci from {args.catalog_json_path}")
    print(f"  - {locus_with_vc_annotation_counter:,d} ({locus_with_vc_annotation_counter/input_locus_counter:.1%}) got VariationCluster annotation")
    print(f"  - {locus_with_filter_annotation_counter:,d} ({locus_with_filter_annotation_counter/input_locus_counter:.1%}) got VariationClusterFilterReason annotation")
    print(f"  - {locus_without_annotation_counter:,d} ({locus_without_annotation_counter/input_locus_counter:.1%}) got no variation cluster annotation")
    print_annotation_replacement_summary("variation cluster",
                                         "no longer belong to a variation cluster and are not filtered",
                                         locus_with_previous_annotation_counter,
                                         locus_that_lost_annotation_counter)
    print(f"Wrote output to {args.output_catalog_json_path}")

    # Validate that all locus IDs in the variation clusters TSV have an exact match in the catalog.
    # A small number of missing IDs is tolerated as an upstream-data quirk; only raise above the threshold.
    # Measured on 2026-08-28: 854 of the 5,543,666 locus IDs in genome_clusters-v2.tsv.gz are absent
    # from the v2 catalog (5,599,658 loci in
    # TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz). The threshold
    # below therefore sits only about 5% above the count it was chosen against, so re-measure it
    # whenever either input is regenerated: a modest shift upstream crosses it and turns this warning
    # back into a hard failure.
    MAX_TOLERATED_MISSING = 900
    vc_locus_ids_not_in_catalog = all_tsv_locus_ids - catalog_locus_ids
    if len(vc_locus_ids_not_in_catalog) >= MAX_TOLERATED_MISSING:
        raise ValueError(
            f"{len(vc_locus_ids_not_in_catalog):,d} locus ID(s) in {args.variation_clusters_tsv_path} were not found "
            f"in {args.catalog_json_path}. Examples: {sorted(vc_locus_ids_not_in_catalog)[:10]}")
    elif vc_locus_ids_not_in_catalog:
        print(f"WARNING: {len(vc_locus_ids_not_in_catalog):,d} locus ID(s) in {args.variation_clusters_tsv_path} were not found "
              f"in {args.catalog_json_path} (under threshold of {MAX_TOLERATED_MISSING:,d}). "
              f"Examples: {sorted(vc_locus_ids_not_in_catalog)[:10]}")

    if args.generate_plot and size_diff_histogram:
        print(f"Generating VC size diff histograms")
        import seaborn as sns
        import matplotlib.pyplot as plt
        plt.figure(figsize=(12, 6))
        sns.barplot(x=list(size_diff_histogram.keys()), y=list(size_diff_histogram.values()))
        plt.xlabel("Size difference")
        plt.ylabel("Count")
        plt.title("Size difference between variation clusters and original loci")
        output_prefix = args.output_catalog_json_path.replace(".json", "").replace(".gz", "") + ".size_diff_histogram"
        plt.savefig(f"{output_prefix}.png")
        plt.yscale("log")
        plt.savefig(f"{output_prefix}.log.png")
        print(f"Wrote VC size diff histograms to {output_prefix}.png and {output_prefix}.log.png")


if __name__ == "__main__":
    main()
