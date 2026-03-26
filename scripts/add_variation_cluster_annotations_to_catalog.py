"""Add variation cluster annotations to catalog.

This script reads a variation clusters TSV file and adds the following annotations to the catalog:
- VariationCluster: The genomic interval of the variation cluster (if offsets are non-zero)
- VariationClusterId: Comma-separated locus IDs of all loci in the same variation cluster (if offsets are non-zero)
- VariationClusterMotifs: Comma-separated unique motifs from all loci in the same variation cluster (if offsets are non-zero)
- VariationClusterSizeDiff: Sum of start and end offsets (if offsets are non-zero)
- VariationClusterFilterReason: "DEPTH" or "EXTENSION" if the locus was filtered from variation clusters
"""

import argparse
import collections
import gzip
import ijson
import os
import simplejson as json
import tqdm

from str_analysis.utils.misc_utils import parse_interval


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
    parser.add_argument("--known-pathogenic-loci-json-path", required=True, help="Path of ExpansionHunter catalog "
                        "containing known pathogenic loci. This is used to retrieve the original locus boundaries for "
                        "these loci since their IDs don't contain these coordinates the way that IDs of other loci do.")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    parser.add_argument("--output-catalog-json-path", help="Path of the output catalog JSON file that includes variation cluster annotations")
    parser.add_argument("--generate-plot", action="store_true", help="Generate a plot of the size differences between variation clusters and simple repeats")
    parser.add_argument("variation_clusters_tsv_path", help="Path of the variation clusters TSV file")
    parser.add_argument("catalog_json_path", help="Path of the JSON catalog to annotate")
    args = parser.parse_args()

    for path in args.variation_clusters_tsv_path, args.catalog_json_path, args.known_pathogenic_loci_json_path:
        if not os.path.isfile(path):
            parser.error(f"File not found: {path}")

    if not args.output_catalog_json_path:
        args.output_catalog_json_path = args.catalog_json_path.replace(".json", ".with_variation_clusters.json")

    print(f"Parsing {args.known_pathogenic_loci_json_path}")
    fopen = gzip.open if args.known_pathogenic_loci_json_path.endswith("gz") else open
    with fopen(args.known_pathogenic_loci_json_path, "rt") as f:
        known_pathogenic_loci = json.load(f)
        known_pathogenic_reference_regions_lookup = {}
        for locus in known_pathogenic_loci:
            if isinstance(locus["ReferenceRegion"], list):
                assert isinstance(locus["VariantId"], list)
                assert len(locus["ReferenceRegion"]) == len(locus["VariantId"])
                for variant_id, reference_region in zip(locus["VariantId"], locus["ReferenceRegion"]):
                    known_pathogenic_reference_regions_lookup[variant_id] = reference_region
            else:
                known_pathogenic_reference_regions_lookup[locus["LocusId"]] = locus["ReferenceRegion"]

    # Data structures to store variation cluster info for each locus
    locus_id_to_variation_cluster_interval = {}
    locus_id_to_variation_cluster_size_diff = {}
    locus_id_to_filter_reason = {}
    vc_region_to_locus_ids = collections.defaultdict(list)
    vc_region_to_motifs = collections.defaultdict(list)

    # Counters for statistics
    size_diff_histogram = collections.Counter()
    all_tsv_locus_ids = set()
    input_locus_counter = 0
    loci_with_variation_cluster = 0
    loci_with_filter_reason = 0
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
            else:
                # Parse offsets as floats, then convert to int for size diff
                try:
                    start_offset = float(vc_start_offset) if vc_start_offset else 0.0
                    end_offset = float(vc_end_offset) if vc_end_offset else 0.0
                except ValueError:
                    print(f"WARNING: Could not parse offsets for locus {locus_id}: start='{vc_start_offset}', end='{vc_end_offset}'")
                    continue

                # Only add variation cluster annotation if offsets are non-zero
                if start_offset != 0 or end_offset != 0:
                    size_diff = int(abs(start_offset) + abs(end_offset))
                    locus_id_to_variation_cluster_interval[locus_id] = vc_region
                    locus_id_to_variation_cluster_size_diff[locus_id] = size_diff
                    vc_region_to_locus_ids[vc_region].append(locus_id)
                    vc_region_to_motifs[vc_region].append(motifs)
                    size_diff_histogram[size_diff] += 1
                    loci_with_variation_cluster += 1
                else:
                    loci_with_zero_offset += 1

    # Build locus_id -> variation cluster ID (comma-joined locus IDs sharing the same vc_region)
    locus_id_to_variation_cluster_id = {}
    locus_id_to_variation_cluster_motifs = {}
    for vc_region, locus_ids in vc_region_to_locus_ids.items():
        variation_cluster_id = ",".join(locus_ids)
        # Get unique motifs while preserving order
        seen_motifs = set()
        unique_motifs = []
        for motifs in vc_region_to_motifs[vc_region]:
            for m in motifs.split(","):
                if m and m not in seen_motifs:
                    seen_motifs.add(m)
                    unique_motifs.append(m)
        variation_cluster_motifs = ",".join(unique_motifs)
        for locus_id in locus_ids:
            locus_id_to_variation_cluster_id[locus_id] = variation_cluster_id
            locus_id_to_variation_cluster_motifs[locus_id] = variation_cluster_motifs

    if args.verbose:
        print(f"Parsed {input_locus_counter:,d} loci from {args.variation_clusters_tsv_path}")
        print(f"  - {loci_with_variation_cluster:,d} ({loci_with_variation_cluster/input_locus_counter:.1%}) have non-zero offsets (will get VariationCluster annotation)")
        print(f"  - {loci_with_filter_reason:,d} ({loci_with_filter_reason/input_locus_counter:.1%}) were filtered (DEPTH or EXTENSION)")
        print(f"  - {loci_with_zero_offset:,d} ({loci_with_zero_offset/input_locus_counter:.1%}) have zero offsets (no annotation)")

    print(f"Annotating {args.catalog_json_path} with variation cluster annotations")
    fopen = gzip.open if args.catalog_json_path.endswith("gz") else open
    with fopen(args.catalog_json_path, "rt") as f:
        f2open = gzip.open if args.output_catalog_json_path.endswith("gz") else open
        with f2open(args.output_catalog_json_path, "wt") as f2:
            input_locus_counter = 0
            locus_with_vc_annotation_counter = 0
            locus_with_filter_annotation_counter = 0
            locus_without_annotation_counter = 0
            catalog_locus_ids = set()

            iterator = ijson.items(f, "item")
            if args.show_progress_bar:
                iterator = tqdm.tqdm(iterator, unit=" records", unit_scale=True)

            f2.write("[")
            for i, record in enumerate(iterator):
                locus_id = record["LocusId"]
                input_locus_counter += 1
                catalog_locus_ids.add(locus_id)

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

                if i > 0:
                    f2.write(", ")
                f2.write(json.dumps(record, f2, use_decimal=True, indent=4))

            f2.write("]")

    print(f"Annotated {input_locus_counter:,d} loci from {args.catalog_json_path}")
    print(f"  - {locus_with_vc_annotation_counter:,d} ({locus_with_vc_annotation_counter/input_locus_counter:.1%}) got VariationCluster annotation")
    print(f"  - {locus_with_filter_annotation_counter:,d} ({locus_with_filter_annotation_counter/input_locus_counter:.1%}) got VariationClusterFilterReason annotation")
    print(f"  - {locus_without_annotation_counter:,d} ({locus_without_annotation_counter/input_locus_counter:.1%}) got no variation cluster annotation")
    print(f"Wrote output to {args.output_catalog_json_path}")

    # Validate that all locus IDs in the variation clusters TSV have an exact match in the catalog
    vc_locus_ids_not_in_catalog = all_tsv_locus_ids - catalog_locus_ids
    if vc_locus_ids_not_in_catalog:
        raise ValueError(
            f"{len(vc_locus_ids_not_in_catalog):,d} locus ID(s) in {args.variation_clusters_tsv_path} were not found "
            f"in {args.catalog_json_path}. Examples: {sorted(vc_locus_ids_not_in_catalog)[:10]}")

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
