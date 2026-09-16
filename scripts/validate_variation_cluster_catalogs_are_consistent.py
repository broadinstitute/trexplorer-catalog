"""Checks that the variation_clusters_and_all_repeats catalog and the legacy
variation_clusters_and_isolated_TRs catalog agree with each other.

The two files are built by separate scripts from the same variation clusters TSV and repeat
catalog, so nothing keeps them in sync except this check. The all_repeats catalog gives every
repeat its own row, plus one row per variation cluster; the legacy catalog gives a row only to
repeats that are not inside a cluster, plus the same per-cluster rows. So a repeat should have a
standalone row in the legacy catalog if and only if it has an unflagged row in the all_repeats
catalog and is not listed as a member of any variation cluster there. A repeat flagged VC_FILTER
in the all_repeats catalog is expected to be absent from the legacy catalog entirely, since the
legacy generator drops DEPTH/EXTENSION loci rather than keeping a flagged row for them
(https://github.com/PacificBiosciences/trgt-lps/issues/5).
"""

import argparse
import collections
import gzip
import re


def open_maybe_gzipped(path):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def locus_id_of(info_field):
    return re.search(r"ID=([^;]+)", info_field).group(1)


def motifs_of(info_field):
    return re.search(r"MOTIFS=([^;]*)", info_field).group(1)


def vc_filter_of(info_field):
    match = re.search(r";VC_FILTER=([^;]*)", info_field)
    return match.group(1) if match else None


def is_variation_cluster(info_field):
    return locus_id_of(info_field).startswith("VC:")


def members_of_cluster(info_field):
    """Returns the list of member locus ids named in a variation cluster row's STRUC field."""
    return re.search(r"STRUC=<VC:([^>]*)>", info_field).group(1).split(",")


def read_catalog(path):
    """Reads a TRGT catalog BED file.

    Returns:
        tuple: (dict of cluster id -> (chrom, start, end, motifs, struc),
            dict of repeat id -> (chrom, start, end, motifs, vc_filter or None))
    """
    clusters = {}
    repeats = {}
    with open_maybe_gzipped(path) as f:
        for line in f:
            chrom, start, end, info_field = line.rstrip("\n").split("\t")[:4]
            locus_id = locus_id_of(info_field)
            motifs = motifs_of(info_field)
            if is_variation_cluster(info_field):
                struc = re.search(r"STRUC=(.*)$", info_field).group(1)
                clusters[locus_id] = (chrom, start, end, motifs, struc)
            else:
                repeats[locus_id] = (chrom, start, end, motifs, vc_filter_of(info_field))
    return clusters, repeats


def check_clusters_match(clusters_a, clusters_b, label_a, label_b):
    failures = []
    only_in_a = set(clusters_a) - set(clusters_b)
    only_in_b = set(clusters_b) - set(clusters_a)
    if only_in_a:
        failures.append(f"{len(only_in_a):,d} variation clusters are in {label_a} but not "
                         f"{label_b}, e.g. {sorted(only_in_a)[:3]}")
    if only_in_b:
        failures.append(f"{len(only_in_b):,d} variation clusters are in {label_b} but not "
                         f"{label_a}, e.g. {sorted(only_in_b)[:3]}")

    mismatched = [locus_id for locus_id in clusters_a.keys() & clusters_b.keys()
                  if clusters_a[locus_id] != clusters_b[locus_id]]
    if mismatched:
        example = mismatched[0]
        failures.append(f"{len(mismatched):,d} variation clusters differ between {label_a} and "
                         f"{label_b}, e.g. {example}: {clusters_a[example]} vs {clusters_b[example]}")
    return failures


def check_isolated_repeats_match_expectation(clusters_a, repeats_a, repeats_b, label_a, label_b):
    """Returns failures if the legacy catalog's standalone repeat rows don't match what the
    all_repeats catalog implies they should be.
    """
    member_ids = {member_id for cluster in clusters_a.values() for member_id in
                  members_of_cluster(f"STRUC={cluster[4]}")}

    flagged_ids = {locus_id for locus_id, fields in repeats_a.items() if fields[4]}
    expected_isolated_ids = set(repeats_a) - flagged_ids - member_ids
    actual_isolated_ids = set(repeats_b)

    failures = []
    missing = expected_isolated_ids - actual_isolated_ids
    if missing:
        failures.append(f"{len(missing):,d} repeats are unflagged, non-member rows in {label_a} "
                         f"but have no row in {label_b}, e.g. {sorted(missing)[:3]}")

    extra = actual_isolated_ids - expected_isolated_ids
    if extra:
        # Distinguish the two ways a legacy row can be unexpected, since they point at different bugs.
        extra_members = extra & member_ids
        extra_other = extra - member_ids
        if extra_members:
            failures.append(f"{len(extra_members):,d} rows in {label_b} are variation cluster "
                             f"members in {label_a}, so they should not also have a standalone row, "
                             f"e.g. {sorted(extra_members)[:3]}")
        if extra_other:
            failures.append(f"{len(extra_other):,d} rows in {label_b} have no corresponding "
                             f"unflagged row in {label_a} at all, e.g. {sorted(extra_other)[:3]}")

    return failures, len(member_ids), len(flagged_ids)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("variation_clusters_and_all_repeats_bed_path",
                        help="Path of the variation_clusters_and_all_repeats TRGT catalog")
    parser.add_argument("legacy_variation_clusters_and_isolated_trs_bed_path",
                        help="Path of the legacy variation_clusters_and_isolated_TRs TRGT catalog")
    args = parser.parse_args()

    label_a = "the all_repeats catalog"
    label_b = "the legacy isolated_TRs catalog"

    clusters_a, repeats_a = read_catalog(args.variation_clusters_and_all_repeats_bed_path)
    clusters_b, repeats_b = read_catalog(args.legacy_variation_clusters_and_isolated_trs_bed_path)

    isolated_failures, member_count, flagged_count = check_isolated_repeats_match_expectation(
        clusters_a, repeats_a, repeats_b, label_a, label_b)

    print(f"{label_a}: {len(clusters_a):,d} clusters, {len(repeats_a):,d} repeat rows "
          f"({member_count:,d} cluster members, {flagged_count:,d} flagged VC_FILTER)")
    print(f"{label_b}: {len(clusters_b):,d} clusters, {len(repeats_b):,d} repeat rows")

    failures = check_clusters_match(clusters_a, clusters_b, label_a, label_b) + isolated_failures

    for failure in failures:
        print(f"FAILED: {failure}")
    if failures:
        raise SystemExit(1)
    print("All checks passed")


if __name__ == "__main__":
    main()
