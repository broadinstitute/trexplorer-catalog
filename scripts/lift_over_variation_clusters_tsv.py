"""Lift the variation clusters table over to another reference genome, for a catalog lifted over from the original one.

The variation clusters TSV (genome_clusters-v2.tsv.gz, shared by Egor Dolzhenko) has one row per locus of the hg38
catalog, keyed by its LocusId, with the variation cluster (vc_region) the locus belongs to, or with the reason (DEPTH or
EXTENSION) that the clustering couldn't process the locus. This writes the same table for the CHM13 catalog, so that the
hg38 pipeline's variation cluster scripts (add_variation_cluster_annotations_to_catalog.py and
generate_TRGT_catalog_with_variation_clusters_and_all_repeats.py) can use it unchanged:

- each locus is replaced by the CHM13 catalog locus whose Hg38LocusId is the locus's LocusId. When two CHM13 loci have
  the same Hg38LocusId, the one with the same boundaries (Hg38LocusMatchType SameBoundaries) is used. Loci with no CHM13
  counterpart are left out.
- each variation cluster's vc_region is lifted with lift_over_catalog.lift_interval, which maps the bases on either side
  of it through the chain file, and then widened, if needed, so that it contains the CHM13 loci of all of its members.
  Its members are written as loci that aren't in a variation cluster (with their own interval as vc_region, and
  offsets of 0, which is how the TSV describes such loci) when the vc_region doesn't lift over, or when it lifts to a
  region more than MAX_VC_REGION_SIZE_INCREASE larger than itself (like a 246bp region that lifts to 82kb, since the
  chain has a CHM13-specific insertion between its flanks). So is a member whose CHM13 locus doesn't overlap the lifted
  vc_region, like one that lifted through another chain to a paralogous copy megabases away, since widening the
  cluster to reach it would turn it into a cluster of that size.
- the offsets (how far the vc_region extends beyond the locus on each side) are recomputed from the CHM13 coordinates.

CHM13 loci with no hg38 counterpart, like most of those in CHM13's satellite arrays, aren't in the table. Like the hg38
catalog's loci from sources added after the table was made, they get no variation cluster annotations, and the TRGT
catalog with variation clusters writes each of them as its own row only.
"""

import argparse
import collections
import gzip

import ijson

from lift_over_catalog import lift_interval, parse_chain_file
from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure
from str_analysis.utils.misc_utils import parse_interval

FILTER_REASONS = ("DEPTH", "EXTENSION")

# How far lift_interval may move an anchor that falls in a chain gap, as in lift_over_catalog.py
MAX_ANCHOR_SEARCH_DISTANCE = 20

# How much larger than the source vc_region its lifted region can be. Of the 827,968 hg38 variation clusters that lift to
# CHM13, 69 lift to a region more than this much larger (up to 82kb larger), while 99.99% lift to a region at most about
# 9.7kb larger.
MAX_VC_REGION_SIZE_INCREASE = 10_000

TargetLocus = collections.namedtuple("TargetLocus", ["locus_id", "chrom", "start_0based", "end_1based", "motifs"])
TableRow = collections.namedtuple("TableRow", ["locus_id", "original_region", "vc_start_offset", "vc_end_offset",
                                               "vc_region"])


def read_target_loci_by_source_locus_id(catalog_records, source_locus_id_field="Hg38LocusId"):
    """Returns the catalog's loci keyed by the LocusId of their counterpart in the other reference genome.

    Args:
        catalog_records (iterable): records of the lifted catalog
        source_locus_id_field (str): field that holds the LocusId of each locus's counterpart

    Returns:
        dict: source LocusId => TargetLocus. When several loci have the same source LocusId, the first one whose
            Hg38LocusMatchType is SameBoundaries, or else the first one, is kept.
    """
    target_loci = {}
    for record in catalog_records:
        source_locus_id = record.get(source_locus_id_field)
        if not source_locus_id:
            continue
        if source_locus_id in target_loci and (
                target_loci[source_locus_id][1] or record.get("Hg38LocusMatchType") != "SameBoundaries"):
            continue

        reference_regions = record["ReferenceRegion"]
        if not isinstance(reference_regions, list):
            reference_regions = [reference_regions]
        intervals = [parse_interval(region) for region in reference_regions]
        target_locus = TargetLocus(
            record["LocusId"], intervals[0][0], min(start for _, start, _ in intervals),
            max(end for _, _, end in intervals), parse_motifs_from_locus_structure(record["LocusStructure"]))
        target_loci[source_locus_id] = (target_locus, record.get("Hg38LocusMatchType") == "SameBoundaries")

    return {source_locus_id: target_locus for source_locus_id, (target_locus, _) in target_loci.items()}


def parse_table_row(line):
    """Parses one row of the variation clusters TSV.

    Args:
        line (str): the row, eg. "ID=1-19175-19184-TCC;MOTIFS=TCC;STRUC=(TCC)n\tchr1:19175-19184\t4.0\t0\tchr1:19171-19184"

    Returns:
        TableRow: the row's fields, with the LocusId taken from the ID key of its first column
    """
    fields = line.rstrip("\n").split("\t")
    fields += [""] * (5 - len(fields))
    info = dict(key_value.split("=", 1) for key_value in fields[0].split(";"))
    return TableRow(info["ID"], *fields[1:5])


def format_table_row(target_locus, vc_start_offset, vc_end_offset, vc_region):
    """Returns a row of the variation clusters TSV for a locus of the lifted catalog, without a trailing newline."""
    region_info = (f"ID={target_locus.locus_id};MOTIFS={','.join(target_locus.motifs)};"
                   f"STRUC={''.join(f'({motif})n' for motif in target_locus.motifs)}")
    original_region = f"{target_locus.chrom}:{target_locus.start_0based}-{target_locus.end_1based}"
    return "\t".join([region_info, original_region, vc_start_offset, vc_end_offset, vc_region])


def format_row_without_variation_cluster(target_locus):
    """Returns the row of a locus that isn't in a variation cluster: its vc_region is its own interval."""
    own_region = f"{target_locus.chrom}:{target_locus.start_0based}-{target_locus.end_1based}"
    return format_table_row(target_locus, "0.0", "0", own_region)


def lift_variation_cluster(vc_region, member_target_loci, chain_index, counters):
    """Lifts a variation cluster, and returns the table rows of its members.

    Args:
        vc_region (str): the variation cluster's interval in the source reference, eg. "chr1:19171-19184"
        member_target_loci (list): TargetLocus of each of the cluster's members that has a counterpart in the lifted
            catalog
        chain_index (dict): output of lift_over_catalog.parse_chain_file
        counters (collections.Counter): counts of what happened, updated in place

    Returns:
        list: table rows (str) of the members
    """
    chrom, start_0based, end_1based = parse_interval(vc_region)
    lifted_interval, failure_reason = lift_interval(
        chain_index, chrom, start_0based, end_1based, max_anchor_search_distance=MAX_ANCHOR_SEARCH_DISTANCE)
    if lifted_interval is None:
        counters[f"variation clusters whose vc_region didn't lift over: {failure_reason}"] += 1
        return [format_row_without_variation_cluster(locus) for locus in member_target_loci]

    target_chrom, target_start_0based, target_end_1based, _ = lifted_interval
    if (target_end_1based - target_start_0based) - (end_1based - start_0based) > MAX_VC_REGION_SIZE_INCREASE:
        counters[f"variation clusters whose vc_region lifted to a region more than {MAX_VC_REGION_SIZE_INCREASE:,d}bp "
                 f"larger"] += 1
        return [format_row_without_variation_cluster(locus) for locus in member_target_loci]

    def overlaps_lifted_vc_region(locus):
        return locus.chrom == target_chrom and locus.start_0based < target_end_1based and \
            locus.end_1based > target_start_0based

    members_in_lifted_vc_region = [locus for locus in member_target_loci if overlaps_lifted_vc_region(locus)]
    if len(members_in_lifted_vc_region) < len(member_target_loci):
        counters["members whose CHM13 locus doesn't overlap their lifted variation cluster"] += (
            len(member_target_loci) - len(members_in_lifted_vc_region))
    if not members_in_lifted_vc_region:
        counters["variation clusters with no member that overlaps their lifted vc_region"] += 1
        return [format_row_without_variation_cluster(locus) for locus in member_target_loci]

    vc_start_0based = min([target_start_0based] + [locus.start_0based for locus in members_in_lifted_vc_region])
    vc_end_1based = max([target_end_1based] + [locus.end_1based for locus in members_in_lifted_vc_region])
    if (vc_start_0based, vc_end_1based) != (target_start_0based, target_end_1based):
        counters["variation clusters widened to contain the CHM13 loci of their members"] += 1
    counters["variation clusters lifted over"] += 1

    lifted_vc_region = f"{target_chrom}:{vc_start_0based}-{vc_end_1based}"
    rows = []
    for locus in member_target_loci:
        if not overlaps_lifted_vc_region(locus):
            rows.append(format_row_without_variation_cluster(locus))
            continue
        rows.append(format_table_row(
            locus, f"{float(locus.start_0based - vc_start_0based)}", f"{vc_end_1based - locus.end_1based}",
            lifted_vc_region))
    return rows


def is_in_a_variation_cluster(rows):
    """Whether the table rows that share a vc_region describe a variation cluster, using the same rule as
    add_variation_cluster_annotations_to_catalog.py: more than one locus, or one locus with a non-zero offset."""
    return len(rows) > 1 or any(float(row.vc_start_offset or 0) != 0 or float(row.vc_end_offset or 0) != 0
                                for row in rows)


def lift_variation_clusters_table(table_rows, target_loci_by_source_locus_id, chain_index, counters):
    """Returns the rows of the lifted variation clusters table.

    Args:
        table_rows (iterable): TableRow of each row of the source table, in order
        target_loci_by_source_locus_id (dict): output of read_target_loci_by_source_locus_id
        chain_index (dict): output of lift_over_catalog.parse_chain_file
        counters (collections.Counter): counts of what happened, updated in place

    Returns:
        list: rows (str) of the lifted table, in the order of the source table's first row of each variation cluster
    """
    # The output follows the source table's order: a filtered locus's row where it was, and a vc_region's rows where
    # its first member was. So this lists, in order, ("filtered row", row text) and ("vc_region", vc_region) entries.
    rows_by_vc_region = collections.defaultdict(list)
    output_entries = []
    for row in table_rows:
        counters["source table rows"] += 1
        target_locus = target_loci_by_source_locus_id.get(row.locus_id)
        if target_locus is None:
            counters["source table rows whose locus isn't in the lifted catalog"] += 1
            continue
        if row.vc_end_offset in FILTER_REASONS:
            counters[f"loci marked {row.vc_end_offset}"] += 1
            output_entries.append(("filtered row", format_table_row(target_locus, "", row.vc_end_offset, "")))
            continue
        if row.vc_region not in rows_by_vc_region:
            output_entries.append(("vc_region", row.vc_region))
        rows_by_vc_region[row.vc_region].append((row, target_locus))

    lifted_rows = []
    for entry_type, value in output_entries:
        if entry_type == "filtered row":
            lifted_rows.append(value)
            continue
        members = rows_by_vc_region[value]
        member_target_loci = [target_locus for _, target_locus in members]
        if not is_in_a_variation_cluster([row for row, _ in members]):
            counters["loci not in a variation cluster"] += len(members)
            lifted_rows.extend(format_row_without_variation_cluster(locus) for locus in member_target_loci)
            continue
        counters["variation clusters"] += 1
        lifted_rows.extend(lift_variation_cluster(value, member_target_loci, chain_index, counters))

    return lifted_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--chain-file", required=True, help="UCSC chain file from the table's reference genome to the "
                        "lifted catalog's reference genome, eg. hg38ToHs1.over.chain.gz")
    parser.add_argument("--lifted-catalog", required=True, help="JSON catalog whose loci have an Hg38LocusId field "
                        "with the LocusId of their counterpart in the table's catalog")
    parser.add_argument("-o", "--output-tsv-path", required=True, help="Path of the lifted variation clusters TSV")
    parser.add_argument("variation_clusters_tsv", help="Variation clusters TSV of the table's reference genome")
    args = parser.parse_args()

    with gzip.open(args.lifted_catalog, "rb") as f:
        target_loci_by_source_locus_id = read_target_loci_by_source_locus_id(ijson.items(f, "item"))
    print(f"Read {len(target_loci_by_source_locus_id):,d} loci with an Hg38LocusId from {args.lifted_catalog}")

    chain_index = parse_chain_file(args.chain_file)

    counters = collections.Counter()
    with (gzip.open if args.variation_clusters_tsv.endswith("gz") else open)(args.variation_clusters_tsv, "rt") as f:
        header = f.readline()
        lifted_rows = lift_variation_clusters_table(
            (parse_table_row(line) for line in f), target_loci_by_source_locus_id, chain_index, counters)

    with (gzip.open if args.output_tsv_path.endswith("gz") else open)(args.output_tsv_path, "wt") as f:
        f.write(header)
        for row in lifted_rows:
            f.write(row + "\n")

    for label, count in sorted(counters.items()):
        print(f"  {count:>12,d}  {label}")
    print(f"Wrote {len(lifted_rows):,d} rows to {args.output_tsv_path}")


if __name__ == "__main__":
    main()
