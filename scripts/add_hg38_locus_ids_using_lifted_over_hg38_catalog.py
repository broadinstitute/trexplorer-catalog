"""Add the Hg38LocusId of each CHM13 catalog locus's hg38 counterpart, using the hg38 catalog lifted over to CHM13.

The CHM13 pipeline merges the TRExplorer hg38 catalog, lifted over to CHM13, after the known disease-associated loci and
the perfect repeats in CHM13. When a lifted hg38 locus matches one of those earlier loci, the merge keeps the earlier
locus and discards the lifted one together with its Hg38LocusId, which leaves most CHM13 loci (even HTT) without an
Hg38LocusId. This script finds each such locus's hg38 counterpart again among the lifted hg38 loci, using the test that
the merge uses to decide that two loci are the same: str_analysis.merge_loci's check_for_sufficient_overlap_and_motif_match
with the merge's overlap fraction and Jaccard similarity thresholds. When several lifted loci pass, it picks the one with
the same canonical motif (the merge also lets VNTRs with different motifs of the same length match), then the highest
Jaccard similarity, then the smallest start, end and Hg38LocusId, so a lifted locus with exactly the same boundaries
and motif always wins. Loci that already have an Hg38LocusId (the lifted loci that the merge kept)
keep it.

A few known disease-associated loci still have no match after that, either because the liftover dropped their hg38
locus (eg. NOTCH2NLC, whose paralog NOTCH2NLA lifts onto the same CHM13 repeat, and XYLT1, whose 15bp hg38 repeat
lifts onto 256bp of CHM13) or because their CHM13 definition differs too much from the lifted one (eg. the primary
SAMD12 definition). Those that have a KnownDiseaseAssociatedLocus label get the LocusId that the known
disease-associated loci catalog's TRExplorerV2 field gives for that label, which is the hg38 TRExplorer catalog's
definition of the same disease locus.

It also adds an Hg38LocusMatchType field that compares the lifted hg38 locus's CHM13 boundaries with the locus's own
boundaries, the way str_analysis.merge_loci's FoundIn fields do:
- SameBoundaries: the lifted hg38 locus has exactly the locus's boundaries
- Hg38LocusIsWider or Hg38LocusIsNarrower: the lifted hg38 locus spans more or fewer CHM13 bases than the locus
- Hg38LocusIsShifted: the lifted hg38 locus spans the same number of CHM13 bases, at a different position
- SameKnownDiseaseAssociatedLocus: the Hg38LocusId comes from the locus's KnownDiseaseAssociatedLocus label rather
  than from a lifted hg38 locus
Loci without an hg38 counterpart, most of them in CHM13 satellite arrays that hg38 lacks, get None in both fields.
"""

import argparse
import collections
import gzip
import ijson
import intervaltree
import simplejson as json

from str_analysis.merge_loci import check_for_sufficient_overlap_and_motif_match
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure
from str_analysis.utils.misc_utils import parse_interval

HG38_LOCUS_ID_FIELD = "Hg38LocusId"
HG38_LOCUS_MATCH_TYPE_FIELD = "Hg38LocusMatchType"


def read_lifted_over_hg38_loci(lifted_over_hg38_catalog_json_path):
    """Reads the hg38 catalog's loci lifted over to CHM13.

    Args:
        lifted_over_hg38_catalog_json_path (str): path of the lifted over catalog, whose records have CHM13 coordinates
            in their ReferenceRegion and the hg38 LocusId in their Hg38LocusId field

    Returns:
        dict: maps each chromosome to an IntervalTree of the lifted loci. Each interval's data holds only the fields
            that matching needs (ReferenceRegion, LocusStructure and Hg38LocusId), to keep the memory use down.
    """
    intervals_by_chrom = collections.defaultdict(list)
    fopen = gzip.open if lifted_over_hg38_catalog_json_path.endswith("gz") else open
    with fopen(lifted_over_hg38_catalog_json_path, "rb") as f:
        for record in ijson.items(f, "item"):
            chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
            intervals_by_chrom[chrom].append(intervaltree.Interval(start_0based, end_1based, {
                "ReferenceRegion": record["ReferenceRegion"],
                "LocusStructure": record["LocusStructure"],
                HG38_LOCUS_ID_FIELD: record[HG38_LOCUS_ID_FIELD],
            }))

    return {chrom: intervaltree.IntervalTree(intervals) for chrom, intervals in intervals_by_chrom.items()}


def read_hg38_locus_ids_of_known_disease_associated_loci(known_disease_associated_loci_json_path):
    """Reads the hg38 TRExplorer catalog LocusId of each known disease-associated locus.

    Args:
        known_disease_associated_loci_json_path (str): path of the known disease-associated loci catalog whose LocusIds
            are the CHM13 catalog's KnownDiseaseAssociatedLocus labels, and whose TRExplorerV2 field holds the LocusId
            of the hg38 TRExplorer catalog's definition of each locus

    Returns:
        dict: maps the LocusId of each known disease-associated locus that has a TRExplorerV2 value to that value
    """
    fopen = gzip.open if known_disease_associated_loci_json_path.endswith("gz") else open
    with fopen(known_disease_associated_loci_json_path, "rt") as f:
        return {record["LocusId"]: record["TRExplorerV2"] for record in json.load(f) if record.get("TRExplorerV2")}


def get_canonical_motif(locus_structure):
    """Returns the canonical motif of a single-repeat LocusStructure, computed the way
    str_analysis.merge_loci.check_for_sufficient_overlap_and_motif_match computes it for step 11's merge, which doesn't
    use the reverse complement.

    Args:
        locus_structure (str): LocusStructure of a locus with one repeat, eg. (CAG)*

    Returns:
        str: the canonical motif
    """
    return compute_canonical_motif(parse_motifs_from_locus_structure(locus_structure)[0], include_reverse_complement=False)


def compute_jaccard_similarity(start1_0based, end1_1based, start2_0based, end2_1based):
    """Returns the Jaccard similarity of two intervals: the size of their intersection divided by the size of the
    interval that spans both, as str_analysis.merge_loci computes it.

    Args:
        start1_0based (int): start of the first interval
        end1_1based (int): end of the first interval
        start2_0based (int): start of the second interval
        end2_1based (int): end of the second interval

    Returns:
        float: the Jaccard similarity, between 0 and 1
    """
    intersection_size = max(0, min(end1_1based, end2_1based) - max(start1_0based, start2_0based))
    return intersection_size / (max(end1_1based, end2_1based) - min(start1_0based, start2_0based))


def get_hg38_locus_match_type(start_0based, end_1based, lifted_start_0based, lifted_end_1based):
    """Compares the CHM13 boundaries of a lifted hg38 locus with the boundaries of the locus it matched.

    Args:
        start_0based (int): start of the CHM13 catalog locus
        end_1based (int): end of the CHM13 catalog locus
        lifted_start_0based (int): CHM13 start of the lifted hg38 locus
        lifted_end_1based (int): CHM13 end of the lifted hg38 locus

    Returns:
        str: SameBoundaries, Hg38LocusIsWider, Hg38LocusIsNarrower or Hg38LocusIsShifted
    """
    if (lifted_start_0based, lifted_end_1based) == (start_0based, end_1based):
        return "SameBoundaries"

    size = end_1based - start_0based
    lifted_size = lifted_end_1based - lifted_start_0based
    if lifted_size > size:
        return "Hg38LocusIsWider"
    if lifted_size < size:
        return "Hg38LocusIsNarrower"
    return "Hg38LocusIsShifted"


def find_best_matching_lifted_hg38_locus(
        lifted_hg38_loci_by_chrom,
        record,
        min_overlap_fraction,
        min_jaccard_similarity,
        motif_length_match_sufficient_for_VNTRs):
    """Finds the lifted hg38 locus that best matches a CHM13 catalog locus.

    Args:
        lifted_hg38_loci_by_chrom (dict): lifted loci returned by read_lifted_over_hg38_loci
        record (dict): CHM13 catalog record with ReferenceRegion and LocusStructure fields
        min_overlap_fraction (float): see str_analysis.merge_loci's --overlap-fraction option
        min_jaccard_similarity (float): see str_analysis.merge_loci's --min-jaccard-similarity option
        motif_length_match_sufficient_for_VNTRs (bool): see str_analysis.merge_loci's
            --motif-length-match-sufficient-for-VNTRs option

    Returns:
        intervaltree.Interval: among the lifted loci that pass the merge's test, the one with the same canonical motif
            and the highest Jaccard similarity, with ties going to the smallest start, end and Hg38LocusId, or None if
            no lifted locus passes. A lifted locus with a different motif is only picked when none has the same motif.
    """
    chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
    if chrom not in lifted_hg38_loci_by_chrom:
        return None

    catalog_locus_interval = intervaltree.Interval(start_0based, end_1based, record)
    matching_lifted_intervals = [
        lifted_interval for lifted_interval in lifted_hg38_loci_by_chrom[chrom].overlap(start_0based, end_1based)
        if check_for_sufficient_overlap_and_motif_match(
            catalog_locus_interval,
            lifted_interval,
            min_overlap_fraction=min_overlap_fraction,
            min_jaccard_similarity=min_jaccard_similarity,
            motif_length_match_sufficient_for_VNTRs=motif_length_match_sufficient_for_VNTRs,
        ) is not None
    ]
    if not matching_lifted_intervals:
        return None

    # with motif_length_match_sufficient_for_VNTRs, a VNTR with a different motif of the same length also passes, so
    # rank the lifted loci with the same canonical motif first
    catalog_locus_canonical_motif = get_canonical_motif(record["LocusStructure"])
    return min(matching_lifted_intervals, key=lambda lifted_interval: (
        get_canonical_motif(lifted_interval.data["LocusStructure"]) != catalog_locus_canonical_motif,
        -compute_jaccard_similarity(start_0based, end_1based, lifted_interval.begin, lifted_interval.end),
        lifted_interval.begin,
        lifted_interval.end,
        lifted_interval.data[HG38_LOCUS_ID_FIELD]))


def find_lifted_hg38_locus_with_hg38_locus_id(lifted_hg38_loci_by_chrom, record):
    """Finds the lifted hg38 locus that a CHM13 catalog locus's existing Hg38LocusId refers to.

    Args:
        lifted_hg38_loci_by_chrom (dict): lifted loci returned by read_lifted_over_hg38_loci
        record (dict): CHM13 catalog record with ReferenceRegion and Hg38LocusId fields

    Returns:
        intervaltree.Interval: the overlapping lifted locus with the same Hg38LocusId, or None if there isn't one
    """
    chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
    for lifted_interval in lifted_hg38_loci_by_chrom.get(chrom, intervaltree.IntervalTree()).overlap(
            start_0based, end_1based):
        if lifted_interval.data[HG38_LOCUS_ID_FIELD] == record[HG38_LOCUS_ID_FIELD]:
            return lifted_interval

    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--lifted-over-hg38-catalog", required=True,
                        help="Path of the hg38 catalog lifted over to CHM13, as filtered for step 11's merge, whose "
                             "records have the hg38 LocusId in their Hg38LocusId field")
    parser.add_argument("--min-overlap-fraction", type=float, required=True,
                        help="The overlap fraction that step 11's merge used (str_analysis.merge_loci's "
                             "--overlap-fraction option)")
    parser.add_argument("--min-jaccard-similarity", type=float, required=True,
                        help="The minimum Jaccard similarity that step 11's merge used (str_analysis.merge_loci's "
                             "--min-jaccard-similarity option)")
    parser.add_argument("--motif-length-match-sufficient-for-VNTRs", action="store_true",
                        help="Set this if step 11's merge used str_analysis.merge_loci's option of the same name")
    parser.add_argument("--known-disease-associated-loci", required=True,
                        help="Path of the known disease-associated loci catalog that the KnownDiseaseAssociatedLocus "
                             "labels come from, whose TRExplorerV2 field holds hg38 TRExplorer catalog LocusIds")
    parser.add_argument("-o", "--output-catalog-json-path", required=True, help="Path of the output catalog JSON file")
    parser.add_argument("catalog_json_path", help="Path of the input CHM13 catalog JSON file")
    args = parser.parse_args()

    lifted_hg38_loci_by_chrom = read_lifted_over_hg38_loci(args.lifted_over_hg38_catalog)
    print(f"Read {sum(len(t) for t in lifted_hg38_loci_by_chrom.values()):,d} lifted hg38 loci from "
          f"{args.lifted_over_hg38_catalog}")

    hg38_locus_ids_by_known_disease_associated_locus = read_hg38_locus_ids_of_known_disease_associated_loci(
        args.known_disease_associated_loci)

    # stream the records one at a time, since loading a whole annotated catalog takes about 20GB of memory
    counters = collections.Counter()
    counts_by_match_type = collections.Counter()
    fopen_input = gzip.open if args.catalog_json_path.endswith("gz") else open
    fopen_output = gzip.open if args.output_catalog_json_path.endswith("gz") else open
    with fopen_input(args.catalog_json_path, "rb") as f, fopen_output(args.output_catalog_json_path, "wt") as f2:
        f2.write("[")
        for i, record in enumerate(ijson.items(f, "item")):
            _, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
            if record.get(HG38_LOCUS_ID_FIELD):
                lifted_interval = find_lifted_hg38_locus_with_hg38_locus_id(lifted_hg38_loci_by_chrom, record)
                counters["already had an Hg38LocusId"] += 1
                if lifted_interval is None:
                    counters["already had an Hg38LocusId that no overlapping lifted locus has"] += 1
                    record[HG38_LOCUS_MATCH_TYPE_FIELD] = None
                else:
                    record[HG38_LOCUS_MATCH_TYPE_FIELD] = get_hg38_locus_match_type(
                        start_0based, end_1based, lifted_interval.begin, lifted_interval.end)
            else:
                lifted_interval = find_best_matching_lifted_hg38_locus(
                    lifted_hg38_loci_by_chrom,
                    record,
                    min_overlap_fraction=args.min_overlap_fraction,
                    min_jaccard_similarity=args.min_jaccard_similarity,
                    motif_length_match_sufficient_for_VNTRs=args.motif_length_match_sufficient_for_VNTRs)
                if lifted_interval is not None:
                    record[HG38_LOCUS_ID_FIELD] = lifted_interval.data[HG38_LOCUS_ID_FIELD]
                    record[HG38_LOCUS_MATCH_TYPE_FIELD] = get_hg38_locus_match_type(
                        start_0based, end_1based, lifted_interval.begin, lifted_interval.end)
                    counters["got an Hg38LocusId from a matching lifted hg38 locus"] += 1
                elif record.get("KnownDiseaseAssociatedLocus") in hg38_locus_ids_by_known_disease_associated_locus:
                    record[HG38_LOCUS_ID_FIELD] = hg38_locus_ids_by_known_disease_associated_locus[
                        record["KnownDiseaseAssociatedLocus"]]
                    record[HG38_LOCUS_MATCH_TYPE_FIELD] = "SameKnownDiseaseAssociatedLocus"
                    counters["got an Hg38LocusId from its KnownDiseaseAssociatedLocus label"] += 1
                else:
                    record[HG38_LOCUS_ID_FIELD] = None
                    record[HG38_LOCUS_MATCH_TYPE_FIELD] = None
                    counters["have no hg38 counterpart"] += 1

            counts_by_match_type[record[HG38_LOCUS_MATCH_TYPE_FIELD]] += 1

            if i > 0:
                f2.write(", ")
            f2.write(json.dumps(record, use_decimal=True, indent=4))
        f2.write("]")

    print(f"Wrote {sum(counts_by_match_type.values()):,d} loci to {args.output_catalog_json_path}. Of these:")
    for description, count in sorted(counters.items()):
        print(f"{count:12,d}  {description}")
    print(f"{HG38_LOCUS_MATCH_TYPE_FIELD} counts:")
    for match_type, count in sorted(counts_by_match_type.items(), key=lambda x: -x[1]):
        print(f"{count:12,d}  {match_type}")


if __name__ == "__main__":
    main()
