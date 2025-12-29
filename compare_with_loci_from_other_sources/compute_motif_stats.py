"""Evaluate various properties of the motifs specified in a catalog, including:

- Motif size distribution
- (Locus size % Motif size) distribution
- Motif vs. the motif that maximizes purity of the reference sequence
- Motif vs. reference sequence TRF annotation
- Motif vs. motif definition at the exact same locus in other catalogs

These properties can be stratified by motif length and reference repeat count.
"""

import argparse
import collections
import intervaltree
import os
import pandas as pd
import pyfaidx
import re
import tqdm

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator
from str_analysis.utils.find_motif_utils import find_highest_purity_motif_length_for_interval, compute_motif_length_purity_for_interval, compute_motif_purity_for_interval


def compute_purity_stats_for_interval(reference_fasta, chrom, start_0based, end_1based, motif):
    motif_length_purity, most_common_motif = compute_motif_length_purity_for_interval(
        reference_fasta, chrom, start_0based, end_1based, len(motif))
    fraction_pure_bases, distance = compute_motif_purity_for_interval(
        reference_fasta, chrom, start_0based, end_1based, motif)

    return fraction_pure_bases, motif_length_purity, most_common_motif

"""
# Motif  Length Purity Compared To Optimal

IF [Motif Length Purity] > [Optimal Motif Purity] + 15 THEN
    "MOTIF MUCH PURITY HIGHER THAN OPTIMAL"
ELSEIF [Motif Length Purity] > [Optimal Motif Purity] THEN
    "MOTIF PURITY HIGHER THAN OPTIMAL"
ELSEIF [Motif Length Purity] == [Optimal Motif Purity] THEN
    "SAME MOTIF PURITY"
ELSEIF [Motif Length Purity] < [Optimal Motif Purity] - 0.15 THEN
    "MOTIF PURITY IS A LOT LOWER"
ELSE
    "MOTIF PURITY IS LOWER"
END

# Motif Length Compared To Optimal
IF [Reference Repeat Count] < 2 THEN
   "Reference Repeat Count < 2"
ELSEIF [Optimal Motif Reference Repeat Count] < 2 THEN
   "Optimal Motif Reference Repeat Count < 2"
ELSEIF [Motif Length] > [Optimal Motif Length] THEN
   IF [Motif Length] % [Optimal Motif Length] == 0 THEN
      "Larger Multiple Of Optimal Motif Length"
   ELSE
      "Longer Motif Length"
   END
ELSEIF [Motif Length] == [Optimal Motif Length] THEN
   "SAME MOTIF LENGTH"
ELSEIF [Motif Length] < [Optimal Motif Length] THEN
   IF [Optimal Motif Length] % [Motif Length] == 0 THEN
      "Smaller Factor Of Optimal Motif Length"
   ELSE
      "Shorter Motif Length"
   END
END


# Motif purity compared to optimal
IF [Motif Purity] > [Optimal Motif Purity] + 15 THEN
    "MOTIF MUCH PURITY HIGHER THAN OPTIMAL"
ELSEIF [Motif Purity] > [Optimal Motif Purity] THEN
    "MOTIF PURITY HIGHER THAN OPTIMAL"
ELSEIF [Motif Purity] == [Optimal Motif Purity] THEN
    "SAME MOTIF PURITY"
ELSEIF [Motif Purity] < [Optimal Motif Purity] - 0.15 THEN
    "MOTIF PURITY IS A LOT LOWER"
ELSE
    "MOTIF PURITY IS LOWER"
END
"""
def process_record(record, reference_fasta=None, verbose=False):
    motif_length = len(record["Motif"])
    chrom = record["Chrom"]
    start_0based = record["Start_0based"]
    end_1based = record["End_1based"]
    locus_width = end_1based - start_0based
    motif_purity, motif_length_purity, _ = compute_purity_stats_for_interval(reference_fasta, chrom, start_0based, end_1based, record["Motif"])
    optimal_motif, optimal_motif_purity, optimal_motif_quality_score = find_highest_purity_motif_length_for_interval(
        reference_fasta, chrom, start_0based, end_1based, max_motif_length=max(
            motif_length, min(motif_length * 2 + 1, (end_1based - start_0based - 1))),
        negligible_change_in_purity=0.02,
    )

    canonical_motif = compute_canonical_motif(record["Motif"])
    canonical_optimal_motif = compute_canonical_motif(optimal_motif)
    optimal_motif_length = len(optimal_motif)
    record.update({
        "LocusWidth": locus_width,
        "IsLocusWidthExactMultipleOfMotif": (locus_width >= motif_length) and (locus_width % motif_length == 0),

        "CanonicalMotif": canonical_motif,
        "MotifLength": motif_length,
        "MotifRepeatCount": locus_width / motif_length,

        "MotifPurity": motif_purity,
        "MotifLengthPurity": motif_length_purity,

        "MotifVsOptimalMotif": compute_motif_match(canonical_motif, canonical_optimal_motif),
        "OptimalMotif": optimal_motif,
        "OptimalCanonicalMotif": canonical_optimal_motif,
        "OptimalMotifLength": optimal_motif_length,
        "OptimalMotifRepeatCount": locus_width / optimal_motif_length,

        "OptimalMotifPurity": optimal_motif_purity,
        "OptimalMotifQualityScore": optimal_motif_quality_score,
    })

    if verbose:
        print(record)


def compute_motif_match(canonical_motif1, canonical_motif2):
    if canonical_motif1 == canonical_motif2:
        return "same motif"
    elif len(canonical_motif1) == len(canonical_motif2):
        return "same motif length"
    elif len(canonical_motif1) < len(canonical_motif2):
        if len(canonical_motif2) % len(canonical_motif1) == 0:
            return "motif2 is longer multiple of motif1"
        else:
            return "motif2 is longer"
    elif len(canonical_motif1) > len(canonical_motif2):
        if len(canonical_motif1) % len(canonical_motif2) == 0:
            return "motif1 is longer multiple of motif2"
        else:
            return "motif1 is longer"
    else:
        raise Exception(f"Logic error")


def compute_overlap_score(interval1, start_0based, end_1based, min_motif_size):
    if not interval1.overlaps(start_0based, end_1based):
        raise ValueError(f"Interval {interval1} does not overlap with new locus {start_0based}-{end_1based}")

    if interval1.begin == start_0based and interval1.end == end_1based:
        return "exact match"

    intersection_size = interval1.overlap_size(start_0based, end_1based)
    union_size = max(interval1.end, end_1based) - min(interval1.begin, start_0based)
    if abs(intersection_size - union_size) <= 2*min_motif_size:
        return "diff <= 2 repeats"

    jaccard_similarity = intersection_size / union_size
    if jaccard_similarity > 2/3.0:
        return "Jaccard > 0.66"
    elif jaccard_similarity > 0.2:
        return "0.2 < Jaccard <= 0.66"
    else:
        return "Jaccard < 0.2"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-fasta", help="Reference FASTA file path", default="~/hg38.fa")
    parser.add_argument("-o", "--output-prefix", help="Output prefix")
    parser.add_argument("--min-motif-size", type=int, help="Minimum motif size", default=1)
    parser.add_argument("--max-motif-size", type=int, help="Maximum motif size", default=1_000)
    parser.add_argument("-n", "--number-of-loci", type=int, help="Number of loci to process")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    parser.add_argument("--known-loci", action="append", help="Optional path of BED file with known loci. "
                        "Input catalog loci that overlap these intervals with Jaccard > 0.66 will have "
                        "'KnownLocus' set to True")
    parser.add_argument("catalog_path", help="Catalog BED or JSON file path")
    args = parser.parse_args()

    if not os.path.isfile(args.catalog_path):
        parser.error(f"File not found: {args.catalog_path}")
    
    catalog_iterator = get_variant_catalog_iterator(args.catalog_path)

    reference_fasta = pyfaidx.Fasta(os.path.expanduser(args.reference_fasta), as_raw=True, one_based_attributes=False, sequence_always_upper=True)

    if not args.output_prefix:
        args.output_prefix = re.sub("(.bed|.json)(.b?gz)?$", "", os.path.basename(args.catalog_path))

    # parse any specified known locus intervals
    known_loci_interval_trees = collections.defaultdict(intervaltree.IntervalTree)
    if args.known_loci:
        for known_loci_path in args.known_loci:
            if not os.path.isfile(known_loci_path):
                parser.error(f"File not found: {known_loci_path}")

        for known_loci_path in args.known_loci:
            for record in get_variant_catalog_iterator(known_loci_path):
                chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
                chrom = chrom.replace("chr", "")
                known_locus_motif = record["LocusStructure"].strip("()+*")
                known_loci_interval_trees[chrom].add(intervaltree.Interval(start_0based, end_1based, data={
                    "CanonicalMotif": compute_canonical_motif(known_locus_motif),
                    "LocusId": record["LocusId"],
                }))

    results = []
    for i, record in tqdm.tqdm(enumerate(catalog_iterator), unit=" loci", unit_scale=True):
        if args.number_of_loci is not None and i >= args.number_of_loci:
            break

        chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
        chrom = chrom.replace("chr", "")
        record["Chrom"] = chrom
        record["Start_0based"] = start_0based
        record["End_1based"] = end_1based
        if start_0based >= end_1based:        
            continue

        record["Motif"] = record["LocusStructure"].strip("()*+")
        motif_length = len(record["Motif"])
        if motif_length < args.min_motif_size or motif_length > args.max_motif_size:
            continue

        if args.known_loci:
            for interval in known_loci_interval_trees[chrom].overlap(start_0based, end_1based):
                # check if Jaccard similarity is > 0.66
                intersection_size = min(interval.end, end_1based) - max(interval.begin, start_0based)
                union_size = max(interval.end, end_1based) - min(interval.begin, start_0based)
                jaccard_similarity = intersection_size / union_size
                if jaccard_similarity > 0.66:
                    record["KnownLocusOverlap"] = compute_overlap_score(
                        interval, start_0based, end_1based, min(len(interval.data["CanonicalMotif"]), len(record["Motif"])))
                    record["KnownLocusMotifMatch"] = compute_motif_match(
                        interval.data["CanonicalMotif"], compute_canonical_motif(record["Motif"]))
                    break

        # Example record: {'LocusId': '1-50537-50545-CTCC', 'ReferenceRegion': '1:50537-50545', 'LocusStructure': '(CTCC)*', 'VariantType': 'Repeat', 'Motif': 'CTCC', 'CanonicalMotif': 'AGGG'}
        process_record(record, reference_fasta=reference_fasta, verbose=args.verbose)
        results.append(record)

    output_path = f"motif_stats.{args.output_prefix}.tsv"
    print(f"Writing results to {output_path}")
    df = pd.DataFrame(results)
    df.to_csv(output_path, sep="\t", index=False, header=True)
    print(f"Wrote {len(df):,d} rows to {output_path}")

if __name__ == "__main__":
    main()
