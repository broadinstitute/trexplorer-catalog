"""This script takes a repeat catalog in BED format and determines how well its loci are represented in the TRExplorer catalog"""

import argparse
import collections
import intervaltree
import gzip
import os
import pandas as pd
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.find_repeat_unit import find_repeat_unit_without_allowing_interruptions

MOTIF_MATCH_SCORE_MAP = {
    3: "same motif",
    2: "same motif length",
    1: "different motif length",
    0: "absent from main catalog",
    -1: "absent from new catalog",
}


OVERLAP_SCORE_MAP = {
    7: "exact match",
    6: "diff ≤ 2 repeats",
    5: "Jaccard similarity > 0.66",
    4: "0.5 < Jaccard similarity ≤ 0.66",
    3: "0.33 < Jaccard similarity ≤ 0.5",
    2: "0.2 < Jaccard similarity ≤ 0.33",
    1: "Jaccard similarity ≤ 0.2",
    0: "absent from main catalog",
    -1: "absent from new catalog",
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-n", type=int, help="Number of loci to process")
    p.add_argument("--write-loci-absent-from-new-catalog", action="store_true", help="When generating the output table, include loci that are absent from the new catalog")
    p.add_argument("--catalog-name", default="TRExplorer_v2", help="Name of the catalog to compare to")
    p.add_argument("--catalog-bed-path", default="../results__2025-11-03/release_draft_2025-11-03/repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz")
    p.add_argument("new_catalog", nargs="+", help="Path of new catalog in BED format. If multiple paths are provided, they will be compared to the TRExplorer catalog in parallel.")
    args = p.parse_args()

    OVERLAP_SCORE_MAP[0] = f"absent from {args.catalog_name} catalog"
    MOTIF_MATCH_SCORE_MAP[0] = f"absent from {args.catalog_name} catalog"
    OVERLAP_SCORE_MAP[-1] = f"only in {args.catalog_name} catalog"
    MOTIF_MATCH_SCORE_MAP[-1] = f"only in {args.catalog_name} catalog"

    main_catalog_loci = load_main_catalog_loci(args)

    for new_catalog_path in args.new_catalog:
        output_tsv = new_catalog_path.replace(".bed", "").replace(".gz", "") + f".overlap_with_{args.catalog_name}.tsv.gz"
        print(f"Comparing {output_tsv}")
        df = compare_loci(main_catalog_loci, new_catalog_path, main_catalog_name=args.catalog_name, write_loci_absent_from_new_catalog=args.write_loci_absent_from_new_catalog)
        df.to_csv(output_tsv, sep="\t", index=False)
        print(f"Wrote {len(df):,d} rows to {output_tsv}")


def load_main_catalog_loci(args):   
    print(f"Parsing {args.catalog_bed_path} to interval tree")
    main_catalog_loci = collections.defaultdict(intervaltree.IntervalTree)
    with gzip.open(args.catalog_bed_path, "rt") as f:
        for counter, line in enumerate(tqdm.tqdm(f, unit=" records", unit_scale=True)):
            if args.n is not None and counter >= args.n:
                break
            counter += 1
            fields = line.strip().split("\t")
            chrom = fields[0].replace("chr", "")
            start = int(fields[1])
            end = int(fields[2])
            motif = fields[3]
            simplified_motif, _, _ = find_repeat_unit_without_allowing_interruptions(motif, allow_partial_repeats=False)
            canonical_motif = compute_canonical_motif(simplified_motif)
            main_catalog_loci[chrom].add(intervaltree.Interval(start, end, data=canonical_motif))

    return main_catalog_loci


def compute_motif_match_score(main_catalog_canonical_motif, other_catalog_canonical_motif):
    if main_catalog_canonical_motif == other_catalog_canonical_motif:
        return 3   # same motif
    elif len(main_catalog_canonical_motif) == len(other_catalog_canonical_motif):
        return 2   # same length
    else:
        return 1   # different length


def compute_overlap_score(main_catalog_interval, start_0based, end_1based, min_motif_size):
    """Computes a numerical score for the overlap between the two intervals

    Args:
        main_catalog_interval (intervaltree.Interval): The interval from the main catalog that overlaps the new locus
        start_0based (int): The start of the new locus (0-based)
        end_1based (int): The end of the new locus (1-based)
        min_motif_size (int): The motif size at this locus, or the smaller of the two motif sizes if they differ between the two locus definitions

    Returns:
        int: A numerical score for the overlap between the two intervals
    """
    if not main_catalog_interval.overlaps(start_0based, end_1based):
        raise ValueError(f"Main catalog interval {main_catalog_interval} does not overlap with new locus {start_0based}-{end_1based}")
        # return 0 

    if main_catalog_interval.begin == start_0based and main_catalog_interval.end == end_1based:
        return 7
    
    union_size = max(main_catalog_interval.end, end_1based) - min(main_catalog_interval.begin, start_0based)
    intersection_size = main_catalog_interval.overlap_size(start_0based, end_1based)
    if abs(intersection_size - union_size) <= 2*min_motif_size:
        return 6
    jaccard_similarity = intersection_size / union_size
    if jaccard_similarity > 2/3.0:
        return 5
    if jaccard_similarity > 0.5:
        return 4
    if jaccard_similarity > 1/3.0:
        return 3
    if jaccard_similarity > 0.2:
        return 2

    return 1
    
    

def compare_loci(main_catalog_loci, new_catalog, main_catalog_name="trexplorer", write_loci_absent_from_new_catalog=False):
    """Compare a new catalog to the main catalog and output a TSV file of the results"""

    # parse new catalog and compare to main_catalog loci
    f = None

    if isinstance(new_catalog, (list, tuple)):
        fields_iterator = new_catalog
    elif isinstance(new_catalog, str) and os.path.isfile(new_catalog):
        f = gzip.open(new_catalog, "rt")
        fields_iterator = (line.strip().split("\t") for line in f)
    else:
        raise ValueError(f"Invalid new_catalog arg type: {type(new_catalog)}: {new_catalog}")

    output_rows = []
    main_catalog_loci_that_overlap_new_catalog = set()
    for fields in tqdm.tqdm(fields_iterator, unit=" records", unit_scale=True):
        chrom = fields[0].replace("chr", "")
        start_0based = int(fields[1])
        end_1based = int(fields[2])
        motif = fields[3]
        simplified_motif, _, _ = find_repeat_unit_without_allowing_interruptions(motif, allow_partial_repeats=False)
        canonical_motif = compute_canonical_motif(simplified_motif)

        overlap_score = 0
        motif_match_score = 0
        closest_match_main_catalog_interval = None

        # find the overlapping main catalog locus that overlaps the largest percentage of the new locus and is closest in size
        main_catalog_overlap = main_catalog_loci[chrom].overlap(start_0based, end_1based)
        for main_catalog_interval in main_catalog_overlap:
            main_catalog_canonical_motif = main_catalog_interval.data
            if write_loci_absent_from_new_catalog:
                main_catalog_loci_that_overlap_new_catalog.add((chrom, main_catalog_interval.begin, main_catalog_interval.end, main_catalog_canonical_motif))

            current_overlap_score = compute_overlap_score(
                main_catalog_interval,
                start_0based,
                end_1based,
                min_motif_size=min(len(canonical_motif), len(main_catalog_canonical_motif)))

            if current_overlap_score < overlap_score:
                continue

            current_motif_match_score = compute_motif_match_score(main_catalog_canonical_motif, canonical_motif)
            if (current_overlap_score, current_motif_match_score) < (overlap_score, motif_match_score):
                continue

            overlap_score = current_overlap_score
            motif_match_score = current_motif_match_score
            closest_match_main_catalog_interval = main_catalog_interval

        output_row = {
            "chrom": chrom,
            "start_0based": start_0based,
            "end_1based": end_1based,
            "motif": motif,
            "canonical_motif": canonical_motif,
            "reference_repeat_count": (end_1based - start_0based) // len(canonical_motif),
            "reference_region_size": end_1based - start_0based,
            "overlap_score": overlap_score,
            "motif_match_score": motif_match_score,
            "overlap": OVERLAP_SCORE_MAP[overlap_score],
            "motif_match": MOTIF_MATCH_SCORE_MAP[motif_match_score],
        }

        if closest_match_main_catalog_interval is not None:
            closest_match_start_0based = closest_match_main_catalog_interval.begin
            closest_match_end_1based = closest_match_main_catalog_interval.end
            closest_match_canonical_motif = closest_match_main_catalog_interval.data    
            output_row.update({
                f"{main_catalog_name}_start_0based": closest_match_start_0based,
                f"{main_catalog_name}_end_1based": closest_match_end_1based,
                f"{main_catalog_name}_canonical_motif": closest_match_canonical_motif,
                f"{main_catalog_name}_reference_repeat_count": (closest_match_end_1based - closest_match_start_0based) // len(closest_match_canonical_motif),
                f"{main_catalog_name}_reference_region_size": closest_match_end_1based - closest_match_start_0based,
            })

        output_rows.append(output_row)

    if write_loci_absent_from_new_catalog:
        for chrom, interval_tree in tqdm.tqdm(main_catalog_loci.items(), unit=" chromosomes", unit_scale=True):
            for main_catalog_interval in interval_tree:
                main_catalog_canonical_motif = main_catalog_interval.data
                if (chrom, main_catalog_interval.begin, main_catalog_interval.end, main_catalog_canonical_motif) in main_catalog_loci_that_overlap_new_catalog:
                    continue
                
                output_rows.append({
                    "chrom": chrom,
                    "start_0based": None,
                    "end_1based": None,
                    "motif": None,
                    "canonical_motif": None,
                    "reference_repeat_count": None,
                    "reference_region_size": None,
                    "overlap_score": -1,
                    "motif_match_score": -1,
                    "overlap": OVERLAP_SCORE_MAP[-1],
                    "motif_match": MOTIF_MATCH_SCORE_MAP[-1],
                    f"{main_catalog_name}_start_0based": main_catalog_interval.begin,
                    f"{main_catalog_name}_end_1based": main_catalog_interval.end,
                    f"{main_catalog_name}_canonical_motif": main_catalog_canonical_motif,
                    f"{main_catalog_name}_reference_repeat_count": (main_catalog_interval.end - main_catalog_interval.begin) // len(main_catalog_canonical_motif),
                    f"{main_catalog_name}_reference_region_size": main_catalog_interval.end - main_catalog_interval.begin,
                })

    if f is not None:
        f.close()
    
    return pd.DataFrame(output_rows)



if __name__ == "__main__":
    main()