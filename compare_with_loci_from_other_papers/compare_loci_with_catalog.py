"""This script takes a repeat catalog in BED format and determines how well its loci are represented in
the TRExplorer catalog (or some other main TR catalog like Adotto or Platinum).
"""

import argparse
import collections
import intervaltree
import gzip
import os
import pandas as pd
import re
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

NO_MATCH_FOUND = "no (no overlapping definitions)"

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-n", type=int, help="Number of loci to process")
    p.add_argument("--print-stats", type=int, default=1, choices=[0, 1, 2, 3],
                   help="At the end, output some stats at this level of detail. The higher the number, the more detailed the stats")
    p.add_argument("--write-loci-absent-from-new-catalog", action="store_true",
                   help="When generating the output table, include loci that are absent from the new catalog")
    p.add_argument("--catalog-bed-path",
                   default="TRExplorer_v2:../results__2025-11-03/release_draft_2025-11-03/repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz",
                   help="BED file path for the main TR catalog. Optionally, the path can be preceded by a name for the catalog, followed by ':' and then the path")
    p.add_argument("new_catalog", nargs="+",
                   help="BED file path for the new catalog to compare with locus definitions in the main catalog specified by --catalog-bed-path. "
                        "Optionally, the path can be preceded by a name for the catalog, followed by ':' and then the path")
    args = p.parse_args()

    if ":" in args.catalog_bed_path:
        if args.catalog_bed_path.count(":") > 1:
            p.error(f"More than one ':' found in {args.catalog_bed_path}")
        args.catalog_name, args.catalog_bed_path = args.catalog_bed_path.split(":")
    else:
        args.catalog_name = re.sub(".bed(.b?gz)?$", os.path.basename(args.catalog_bed_path))

    if not os.path.isfile(args.catalog_bed_path):
        p.error(f"File not found: {args.catalog_bed_path}")

    new_catalog_names = []
    for i, path in enumerate(args.new_catalog):
        if ":" in path:
            if path.count(":") > 1:
                p.error(f"More than one ':' found in {path}")
            new_catalog_name, path = path.split(":")
            args.new_catalog[i] = path
            new_catalog_names.append(new_catalog_name)
        else:
            new_catalog_names.append(re.sub(".bed(.b?gz)?$", os.path.basename(path)))

        if not os.path.isfile(path):
            p.error(f"File not found: {path}")


    OVERLAP_SCORE_MAP[0] = f"absent from {args.catalog_name} catalog"
    MOTIF_MATCH_SCORE_MAP[0] = f"absent from {args.catalog_name} catalog"
    OVERLAP_SCORE_MAP[-1] = f"only in {args.catalog_name} catalog"
    MOTIF_MATCH_SCORE_MAP[-1] = f"only in {args.catalog_name} catalog"

    main_catalog_loci = load_main_catalog_loci(args)

    for new_catalog_path, new_catalog_name in zip(args.new_catalog, new_catalog_names):
        output_tsv = new_catalog_path.replace(".bed", "").replace(".gz", "") + f".overlap_with_{args.catalog_name}.tsv.gz"
        print(f"Comparing {output_tsv}")
        df = compare_loci(
            main_catalog_loci,
            new_catalog_path,
            main_catalog_name=args.catalog_name,
            new_catalog_name=new_catalog_name,
            write_loci_absent_from_new_catalog=args.write_loci_absent_from_new_catalog)
        df.to_csv(output_tsv, sep="\t", index=False)
        print(f"Wrote {len(df):,d} rows to {output_tsv}")

        if args.print_stats > 0:
            print_stats(args, main_catalog_loci, df, new_catalog_name=new_catalog_name)

def print_stats(args, main_catalog_loci, df, new_catalog_name):
        # print some stats
        if args.print_stats >= 1:
            print("Summary:")
            for key, count in sorted(df["match_found?"].value_counts().items(), reverse=True):
                print(f"{count:10,d}  {key}")

            total_count = collections.Counter()
            for _, tree in main_catalog_loci.items():
                for interval in tree:
                    total_count["TRs"] += 1
                    canonical_motif = interval.data  # this motif was simplified before computing the canonical motif
                    if len(canonical_motif) <= 6:
                        total_count["STRs"] += 1
                    else:
                        total_count["VNTRs"] += 1

            for locus_type in "TRs", "STRs", "VNTRs":
                if locus_type == "TRs":
                    current_df = df
                elif locus_type == "STRs":
                    current_df = df[df[f"{new_catalog_name}_canonical_motif"].str.len() <= 6]  # this motif was simplified before computing the canonical motif
                elif locus_type == "VNTRs":
                    current_df = df[df[f"{new_catalog_name}_canonical_motif"].str.len() > 6]  # this motif was simplified before computing the canonical motif
                else:
                    raise ValueError(f"Unknown locus type: {locus_type}")

                if len(current_df) == 0 or total_count[locus_type] == 0:
                    continue

                print("---")
                print(f"Summary of {locus_type}:")
                yes_captured_count = sum(current_df["match_found?"].str.startswith("yes "))
                sort_of_captured_count = sum(current_df["match_found?"].str.startswith("sort of "))
                not_captured_count = sum(current_df["match_found?"].str.startswith("no "))

                margin = 4
                print(" "*margin, f"{yes_captured_count:10,d} of {new_catalog_name} {locus_type} are captured by {args.catalog_name}")
                print(" "*margin, f"{sort_of_captured_count:10,d} of additional {new_catalog_name} {locus_type} are sort of captured by {args.catalog_name}")
                print(" "*margin, f"{not_captured_count:10,d} of {new_catalog_name} {locus_type} are not in {args.catalog_name}")
                print(" "*margin, f"{len(current_df):10,d} total {locus_type} in {new_catalog_name}")
                print(" "*margin, f"{total_count[locus_type]:10,d} total {locus_type} in {args.catalog_name}")
                if locus_type != "TRs":
                    print(" "*margin, f"{total_count['TRs']:10,d} total TRs in {args.catalog_name}")
                print(" "*margin, f"{yes_captured_count/len(current_df):10.2f} recall            ( = {yes_captured_count:,d}/{len(current_df):,d} = fraction of {new_catalog_name} {locus_type} captured by {args.catalog_name})")
                print(" "*margin, f"{yes_captured_count/total_count[locus_type]:10.2e} precision         ( = {yes_captured_count:,d}/{total_count[locus_type]:,d} = {new_catalog_name} {locus_type} captured by {args.catalog_name} divided by total {locus_type} in {args.catalog_name})")
                print(" "*margin, f"{(yes_captured_count + sort_of_captured_count)/len(current_df):10.2f} sort of recall    ( = {yes_captured_count + sort_of_captured_count:,d}/{len(current_df):,d} = fraction of {new_catalog_name} {locus_type} sort of captured by {args.catalog_name})")
                print(" "*margin, f"{(yes_captured_count + sort_of_captured_count)/total_count[locus_type]:10.2e} sort of precision ( = {yes_captured_count + sort_of_captured_count:,d}/{total_count[locus_type]:,d} = {new_catalog_name} {locus_type} sort of captured by {args.catalog_name} divided by total {locus_type} in {args.catalog_name})")

        if args.print_stats >= 2:
            print()
            print("Details:")
            for (overlap_score, motif_match_score), group_df in sorted(df.groupby(["overlap_score", "motif_match_score"]), reverse=True):
                print(f"{len(group_df):10,d} loci:   {OVERLAP_SCORE_MAP[overlap_score]:<40}    {MOTIF_MATCH_SCORE_MAP[motif_match_score]:<20}")

        if args.print_stats >= 3:
            print()
            print(f"All loci in {new_catalog_name}:")
            for _, row in df[df[f"{new_catalog_name}_reference_region"].notna()].iterrows():
                print(" "*8, "\t".join(map(str, [row[f"{new_catalog_name}_reference_region"], row["match_found?"]])))


def load_main_catalog_loci(args):   
    print(f"Parsing {args.catalog_bed_path} to interval tree")
    main_catalog_loci = collections.defaultdict(intervaltree.IntervalTree)
    fopen = gzip.open if args.catalog_bed_path.endswith("gz") else open
    with fopen(args.catalog_bed_path, "rt") as f:
        for counter, line in enumerate(tqdm.tqdm(f, unit=" records", unit_scale=True)):
            if args.n is not None and counter >= args.n:
                break

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


def compute_match_summary(overlap_score, motif_match_score, new_catalog_motif, main_catalog_motif, main_catalog_name=None):
    high_overlap_score = overlap_score >= 5
    medium_or_high_overlap_score = overlap_score >= 2
    high_motif_similarity = (motif_match_score == 3 or (len(new_catalog_motif) > 6 and motif_match_score == 2))

    if main_catalog_motif is not None and not high_motif_similarity:
        #different_motifs_description = f"different motifs, "
        different_motifs_description = f"{len(main_catalog_motif)}bp motif in {main_catalog_name} instead of {len(new_catalog_motif)}bp, "
    else:
        different_motifs_description = ""

    if high_overlap_score and high_motif_similarity:
        was_match_found = "yes (Jaccard > 0.66 and similar motifs)"
    elif medium_or_high_overlap_score:
        if high_overlap_score:
            was_match_found = f"sort of ({different_motifs_description}Jaccard > 0.66)"
        else:
            was_match_found = f"sort of ({different_motifs_description}0.2 > Jaccard <= 0.66)"
    elif overlap_score > 0 and motif_match_score > 0:
        was_match_found = f"no ({different_motifs_description}Jaccard <= 0.2)"
    else:
        was_match_found = NO_MATCH_FOUND

    return was_match_found

def compare_loci(
        main_catalog_loci,
        new_catalog,
        main_catalog_name="catalog1",
        new_catalog_name="catalog2",
        write_loci_absent_from_new_catalog=False):
    """Compare a new catalog to the main catalog and output a TSV file of the results"""

    # parse new catalog and compare to main_catalog loci
    f = None

    if isinstance(new_catalog, (list, tuple)):
        fields_iterator = new_catalog
    elif isinstance(new_catalog, str):
        if os.path.isfile(new_catalog):
            fopen = gzip.open if new_catalog.endswith("gz") else open
            f = fopen(new_catalog, "rt")
            fields_iterator = (line.strip().split("\t") for line in f)
        else:
            raise ValueError(f"File not found: {new_catalog}")
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

        main_catalog_canonical_motif = closest_match_main_catalog_interval.data if closest_match_main_catalog_interval is not None else None
        was_match_found = compute_match_summary(overlap_score, motif_match_score, canonical_motif, main_catalog_canonical_motif, main_catalog_name)

        output_row = {
            "chrom": chrom,
            f"{new_catalog_name}_start_0based": start_0based,
            f"{new_catalog_name}_end_1based": end_1based,
            f"{new_catalog_name}_motif": motif,
            f"{new_catalog_name}_simplified_motif": simplified_motif,
            f"{new_catalog_name}_canonical_motif": canonical_motif,
            f"{new_catalog_name}_reference_region": f"{chrom}:{start_0based}-{end_1based}",
            f"{new_catalog_name}_reference_repeat_count": (end_1based - start_0based) // len(canonical_motif),
            f"{new_catalog_name}_reference_region_size": end_1based - start_0based,
            "overlap_score": overlap_score,
            "motif_match_score": motif_match_score,
            "overlap": OVERLAP_SCORE_MAP[overlap_score],
            "motif_match": MOTIF_MATCH_SCORE_MAP[motif_match_score],
            "match_found?": was_match_found,
        }

        if closest_match_main_catalog_interval is not None:
            closest_match_start_0based = closest_match_main_catalog_interval.begin
            closest_match_end_1based = closest_match_main_catalog_interval.end
            closest_match_canonical_motif = closest_match_main_catalog_interval.data
            closest_match_reference_region = f"{chrom}:{closest_match_start_0based}-{closest_match_end_1based}"
            output_row.update({
                f"{main_catalog_name}_start_0based": closest_match_start_0based,
                f"{main_catalog_name}_end_1based": closest_match_end_1based,
                f"{main_catalog_name}_canonical_motif": closest_match_canonical_motif,
                f"{main_catalog_name}_reference_region": closest_match_reference_region,
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

                reference_region = f"{chrom}:{main_catalog_interval.begin}-{main_catalog_interval.end}"
                output_rows.append({
                    "chrom": chrom,
                    "overlap_score": -1,
                    "motif_match_score": -1,
                    "overlap": OVERLAP_SCORE_MAP[-1],
                    "motif_match": MOTIF_MATCH_SCORE_MAP[-1],
                    "match_found?": NO_MATCH_FOUND,
                    f"{main_catalog_name}_start_0based": main_catalog_interval.begin,
                    f"{main_catalog_name}_end_1based": main_catalog_interval.end,
                    f"{main_catalog_name}_canonical_motif": main_catalog_canonical_motif,
                    f"{main_catalog_name}_reference_region": reference_region,
                    f"{main_catalog_name}_reference_repeat_count": (main_catalog_interval.end - main_catalog_interval.begin) // len(main_catalog_canonical_motif),
                    f"{main_catalog_name}_reference_region_size": main_catalog_interval.end - main_catalog_interval.begin,
                })

    if f is not None:
        f.close()
    
    return pd.DataFrame(output_rows)



if __name__ == "__main__":
    main()