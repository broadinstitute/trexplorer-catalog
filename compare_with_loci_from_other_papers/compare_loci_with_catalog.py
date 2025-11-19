"""This script takes a repeat catalog in BED format and determines how well its loci are represented in the TRExplorer catalog"""

import argparse
import collections
import intervaltree
import gzip
import os
import pandas as pd
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif



def main():
    p = argparse.ArgumentParser()
    p.add_argument("-n", type=int, help="Number of loci to process")
    p.add_argument("--catalog-name", default="trexplorer", help="Name of the catalog to compare to")
    p.add_argument("--catalog-bed-path", default="../results__2025-11-03/release_draft_2025-11-03/repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz")
    p.add_argument("new_catalog", nargs="+", help="Path of new catalog in BED format. If multiple paths are provided, they will be compared to the TRExplorer catalog in parallel.")
    args = p.parse_args()

    main_catalog_loci = load_main_catalog_loci(args)

    for new_catalog_path in args.new_catalog:
        output_tsv = new_catalog_path.replace(".bed", "").replace(".gz", "") + ".main_catalog_overlap.tsv.gz"
        df = compare_loci(main_catalog_loci, new_catalog_path, main_catalog_name=args.catalog_name)
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
            canonical_motif = compute_canonical_motif(fields[3])
            main_catalog_loci[chrom].add(intervaltree.Interval(start, end, data=canonical_motif))

    return main_catalog_loci



MOTIF_MATCH_SCORE_MAP = {
    3: "same motif",
    2: "same motif length",
    1: "different motif length",
    0: "absent",
}

def compute_motif_match_score(main_catalog_canonical_motif, other_catalog_canonical_motif):
    if main_catalog_canonical_motif == other_catalog_canonical_motif:
        return 3   # same motif
    elif len(main_catalog_canonical_motif) == len(other_catalog_canonical_motif):
        return 2   # same length
    else:
        return 1   # different length


OVERLAP_SCORE_MAP = {
    7: "exact match",
    6: "within 2 repeats",
    5: "union/overlap < 1.5x",
    4: "1.5x <= union/overlap < 2x",
    3: "2x <= union/overlap < 3x",
    2: "3x <= union/overlap < 4x",
    1: "union/overlap >= 4x",
    0: "absent",
}

def compute_overlap_score(main_catalog_interval, start_0based, end_1based, motif_size):
    if main_catalog_interval.begin == start_0based and main_catalog_interval.end == end_1based:
        return 7
    
    union_size = max(main_catalog_interval.end, end_1based) - min(main_catalog_interval.begin, start_0based)
    overlap_size = main_catalog_interval.overlap_size(start_0based, end_1based)
    if abs(overlap_size - union_size) <= 2*motif_size:
        return 6
    union_overlap_ratio = union_size / abs(overlap_size)
    if union_overlap_ratio < 1.5:
        return 5
    if union_overlap_ratio < 2.0:
        return 4
    if union_overlap_ratio < 3.0:
        return 3
    if union_overlap_ratio < 4.0:
        return 2
    return 1
    
    

def compare_loci(main_catalog_loci, new_catalog, main_catalog_name="trexplorer"):
    """Compare a new catalog to the main catalog and output a TSV file of the results"""

    # parse new catalog and compare to main_catalog loci
    f = None

    if isinstance(new_catalog, (list, tuple, iter)):
        fields_iterator = new_catalog
    elif isinstance(new_catalog, str) and os.path.isfile(new_catalog):
        f = gzip.open(new_catalog, "rt")
        fields_iterator = (line.strip().split("\t") for line in f)
    else:
        raise ValueError(f"Invalid new_catalog arg type: {type(new_catalog)}: {new_catalog}")
    
    output_rows = []
    for fields in tqdm.tqdm(fields_iterator, unit=" records", unit_scale=True):
        chrom = fields[0].replace("chr", "")
        start_0based = int(fields[1])
        end_1based = int(fields[2])
        motif = fields[3]
        canonical_motif = compute_canonical_motif(motif)

        overlap_score = 0
        motif_match_score = 0
        closest_match_main_catalog_interval = None

        # find the overlapping main catalog locus that overlaps the largest percentage of the new locus and is closest in size
        main_catalog_overlap = main_catalog_loci[chrom].overlap(start_0based, end_1based)
        for main_catalog_interval in main_catalog_overlap:
            main_catalog_canonical_motif = main_catalog_interval.data

            current_overlap_score = compute_overlap_score(main_catalog_interval, start_0based, end_1based, len(canonical_motif))
            if current_overlap_score < overlap_score:
                continue

            current_motif_match_score = compute_motif_match_score(main_catalog_canonical_motif, canonical_motif)
            if (current_overlap_score, current_motif_match_score) < (overlap_score, motif_match_score):
                continue

            overlap_score = current_overlap_score
            motif_match_score = current_motif_match_score
            closest_match_main_catalog_interval = main_catalog_interval

        output_rows.append({
            "chrom": chrom,
            "start_0based": start_0based,
            "end_1based": end_1based,
            "motif": motif,
            "canonical_motif": canonical_motif,
            "repeat_count": (end_1based - start_0based) // len(canonical_motif),
            "reference_region_size": end_1based - start_0based,
            "overlap_score": overlap_score,
            "motif_match_score": motif_match_score,
            "overlap": OVERLAP_SCORE_MAP[overlap_score],
            "motif_match": MOTIF_MATCH_SCORE_MAP[motif_match_score],
            f"{main_catalog_name}_start_0based": closest_match_main_catalog_interval.begin if closest_match_main_catalog_interval is not None else None,
            f"{main_catalog_name}_end_1based": closest_match_main_catalog_interval.end if closest_match_main_catalog_interval is not None else None,
            f"{main_catalog_name}_canonical_motif": closest_match_main_catalog_interval.data if closest_match_main_catalog_interval is not None else None,
        })

    if f is not None:
        f.close()
    
    return pd.DataFrame(output_rows)



if __name__ == "__main__":
    main()