"""Generate a separate BED file for each Source in the annotated EH catalog JSON.

Each BED file contains the intervals (from ReferenceRegion) and motif (from LocusStructure).
Output files are bgzipped and tabix-indexed.
"""

import argparse
import collections
import gzip
import os
import subprocess

import ijson


DEFAULT_CATALOG_PATH = "results__2026-02-01/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz"


def parse_motif(locus_structure):
    """Extract motif from LocusStructure like '(TAACCC)*'."""
    if locus_structure.startswith("(") and ")*" in locus_structure:
        return locus_structure[1:locus_structure.index(")*")]
    return locus_structure


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--annotated-catalog-json", default=DEFAULT_CATALOG_PATH,
                        help="Path to the annotated EH catalog JSON.gz file")
    parser.add_argument("--output-dir", default="bed_file_per_source_catalog",
                        help="Output directory for per-source BED files")
    args = parser.parse_args()

    if not os.path.isfile(args.annotated_catalog_json):
        parser.error(f"File not found: {args.annotated_catalog_json}")

    os.makedirs(args.output_dir, exist_ok=True)

    # First pass: write unsorted BED lines per source
    file_handles = {}
    counts = collections.Counter()
    print(f"Reading {args.annotated_catalog_json}...")
    with gzip.open(args.annotated_catalog_json, "rt") as f:
        for record in ijson.items(f, "item", use_float=True):
            source = record.get("Source", "UNKNOWN")
            region = record["ReferenceRegion"]
            chrom, coords = region.split(":")
            start, end = coords.split("-")
            motif = parse_motif(record["LocusStructure"])

            if source not in file_handles:
                bed_path = os.path.join(args.output_dir, f"{source}.bed")
                file_handles[source] = open(bed_path, "w")

            file_handles[source].write(f"{chrom}\t{start}\t{end}\t{motif}\n")
            counts[source] += 1

    for fh in file_handles.values():
        fh.close()

    # bgzip and tabix each BED file
    for source in sorted(file_handles):
        bed_path = os.path.join(args.output_dir, f"{source}.bed")
        gz_path = bed_path + ".gz"
        subprocess.run(["bgzip", "-f", bed_path], check=True)
        subprocess.run(["tabix", "-p", "bed", gz_path], check=True)

    # Print stats
    print()
    total = sum(counts.values())
    print(f"{'Source':<50s} {'Loci':>10s}  {'%':>6s}")
    print("-" * 70)
    for source in sorted(counts, key=lambda s: -counts[s]):
        c = counts[source]
        print(f"{source:<50s} {c:>10,d}  {100*c/total:5.1f}%")
    print("-" * 70)
    print(f"{'TOTAL':<50s} {total:>10,d}")
    print(f"\nWrote {len(counts)} BED files to {args.output_dir}/")


if __name__ == "__main__":
    main()
