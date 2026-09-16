"""Compares the locus ids in a released TRGT catalog BED file against the locus ids that actually
appear in a genotyped output file (an HPRC-style per-locus-and-motif LPS TSV, or a TRGT VCF).

This is a manual, run-by-hand diagnostic, not a pipeline gate: the two files can legitimately
disagree, since the variation clusters table that builds the catalog's VC rows is computed against
a different set of ids than whatever catalog originally produced the genotyped output
(scripts/generate_TRGT_catalog_with_variation_clusters_and_all_repeats.py documents the same
caveat). Run this after a release to see whether that gap looks like the expected background level
or like the kind of catalog/output mismatch reported in
https://github.com/PacificBiosciences/trgt-lps/issues/5#issuecomment-5683540748, where an entire
class of loci (variation cluster members) was missing from the catalog.
"""

import argparse
import collections
import gzip
import re


def open_maybe_gzipped(path):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def read_catalog_locus_ids(path, chrom=None):
    locus_ids = set()
    with open_maybe_gzipped(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if chrom and fields[0] != chrom:
                continue
            locus_ids.add(re.search(r"ID=([^;]+)", fields[3]).group(1))
    return locus_ids


def read_vcf_locus_ids(path, chrom=None):
    locus_ids = set()
    with open_maybe_gzipped(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t", 8)
            if chrom and fields[0] != chrom:
                continue
            match = re.search(r"TRID=([^;]+)", fields[7])
            if match:
                locus_ids.add(match.group(1))
    return locus_ids


def read_lps_tsv_locus_ids(path, chrom=None):
    locus_ids = set()
    with open_maybe_gzipped(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        locus_id_column = header.index("locus_id")
        for line in f:
            locus_id = line.rstrip("\n").split("\t", locus_id_column + 1)[locus_id_column]
            if chrom and not (locus_id.startswith(f"{chrom.replace('chr', '')}-")
                               or locus_id.startswith(f"VC:{chrom.replace('chr', '')}:")):
                continue
            locus_ids.add(locus_id)
    return locus_ids


def read_genotyped_locus_ids(path, chrom=None):
    if ".vcf" in path:
        return read_vcf_locus_ids(path, chrom=chrom)
    return read_lps_tsv_locus_ids(path, chrom=chrom)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--chrom", help="Restrict the comparison to one chromosome, e.g. chr19")
    parser.add_argument("--show-examples", type=int, default=10,
                        help="Number of example locus ids to print per category")
    parser.add_argument("catalog_bed_path", help="Path of the released TRGT catalog BED file")
    parser.add_argument("genotyped_path",
                        help="Path of the genotyped output to compare against: an HPRC-style "
                             "per-locus-and-motif LPS TSV, or a TRGT VCF")
    args = parser.parse_args()

    catalog_ids = read_catalog_locus_ids(args.catalog_bed_path, chrom=args.chrom)
    genotyped_ids = read_genotyped_locus_ids(args.genotyped_path, chrom=args.chrom)

    catalog_only = catalog_ids - genotyped_ids
    genotyped_only = genotyped_ids - catalog_ids
    in_both = catalog_ids & genotyped_ids

    genotyped_only_by_kind = collections.Counter(
        "variation cluster" if locus_id.startswith("VC:") else "repeat" for locus_id in genotyped_only)
    catalog_only_by_kind = collections.Counter(
        "variation cluster" if locus_id.startswith("VC:") else "repeat" for locus_id in catalog_only)

    print(f"{len(catalog_ids):,d} locus ids in {args.catalog_bed_path}")
    print(f"{len(genotyped_ids):,d} locus ids in {args.genotyped_path}")
    print(f"{len(in_both):,d} locus ids in both")
    print(f"{len(catalog_only):,d} locus ids only in the catalog: {dict(catalog_only_by_kind)}")
    if catalog_only:
        print("  examples:", list(catalog_only)[:args.show_examples])
    print(f"{len(genotyped_only):,d} locus ids only in the genotyped output: "
          f"{dict(genotyped_only_by_kind)}")
    if genotyped_only:
        print("  examples:", list(genotyped_only)[:args.show_examples])


if __name__ == "__main__":
    main()
