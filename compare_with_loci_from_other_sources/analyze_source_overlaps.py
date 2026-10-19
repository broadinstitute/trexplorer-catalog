"""Analyze self-overlap within each v2 source catalog and overlap with the merged v1 and v2a catalogs.

Two loci are considered overlapping if they overlap by at least 2x the larger motif size
(matching the criterion used in str_analysis.compute_catalog_stats).
"""

import argparse
import collections
import gzip
import json
import os

from ncls import NCLS
import numpy as np


# v2a sources
V2A_SOURCES = [
    ("KnownDiseaseAssociatedLociV2", "known_disease_associated_loci_v2.loci_to_include_in_catalog.filtered.json.gz"),
    ("PolymorphicTRsInT2TAssembliesV2", "combined.321_catalogs.tandem_repeats.filtered.json.gz"),
    ("KnownFunctionalVNTRs", "known_functional_VNTRs.loci_to_include_in_catalog.filtered.json.gz"),
]

# v2b sources
V2B_SOURCES = [
    ("HipSTRCatalog", "hg38.hipstr_reference.loci_to_include_in_catalog.filtered.json.gz"),
    ("AdottoTRsFromDanzi2025", "Danzi_2025.adotto.loci_to_include_in_catalog.filtered.json.gz"),
    ("ClinvarIndelsThatAreTRs2025", "clinvar_2025_11_03.loci_to_include_in_catalog.filtered.json.gz"),
    ("Tanudisastro2025", "Tanudisastro_2025.loci_to_include_in_catalog.filtered.json.gz"),
    ("Manigbas2024", "Manigbas_2024_genotyped_TRs.filtered.json.gz"),
    ("Sulovari2019", "Sulovari_2019_human_specific_STRs.filtered.json.gz"),
    ("Garg2021", "Garg_2021_eVNTRs_mVNTRs.filtered.json.gz"),
    ("Mukamel2021", "Mukamel_2021.loci_to_include_in_catalog.filtered.json.gz"),
    ("Annear2021", "Annear_2021.loci_to_include_in_catalog.filtered.json.gz"),
    ("Gymrek2016", "Gymrek_2016_eSTRs.tandem_repeats.filtered.json.gz"),
    ("Hause2016", "Hause_2016_cancer_MSI_loci.filtered.json.gz"),
    ("VamosV3", "vamosGenomicTR_v3.0.loci_to_include_in_catalog.filtered.json.gz"),
]

MERGED_CATALOGS = [
    ("merged_v1", "TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.merged_v1.annotated.json.gz"),
    ("merged_v2a", "TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.merged_v2a.annotated.json.gz"),
]


def parse_loci(path):
    """Parse a JSON catalog and return a list of (chrom, start, end, motif_size) tuples."""
    loci = []
    with gzip.open(path, "rt") as f:
        for record in json.load(f):
            region = record["ReferenceRegion"]
            if isinstance(region, list):
                region = region[0]
            chrom, coords = region.split(":")
            start, end = coords.split("-")
            motif = record.get("Motif") or record.get("CanonicalMotif", "")
            if not motif:
                ls = record.get("LocusStructure", "")
                if ls.startswith("(") and ")*" in ls:
                    motif = ls[1:ls.index(")*")]
            loci.append((chrom, int(start), int(end), len(motif)))
    return loci


def build_ncls(loci):
    """Build per-chromosome NCLS trees from a list of (chrom, start, end, motif_size) tuples.

    Returns dict mapping chrom to (NCLS, list_of_loci_on_that_chrom).
    """
    by_chrom = collections.defaultdict(list)
    for chrom, start, end, motif_size in loci:
        by_chrom[chrom].append((start, end, motif_size))

    trees = {}
    for chrom, intervals in by_chrom.items():
        starts = np.array([s for s, _, _ in intervals], dtype=np.int64)
        ends = np.array([max(e, s + 1) for s, e, _ in intervals], dtype=np.int64)
        ids = np.arange(len(intervals), dtype=np.int64)
        trees[chrom] = (NCLS(starts, ends, ids), intervals)
    return trees


def count_self_overlaps(loci):
    """Count loci that overlap at least one other locus in the same set by >= 2x the larger motif."""
    trees = build_ncls(loci)
    overlapping = set()

    by_chrom = collections.defaultdict(list)
    for i, (chrom, start, end, motif_size) in enumerate(loci):
        by_chrom[chrom].append((i, start, end, motif_size))

    for chrom, chrom_loci in by_chrom.items():
        if chrom not in trees:
            continue
        ncls_tree, intervals = trees[chrom]
        for global_idx, start, end, motif_size in chrom_loci:
            for _, _, local_idx in ncls_tree.find_overlap(start, max(end, start + 1)):
                other_start, other_end, other_motif_size = intervals[local_idx]
                if other_start == start and other_end == end:
                    continue
                overlap_size = min(end, other_end) - max(start, other_start)
                if overlap_size >= 2 * max(motif_size, other_motif_size):
                    overlapping.add(global_idx)
                    break

    return len(overlapping)


def count_overlaps_with_reference(source_loci, ref_trees):
    """Count source loci that overlap at least one reference locus by >= 2x the larger motif."""
    overlapping = 0
    for chrom, start, end, motif_size in source_loci:
        if chrom not in ref_trees:
            continue
        ncls_tree, intervals = ref_trees[chrom]
        for _, _, local_idx in ncls_tree.find_overlap(start, max(end, start + 1)):
            other_start, other_end, other_motif_size = intervals[local_idx]
            overlap_size = min(end, other_end) - max(start, other_start)
            if overlap_size >= 2 * max(motif_size, other_motif_size):
                overlapping += 1
                break
    return overlapping


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results-dir", default=".", help="Directory containing the filtered .json.gz catalog files")
    args = parser.parse_args()

    results_dir = args.results_dir

    # Check that files exist
    all_sources = V2A_SOURCES + V2B_SOURCES
    for name, filename in all_sources + MERGED_CATALOGS:
        path = os.path.join(results_dir, filename)
        if not os.path.isfile(path):
            parser.error(f"File not found: {path}")

    # Load and index merged catalogs
    ref_trees = {}
    for ref_name, ref_filename in MERGED_CATALOGS:
        print(f"Loading {ref_name} from {ref_filename}...")
        ref_loci = parse_loci(os.path.join(results_dir, ref_filename))
        ref_trees[ref_name] = build_ncls(ref_loci)
        print(f"  {len(ref_loci):,} loci loaded")

    # Analyze each source
    results = []
    for name, filename in all_sources:
        path = os.path.join(results_dir, filename)
        print(f"\nAnalyzing {name} ({filename})...")
        loci = parse_loci(path)
        total = len(loci)
        print(f"  {total:,} total loci")

        self_overlaps = count_self_overlaps(loci)
        print(f"  {self_overlaps:,} ({100*self_overlaps/total:.1f}%) self-overlapping loci")

        overlaps_v1 = count_overlaps_with_reference(loci, ref_trees["merged_v1"])
        print(f"  {overlaps_v1:,} ({100*overlaps_v1/total:.1f}%) overlap merged_v1")

        overlaps_v2a = count_overlaps_with_reference(loci, ref_trees["merged_v2a"])
        print(f"  {overlaps_v2a:,} ({100*overlaps_v2a/total:.1f}%) overlap merged_v2a")

        results.append({
            "source": name,
            "total": total,
            "self_overlaps": self_overlaps,
            "overlaps_v1": overlaps_v1,
            "overlaps_v2a": overlaps_v2a,
        })

    # Print summary table
    print("\n" + "=" * 130)
    print(f"{'Source':<35s} {'Total':>10s}  {'Self-overlap':>14s}  {'Overlap v1':>14s}  {'Overlap v2a':>14s}  {'New vs v1':>14s}  {'New vs v2a':>14s}")
    print("-" * 130)
    for r in results:
        t = r["total"]
        if t == 0:
            continue
        so = r["self_overlaps"]
        ov1 = r["overlaps_v1"]
        ov2a = r["overlaps_v2a"]
        new_v1 = t - ov1
        new_v2a = t - ov2a
        print(f"{r['source']:<35s} {t:>10,d}  {so:>8,d} ({100*so/t:4.1f}%)  {ov1:>8,d} ({100*ov1/t:4.1f}%)  {ov2a:>8,d} ({100*ov2a/t:4.1f}%)  {new_v1:>8,d} ({100*new_v1/t:4.1f}%)  {new_v2a:>8,d} ({100*new_v2a/t:4.1f}%)")
    print("=" * 130)


if __name__ == "__main__":
    main()
