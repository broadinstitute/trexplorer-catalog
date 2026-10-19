"""Compute per-source-catalog stats for the loci that each source contributed to the v2 catalog.

Streams the annotated catalog TSV, groups loci by the "Source" column (which records the source catalog
that the locus definition came from), and for each source records:
  - the number of contributed loci
  - the ReferenceRepeatPurity of each contributed locus
  - whether the reference region starts with the annotated motif, ends with it, both, or neither

Writes source_catalog_stats.tsv (one row per source) and purity_values.npz (one array per source).

Usage:
    python3 compute_source_catalog_stats.py --catalog-tsv <path> --reference-fasta <path>
"""

import argparse
import collections
import gzip
import os
import subprocess

import numpy as np
import pysam

# source catalogs in the same order as source_catalogs_in_order in run_all_steps_to_generate_the_hg38_catalog.py,
# paired with the basename of the URL each one is downloaded from. The merge step writes stats keyed by the
# filtered version of that basename, so this is what lets us recover the "loci in source catalog" counts.
SOURCE_CATALOGS = [
    ("v1",  "TRExplorerV1:KnownDiseaseAssociatedLoci",      "variant_catalog_without_offtargets.GRCh38.split"),
    ("v1",  "TRExplorerV1:Illumina174kPolymorphicTRs",      "illumina_variant_catalog.sorted"),
    ("v1",  "TRExplorerV1:PerfectRepeatsInReference",       "hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp"),
    ("v1",  "TRExplorerV1:PolymorphicTRsInT2TAssemblies",   "merged_expansion_hunter_catalog.78_samples"),
    ("v2a", "TRExplorerV2:KnownDiseaseAssociatedLociV2",    "known_disease_associated_loci_v2.loci_to_include_in_catalog"),
    ("v2a", "TRExplorerV2:PolymorphicTRsInT2TAssembliesV2", "combined.321_catalogs.tandem_repeats"),
    ("v2a", "TRExplorerV2:KnownFunctionalVNTRs",            "known_functional_VNTRs.loci_to_include_in_catalog"),
    ("v2b", "TRExplorerV2:HipSTRCatalog",                   "hg38.hipstr_reference.loci_to_include_in_catalog"),
    ("v2b", "TRExplorerV2:AdottoTRsFromDanzi2025",          "Danzi_2025.adotto.loci_to_include_in_catalog"),
    ("v2b", "TRExplorerV2:ClinvarIndelsThatAreTRs2025",     "clinvar_2025_11_03.loci_to_include_in_catalog"),
    ("v2b", "TRExplorerV2:Tanudisastro2025",                "Tanudisastro_2025.loci_to_include_in_catalog"),
    ("v2b", "TRExplorerV2:Manigbas2024",                    "Manigbas_2024_genotyped_TRs"),
    ("v2b", "TRExplorerV2:Sulovari2021",                    "Sulovari_2019_human_specific_STRs"),
    ("v2b", "TRExplorerV2:Garg2021",                        "Garg_2021_eVNTRs_mVNTRs"),
    ("v2b", "TRExplorerV2:Mukamel2021",                     "Mukamel_2021.loci_to_include_in_catalog"),
    ("v2b", "TRExplorerV2:Annear2021",                      "Annear_2021.loci_to_include_in_catalog"),
    ("v2b", "TRExplorerV2:Gymrek2016",                      "Gymrek_2016_eSTRs.tandem_repeats"),
    ("v2b", "TRExplorerV2:Hause2016",                       "Hause_2016_cancer_MSI_loci"),
    ("v2b", "TRExplorerV2:VamosV3",                         "vamosGenomicTR_v3.0.loci_to_include_in_catalog"),
]

# the source folder in compare_with_loci_from_other_sources/ that produces each TRExplorerV2 source catalog
SOURCE_FOLDERS = {
    "TRExplorerV2:KnownDiseaseAssociatedLociV2":    "loci_from_known_disease_associated_loci_v2",
    "TRExplorerV2:PolymorphicTRsInT2TAssembliesV2": "",
    "TRExplorerV2:KnownFunctionalVNTRs":            "loci_from_functional_VNTRs",
    "TRExplorerV2:HipSTRCatalog":                   "loci_from_HipSTR_catalog",
    "TRExplorerV2:AdottoTRsFromDanzi2025":          "loci_from_Danzi_2025",
    "TRExplorerV2:ClinvarIndelsThatAreTRs2025":     "loci_from_clinvar_2025_11_03",
    "TRExplorerV2:Tanudisastro2025":                "loci_from_Tanudisastro_2024__sc-eTRs",
    "TRExplorerV2:Manigbas2024":                    "loci_from_Manigbas_2024",
    "TRExplorerV2:Sulovari2021":                    "loci_from_Sulovari_2019__human-specific_STRs",
    "TRExplorerV2:Garg2021":                        "loci_from_Garg_2021__eVNTRs_mVNTRs",
    "TRExplorerV2:Mukamel2021":                     "loci_from_Mukamel_2021",
    "TRExplorerV2:Annear2021":                      "loci_from_Annear_2021",
    "TRExplorerV2:Gymrek2016":                      "loci_from_Gymrek_2016__eSTRs",
    "TRExplorerV2:Hause2016":                       "loci_from_Hause_2016__cancer",
    "TRExplorerV2:VamosV3":                         "loci_from_Vamos_v3",
}

MOTIF_ANCHOR_CATEGORIES = ["both", "starts_only", "ends_only", "neither"]

# loci with a motif of at most this many base pairs are counted as STRs, and longer ones as VNTRs. This is
# the same threshold as MAX_STR_MOTIF_SIZE in compare_STRchive_with_TRExplorer_catalog.py.
MAX_STR_MOTIF_SIZE = 6

REPEAT_CLASSES = ["all", "STR", "VNTR"]


def get_column_indices(catalog_tsv_path, column_names):
    """Return the 1-based indices of the given columns in the header line of a gzipped TSV.

    Args:
        catalog_tsv_path (str): path of the gzipped TSV file.
        column_names (list): column names to look up.

    Returns:
        list: 1-based column indices, in the same order as column_names.
    """
    with gzip.open(catalog_tsv_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    return [header.index(c) + 1 for c in column_names]


def parse_merge_stats(results_dir, output_prefix):
    """Parse the merge_stats.tsv files to get the number of loci in each filtered source catalog.

    Args:
        results_dir (str): directory containing the merge stats files.
        output_prefix (str): catalog file prefix, eg. "TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs".

    Returns:
        tuple: (dict of source name -> loci in the filtered source catalog,
                dict of source name -> loci added by the merge)
    """
    loci_in_source = {}
    loci_added = {}
    basename_to_source = {f"{basename}.filtered.json.gz": name for _, name, basename in SOURCE_CATALOGS}
    for tier in ["merged_v1", "merged_v2a", "merged_v2b"]:
        path = os.path.join(results_dir, f"{output_prefix}.{tier}.merge_stats.tsv")
        with open(path) as f:
            next(f)
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if fields[0] not in basename_to_source:
                    continue  # the row for the catalog that the tier was merged into
                loci_added[basename_to_source[fields[0]]] = int(fields[1])
                loci_in_source[basename_to_source[fields[0]]] = int(fields[2])
    return loci_in_source, loci_added


def ends_with_motif(reference_sequence, motif):
    """Check whether a reference sequence ends with a complete motif copy followed by a partial copy.

    A locus whose region spans a whole number of repeats ends with the motif itself, but a region that stops
    partway through a repeat ends with a prefix of the motif instead, eg. CAGCAGCAGCA for motif CAG. Both
    count here. Requiring a complete copy before the partial one keeps a single coincidentally matching base
    from qualifying.

    Args:
        reference_sequence (str): reference sequence spanned by the locus, at least as long as the motif.
        motif (str): the motif annotated for the locus.

    Returns:
        bool: True if the sequence ends with motif + motif[:j] for some j from 0 to len(motif) - 1.
    """
    tail = reference_sequence[-(2*len(motif) - 1):].upper()
    for length in range(len(motif), len(tail) + 1):
        if tail[-length:] == (motif + motif)[:length]:
            return True
    return False


def classify_motif_anchoring(reference_sequence, motif):
    """Classify whether a reference region starts with the motif, ends with it, both, or neither.

    Comparison is case-insensitive so that soft-masked reference bases are handled the same as unmasked ones.
    A partial copy of the motif at the end of the region counts as ending with the motif, see ends_with_motif.

    Args:
        reference_sequence (str): reference sequence spanned by the locus.
        motif (str): the motif annotated for the locus.

    Returns:
        str: one of "both", "starts_only", "ends_only", "neither".
    """
    if len(reference_sequence) < len(motif):
        return "neither"
    starts = reference_sequence[:len(motif)].upper() == motif
    ends = ends_with_motif(reference_sequence, motif)
    if starts and ends:
        return "both"
    if starts:
        return "starts_only"
    if ends:
        return "ends_only"
    return "neither"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--catalog-tsv", required=True, help="Annotated catalog TSV (.tsv.gz) with a Source column")
    p.add_argument("--reference-fasta", required=True, help="hg38 reference genome FASTA file")
    p.add_argument("--output-dir", default=os.path.dirname(os.path.abspath(__file__)))
    p.add_argument("-n", "--num-loci", type=int, help="Only process this many loci (for testing)")
    args = p.parse_args()

    results_dir = os.path.dirname(os.path.abspath(args.catalog_tsv))
    output_prefix = os.path.basename(args.catalog_tsv).replace(".EH.with_annotations.tsv.gz", "")
    loci_in_source, loci_added = parse_merge_stats(results_dir, output_prefix)

    columns = ["ReferenceRegion", "Source", "ReferenceRepeatPurity", "Motif"]
    print(f"Reading {', '.join(columns)} from {args.catalog_tsv}")

    fasta = pysam.FastaFile(args.reference_fasta)
    purities = collections.defaultdict(list)
    anchor_counts = collections.defaultdict(collections.Counter)
    skipped_counter = collections.Counter()
    total = 0

    # gzcat | cut is much faster than parsing all 100 columns of the catalog TSV in python
    reader = subprocess.Popen(
        f"gzcat {args.catalog_tsv} | cut -f {','.join(map(str, get_column_indices(args.catalog_tsv, columns)))}",
        shell=True, stdout=subprocess.PIPE, text=True)
    next(reader.stdout)
    for line in reader.stdout:
        reference_region, source, purity, motif = line.rstrip("\n").split("\t")
        total += 1
        if args.num_loci and total > args.num_loci:
            break
        if total % 1_000_000 == 0:
            print(f"Processed {total:,d} loci")
        if not motif:
            skipped_counter["missing motif"] += 1
            continue

        repeat_class = "STR" if len(motif) <= MAX_STR_MOTIF_SIZE else "VNTR"
        purities[source, repeat_class].append(float(purity) if purity else np.nan)

        chrom, start_and_end = reference_region.split(":")
        start_0based, end_1based = start_and_end.split("-")
        anchor_counts[source, repeat_class][classify_motif_anchoring(
            fasta.fetch(chrom, int(start_0based), int(end_1based)), motif)] += 1

    reader.stdout.close()
    reader.wait()
    print(f"Processed {total:,d} loci total")
    for reason, count in skipped_counter.items():
        print(f"Skipped {count:,d} loci with {reason}")

    unexpected_sources = {source for source, _ in purities} - {name for _, name, _ in SOURCE_CATALOGS}
    if unexpected_sources:
        raise ValueError(f"Found unexpected values in the Source column: {sorted(unexpected_sources)}")

    # one array per source catalog and repeat class, with the "all" class holding the STRs and the VNTRs together
    purity_arrays = {}
    for _, source, _ in SOURCE_CATALOGS:
        for repeat_class in REPEAT_CLASSES:
            purity_arrays[source, repeat_class] = np.array(
                purities[source, "STR"] + purities[source, "VNTR"] if repeat_class == "all"
                else purities[source, repeat_class], dtype=np.float32)

    output_tsv_path = os.path.join(args.output_dir, "source_catalog_stats.tsv")
    with open(output_tsv_path, "wt") as f:
        f.write("\t".join([
            "merge_tier", "source_catalog", "source_folder", "repeat_class", "loci_in_source_catalog",
            "loci_added_by_merge", "loci_with_source_annotation", "purity_mean", "purity_median",
            "purity_25th_percentile", "purity_75th_percentile", "purity_fraction_equal_to_1"]
            + [f"count_{c}" for c in MOTIF_ANCHOR_CATEGORIES]) + "\n")
        for tier, source, _ in SOURCE_CATALOGS:
            for repeat_class in REPEAT_CLASSES:
                purity_values = purity_arrays[source, repeat_class]
                counts = [anchor_counts[source, "STR"][c] + anchor_counts[source, "VNTR"][c]
                          if repeat_class == "all" else anchor_counts[source, repeat_class][c]
                          for c in MOTIF_ANCHOR_CATEGORIES]
                f.write("\t".join(map(str, [
                    tier, source, SOURCE_FOLDERS.get(source, ""), repeat_class,
                    loci_in_source[source] if repeat_class == "all" else "",
                    loci_added[source] if repeat_class == "all" else "",
                    len(purity_values),
                    f"{np.nanmean(purity_values):.6f}" if len(purity_values) else "",
                    f"{np.nanmedian(purity_values):.6f}" if len(purity_values) else "",
                    f"{np.nanpercentile(purity_values, 25):.6f}" if len(purity_values) else "",
                    f"{np.nanpercentile(purity_values, 75):.6f}" if len(purity_values) else "",
                    f"{np.mean(purity_values == 1):.6f}" if len(purity_values) else "",
                ] + counts)) + "\n")
    print(f"Wrote stats for {len(SOURCE_CATALOGS)} source catalogs to {output_tsv_path}")

    output_npz_path = os.path.join(args.output_dir, "purity_values.npz")
    np.savez_compressed(output_npz_path, **{f"{source}|{repeat_class}": values
                                            for (source, repeat_class), values in purity_arrays.items()})
    print(f"Wrote purity values to {output_npz_path}")


if __name__ == "__main__":
    main()
