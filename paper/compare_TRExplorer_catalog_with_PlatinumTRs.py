"""Compares a TRExplorer v2.1 catalog with the Platinum TRs v1.0 catalog of the same reference genome (hg38 or CHM13).

The plots, and the steps that compute their inputs, are the ones the paper uses to compare hg38 catalogs:
- the Figure S1 pie charts, made from str_analysis.merge_loci outer join tables, with the loci in each category counted
  the same way as in compute_pairwise_comparisons_for_pie_charts.py (also split by motif size, like its --by-motif-size
  option)
- the Figure 1 motif size and Figure S2 reference repeat count bar plots, made by plot_motif_and_locus_size_histograms.py
  from str_analysis.compute_catalog_stats tables, as in compare_catalogs.py

As in the paper, both catalogs are first trimmed (annotate_and_filter_str_catalog --trim-loci) so that each locus's
width is an exact multiple of its motif size. Everything runs locally. Like make, a rerun skips each step whose output
is newer than its inputs, so after a catalog is regenerated the steps that depend on it run again.

The gene region bar plots are skipped: the CHM13 Gencode annotation (the T2T consortium's CAT/Liftoff GTF) has no UTR
features, so the 5' UTR and 3' UTR bars would always be empty, and the hg38 comparison skips them too so that the two
comparisons have the same plots.
"""

import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import pandas as pd
import re
import subprocess

mpl.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

PLATINUM_TRS_URLS = {
    "hg38": "https://zenodo.org/records/13178746/files/human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.bed.gz",
    "chm13": "https://zenodo.org/records/13178746/files/chm13v2.0_maskedY_rCRS.platinumTRs-v1.0.trgt.bed.gz",
}

# the reference genome names used in catalog labels and plot titles
GENOME_VERSION_NAMES = {
    "hg38": "hg38",
    "chm13": "CHM13",
}

# the motif size bins of compute_pairwise_comparisons_for_pie_charts.py --by-motif-size
MOTIF_SIZE_RANGES = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (7, 24), (25, None)]

# pie chart categories in clockwise order, with the legend labels and colors of Figure S1
PIE_CHART_CATEGORIES = [
    ("unique1", "Locus is only in catalog1", "#a3acb9"),
    ("widerInCatalog1", "Locus has wider boundaries in catalog1", "#d5d5d5"),
    ("intersection", "Locus is in both catalogs", "#8074a8"),
    ("widerInCatalog2", "Locus has wider boundaries in catalog2", "#b9ddf1"),
    ("unique2", "Locus is only in catalog2", "#80b0d5"),
]


def run(command):
    print("-"*100)
    command = re.sub("[ \\t]{2,}", "  ", command)  # remove extra spaces
    print(command)
    subprocess.run(command, shell=True, check=True)


def run_in_parallel(commands):
    """Runs shell commands at the same time and raises an error if any of them fails."""
    processes = []
    for command in commands:
        print("-"*100)
        command = re.sub("[ \\t]{2,}", "  ", command)  # remove extra spaces
        print(command)
        processes.append((command, subprocess.Popen(command, shell=True)))

    failed_commands = [command for command, process in processes if process.wait() != 0]
    if failed_commands:
        raise RuntimeError("Commands failed:\n" + "\n".join(failed_commands))


def is_up_to_date(output_path, *input_paths):
    """Returns True if output_path exists and is at least as new as each of the input paths, like make."""
    return os.path.isfile(output_path) and all(
        os.path.getmtime(output_path) >= os.path.getmtime(input_path) for input_path in input_paths)


def get_motif_size_label(min_motif_size, max_motif_size):
    if min_motif_size == max_motif_size:
        return f"{min_motif_size}bp"
    if max_motif_size is None:
        return f"{min_motif_size}+bp"
    return f"{min_motif_size}-{max_motif_size}bp"


def count_loci_in_each_pie_chart_category(df, catalog1, catalog2):
    """Counts the loci in each pie chart category the same way as compute_pairwise_comparisons_for_pie_charts.py.

    Args:
        df (pandas.DataFrame): rows of a str_analysis.merge_loci outer join table, which has a column per catalog
        catalog1 (str): column name of the first catalog
        catalog2 (str): column name of the second catalog

    Returns:
        dict: maps each pie chart category to its number of loci
    """
    return {
        "intersection": sum(df[catalog1].isin({"Yes", "YesButShifted"}) & df[catalog2].isin({"Yes", "YesButShifted"})),
        "unique1": sum(df[catalog1].notna() & df[catalog2].isna()),
        "unique2": sum(df[catalog1].isna() & df[catalog2].notna()),
        "widerInCatalog1": sum((df[catalog1] == "YesButWider") | (df[catalog2] == "YesButNarrower")),
        "widerInCatalog2": sum((df[catalog2] == "YesButWider") | (df[catalog1] == "YesButNarrower")),
    }


def compute_pie_chart_counts(outer_join_table_path, catalog1, catalog2, catalog_names):
    """Returns the pie chart counts of one catalog pair for all loci and for each motif size bin.

    Args:
        outer_join_table_path (str): path of the str_analysis.merge_loci outer join table of catalog1 vs. catalog2
        catalog1 (str): label of the first catalog
        catalog2 (str): label of the second catalog
        catalog_names (dict): maps each catalog label to the name shown in the plots

    Returns:
        list: one dict per category and motif size bin, with the motif_size, motif_label, category, catalog1, catalog2
            and count columns of the paper's pairwise_catalog_comparison_results.tsv (but not its catalog size
            columns, and with motif_label values like "all", "1bp" and "25+bp" rather than "None", "1" and "25+")
    """
    df_full = pd.read_table(outer_join_table_path)

    # the same checks as in compute_pairwise_comparisons_for_pie_charts.py
    assert all(df_full[catalog1].isin({"Yes", "YesButShifted", "YesButNarrower", "YesButWider", float("nan")}))
    assert sum(df_full[catalog1].isna() & df_full[catalog2].isna()) == 0
    assert sum(df_full[catalog1].isna() & (df_full[catalog2] != "Yes")) == 0
    assert sum(df_full[catalog2].isna() & (df_full[catalog1] != "Yes")) == 0

    # the outer join table's LocusIds are chrom-start-end-motif
    df_full["motif_size"] = df_full["LocusId"].apply(lambda x: len(x.split("-")[3]))

    rows = []
    for min_motif_size, max_motif_size in [(None, None)] + MOTIF_SIZE_RANGES:
        if min_motif_size is None:
            df = df_full
            motif_size_label = "all"
        else:
            filter_condition = df_full["motif_size"] >= min_motif_size
            if max_motif_size is not None:
                filter_condition &= df_full["motif_size"] <= max_motif_size
            df = df_full[filter_condition]
            motif_size_label = get_motif_size_label(min_motif_size, max_motif_size)

        for category, count in count_loci_in_each_pie_chart_category(df, catalog1, catalog2).items():
            rows.append({
                "motif_size": min_motif_size,
                "motif_label": motif_size_label,
                "category": category,
                "catalog1": catalog_names[catalog1],
                "catalog2": catalog_names[catalog2],
                "count": count,
            })

    return rows


def plot_pie_charts(pie_chart_counts_df, output_path_prefix):
    """Plots one pie chart for all loci, followed by one pie chart per motif size bin.

    Args:
        pie_chart_counts_df (pandas.DataFrame): table returned by compute_pie_chart_counts
        output_path_prefix (str): output path without the .png or .svg suffix

    Returns:
        list: paths of the output plots
    """
    catalog1 = pie_chart_counts_df["catalog1"].iloc[0]
    catalog2 = pie_chart_counts_df["catalog2"].iloc[0]
    motif_size_labels = list(dict.fromkeys(pie_chart_counts_df["motif_label"]))

    fig, all_axes = plt.subplots(1, len(motif_size_labels), figsize=(2.4 * len(motif_size_labels), 3.6))
    for axes, motif_size_label in zip(all_axes, motif_size_labels):
        counts_by_category = pie_chart_counts_df[
            pie_chart_counts_df["motif_label"] == motif_size_label].set_index("category")["count"]
        counts = [counts_by_category[category] for category, _, _ in PIE_CHART_CATEGORIES]
        axes.pie(counts, colors=[color for _, _, color in PIE_CHART_CATEGORIES], startangle=90, counterclock=False,
                 autopct=lambda percent: f"{percent:.0f}%" if percent >= 5 else "", pctdistance=0.7,
                 textprops={"fontsize": 9}, wedgeprops={"linewidth": 0.5, "edgecolor": "white"})

        title = "All loci" if motif_size_label == "all" else f"{motif_size_label} motifs"
        axes.set_title(f"{title}\n{sum(counts):,d} loci", fontsize=11)

    fig.legend(handles=[
        mpl.patches.Patch(color=color, label=label.replace("catalog1", catalog1).replace("catalog2", catalog2))
        for _, label, color in PIE_CHART_CATEGORIES
    ], loc="upper center", ncol=3, fontsize=11, frameon=False)
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.02, top=0.7, wspace=0.1)

    output_paths = []
    for suffix in "png", "svg":
        output_path = f"{output_path_prefix}.{suffix}"
        fig.savefig(output_path, dpi=150)
        print(f"Saved pie charts to {output_path}")
        output_paths.append(output_path)
    plt.close(fig)

    return output_paths


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--genome-version", choices=sorted(PLATINUM_TRS_URLS), required=True,
                        help="Reference genome of the TRExplorer catalog")
    parser.add_argument("-R", "--reference-fasta", required=True,
                        help="Path of the reference genome FASTA, eg. hg38.fa or chm13v2.0.maskedY.fasta")
    parser.add_argument("--trexplorer-catalog", required=True,
                        help="Path of the TRExplorer v2.1 catalog in ExpansionHunter JSON format, eg. the "
                             "TRExplorer.repeat_catalog_v2.1.chm13.1_to_1000bp_motifs.EH.json.gz release file")
    parser.add_argument("--output-dir", help="Directory for the downloaded catalog, intermediate files and plots. "
                        "Defaults to compare_catalogs_<genome version>")
    args = parser.parse_args()

    for path in args.reference_fasta, args.trexplorer_catalog:
        if not os.path.isfile(path):
            parser.error(f"File not found: {path}")

    reference_fasta = os.path.abspath(args.reference_fasta)
    base_dir = os.path.dirname(os.path.abspath(__file__))
    genome_version_name = GENOME_VERSION_NAMES[args.genome_version]
    trexplorer_label = f"TRExplorer_v2.1_{genome_version_name}"
    platinum_trs_label = f"PlatinumTRs_v1.0_{genome_version_name}"

    # maps the catalog labels used in file names and outer join table columns to the names shown in the plots
    catalog_names = {
        trexplorer_label: f"TRExplorer v2.1 {genome_version_name}",
        platinum_trs_label: f"Platinum TRs v1.0 {genome_version_name}",
    }
    catalog_paths = {trexplorer_label: os.path.abspath(args.trexplorer_catalog)}

    output_dir = args.output_dir or f"compare_catalogs_{args.genome_version}"
    os.makedirs(output_dir, exist_ok=True)
    print(f"cd {output_dir}")
    os.chdir(output_dir)

    # download the Platinum TRs catalog and convert it from TRGT to ExpansionHunter format, as in compare_catalogs.py
    platinum_trs_url = PLATINUM_TRS_URLS[args.genome_version]
    platinum_trs_trgt_bed_path = os.path.abspath(os.path.basename(platinum_trs_url))
    if not os.path.isfile(platinum_trs_trgt_bed_path):
        run(f"wget -O {platinum_trs_trgt_bed_path}.tmp -q {platinum_trs_url} && "
            f"mv {platinum_trs_trgt_bed_path}.tmp {platinum_trs_trgt_bed_path}")

    catalog_paths[platinum_trs_label] = platinum_trs_trgt_bed_path.replace(".bed.gz", ".catalog.json.gz")
    if not is_up_to_date(catalog_paths[platinum_trs_label], platinum_trs_trgt_bed_path):
        run(f"python3 -u -m str_analysis.convert_trgt_catalog_to_expansion_hunter_catalog "
            f"-r {reference_fasta} {platinum_trs_trgt_bed_path} -o {catalog_paths[platinum_trs_label]}")

    # trim both catalogs
    trimmed_catalog_paths = {label: os.path.abspath(f"trimmed.{label}.json.gz") for label in catalog_paths}
    run_in_parallel([
        f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog \
            --verbose \
            --trim-loci \
            --reference-fasta {reference_fasta} \
            --skip-gene-annotations \
            --skip-disease-loci-annotations \
            --skip-mappability-annotations \
            --output-path {trimmed_catalog_paths[label]} \
            {catalog_paths[label]} > trimmed.{label}.log 2>&1"""
        for label in catalog_paths if not is_up_to_date(trimmed_catalog_paths[label], catalog_paths[label])
    ])

    # compute the stats for the bar plots, and combine them into one table as in compare_catalogs.py
    stats_dfs = []
    for label, trimmed_catalog_path in trimmed_catalog_paths.items():
        stats_tsv_path = re.sub("(.json|.bed)(.gz)?$", "", os.path.basename(trimmed_catalog_path)) + ".catalog_stats.tsv"
        if not is_up_to_date(stats_tsv_path, trimmed_catalog_path):
            run(f"python3 -u -m str_analysis.compute_catalog_stats --verbose {trimmed_catalog_path} > {stats_tsv_path}.log 2>&1")
        df = pd.read_table(stats_tsv_path)
        df.rename(columns={"catalog": "filename"}, inplace=True)
        df["catalog"] = catalog_names[label]
        stats_dfs.append(df)

    stats_df = pd.concat(stats_dfs)
    stats_df["total"] = stats_df["total"].apply(lambda s: str(s).replace(",", "")).astype(int)
    stats_df.sort_values("total", inplace=True)
    combined_stats_tsv_path = f"combined_catalog_stats.{len(stats_df)}_catalogs.tsv"
    stats_df.to_csv(combined_stats_tsv_path, sep="\t", index=False)
    print(f"Wrote {combined_stats_tsv_path}")

    run(f"python3 -u {base_dir}/plot_motif_and_locus_size_histograms.py -t {combined_stats_tsv_path} "
        f"-W {len(stats_df)} -H 1 --skip-gene-region-plots")

    # compare the catalogs using the merge_loci options of compute_pairwise_comparisons_for_pie_charts.py. Figure S1
    # also compares each pair in the reverse order, but for the CHM13 catalogs that gave the mirror image: the same
    # outer join table rows, with every category count within 0.03% of the swapped category.
    catalog1, catalog2 = trexplorer_label, platinum_trs_label
    output_prefix = f"merged___{catalog1}___vs___{catalog2}"
    outer_join_table_path = f"{output_prefix}.outer_join_overlap_table.tsv.gz"
    if not is_up_to_date(outer_join_table_path, trimmed_catalog_paths[catalog1], trimmed_catalog_paths[catalog2]):
        run(f"""python3 -u -m str_analysis.merge_loci \
            --discard-extra-fields-from-input-catalogs \
            --overlap-fraction 0.05 \
            --write-outer-join-table \
            --output-format JSON \
            --output-prefix {output_prefix} \
            {catalog1}:{trimmed_catalog_paths[catalog1]} \
            {catalog2}:{trimmed_catalog_paths[catalog2]} > {output_prefix}.log 2>&1""")

    pie_chart_counts_df = pd.DataFrame(compute_pie_chart_counts(outer_join_table_path, catalog1, catalog2, catalog_names))
    pie_chart_counts_tsv_path = f"pairwise_catalog_comparison_results.{args.genome_version}.tsv"
    pie_chart_counts_df.to_csv(pie_chart_counts_tsv_path, sep="\t", index=False)
    print(f"Wrote {pie_chart_counts_tsv_path}")

    plot_pie_charts(pie_chart_counts_df, f"catalog_pie_charts.{args.genome_version}")


if __name__ == "__main__":
    main()
