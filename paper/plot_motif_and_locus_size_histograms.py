import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import pandas as pd
from pprint import pprint
import seaborn as sns

mpl.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

GENE_REGIONS = ["5utr", "3utr", "cds", "exon", "intergenic", "intron", "promoter"]


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-t", "--stats-table-path", default="combined_catalog_stats.all_15_catalogs.tsv")
    parser.add_argument("-W", "--grid-width", type=int, default=5, help="Number of columns in the trellis grid")
    parser.add_argument("-H", "--grid-height", type=int, default=3, help="Number of rows in the trellis grid")
    g1 = parser.add_mutually_exclusive_group()
    g1.add_argument("--skip-motif-size-plots", action="store_true", help="Skip plotting motif size histograms")
    g1.add_argument("--only-plot-motif-size", action="store_true", help="Only plot motif size histograms")
    g2 = parser.add_mutually_exclusive_group()
    g2.add_argument("--skip-locus-size-plots", action="store_true", help="Skip plotting locus size histograms")
    g2.add_argument("--only-plot-locus-size", action="store_true", help="Only plot locus size histograms")
    g3 = parser.add_mutually_exclusive_group()
    g3.add_argument("--skip-gene-region-plots", action="store_true", help="Skip plotting gene region histograms")
    g3.add_argument("--only-plot-gene-regions", action="store_true", help="Only plot gene region histograms")
    args = parser.parse_args()

    catalog_stats_table_path = args.stats_table_path
    if not os.path.isfile(catalog_stats_table_path):
        parser.error(f"File not found: {catalog_stats_table_path}")

    #if args.grid_width * args.grid_height != 15:
    #    parser.error("grid_width * grid_height must equal 15")

    df = pd.read_table(catalog_stats_table_path)
    df = df[df["catalog"] != "known_disease_associated_loci"]
    df["catalog"] = df["catalog"].replace({
        "KnownDiseaseAssociatedLoci": "Known disease-associated loci",
        "Illumina174kPolymorphicTRs": "Illumina 174k TRs polymorphic in 1kGP",
        "UCSC_SimpleRepeatTrack": "UCSC simple repeat track",
        "VamosCatalog_v2.1": "vamos catalog v2.1",
        "GangSTR_v17": "GangSTR catalog v17",
        "HipSTR_Catalog": "HipSTR catalog",
        "PolymorphicTRsInT2TAssemblies": "Polymorphic TRs in 78 HPRC assemblies",
        "Adotto_v1.2": "Adotto v1.2",
        "PerfectRepeatsInReference": "All perfect repeats (≥ 3x and ≥ 9bp) in hg38",
        "PopSTR_Catalog": "PopSTR catalog",
        "PlatinumTRs_v1.0": "Platinum TRs v1.0",
        "Chiu_et_al": "Comprehensive catalog from Chiu et al. 2024",
        "TRExplorer_Catalog_v1": "TRExplorer Catalog v1",
    })

    print(f"Parsed {len(df)} catalogs from {catalog_stats_table_path}")

    #pprint(list(df.columns))
    histogram_types = []
    if not args.skip_motif_size_plots and not args.only_plot_locus_size and not args.only_plot_gene_regions:
        histogram_types.append("motif_sizes")
    if not args.skip_locus_size_plots and not args.only_plot_motif_size and not args.only_plot_gene_regions:
        histogram_types.append("locus_sizes")
    if not args.skip_gene_region_plots and not args.only_plot_motif_size and not args.only_plot_locus_size:
        histogram_types.append("gene_regions")
        histogram_types.append("3bp_motif_gene_regions")

    output_paths = []
    for histogram_type in histogram_types:
        output_paths.extend(
            create_plot(histogram_type, df, args))

    print("Done generating", ", ".join(output_paths))


def create_plot(histogram_type, df, args):
    df = df.copy()
    grid_rows = args.grid_height
    grid_columns = args.grid_width

    output_path_prefix = f"variant_catalog_stats.{len(df)}_catalogs.{grid_rows}x{grid_columns}.{histogram_type}"
    print(f"Generating plots for {output_path_prefix}")
    columns_of_interest = []
    column_name_map = {}
    if histogram_type == "motif_sizes":
        for label in list(range(1, 7)) + ["7-24", "25+"]:
            new_column_name = f"{label}bp motifs"
            column_name_map[f"count_{label}bp_motifs"] = new_column_name
            columns_of_interest.append(new_column_name)
    elif histogram_type == "locus_sizes":
        for label in list(range(1, 25)) + ["25-50", "51+"]:
            new_column_name = f"{label}x"
            column_name_map[f"num_repeats_per_locus:{label}x"] = new_column_name
            columns_of_interest.append(new_column_name)
    elif histogram_type == "gene_regions":
        for label in GENE_REGIONS:
            new_column_name = label
            column_name_map[f"count_gene_region_{label}"] = new_column_name
            columns_of_interest.append(new_column_name)
    elif histogram_type == "3bp_motif_gene_regions":
        for label in GENE_REGIONS:
            new_column_name = label
            column_name_map[f"count_3bp_motif_gene_region_{label}"] = new_column_name
            columns_of_interest.append(new_column_name)
    else:
        raise ValueError(f"Invalid histogram_type: {histogram_type}")
    
    df.rename(columns=column_name_map, inplace=True)


    # check if all columns are present
    missing_columns = []
    for column in columns_of_interest:
        if column not in df.columns:
            missing_columns.append(column)

    if missing_columns:
        print(f"WARNING: The following columns are missing from the dataframe:", ", ".join(missing_columns))
        return []


    df.index = range(len(df))

    fig, all_axes = plt.subplots(grid_rows, grid_columns, figsize=(7 * grid_columns, 7 * grid_rows), sharey=True,
                                 sharex=False, constrained_layout=True)  # , dpi=600)

    # check if it's a sequence type or a single axis
    if isinstance(all_axes, plt.Axes):
        all_axes = [[all_axes]]

    i = 0
    for idx_i in range(grid_rows):
        for idx_j in range(grid_columns):
            axes = all_axes[idx_i][idx_j]

            # if idx_i >= 1 or idx_j >= 1:
            #    # don't plot axes
            #    axes.axis("off")
            #    continue

            current_catalog_name = df.loc[i, 'catalog'].replace("_", " ")
            #total = int(df.loc[i, "total"])
            total = sum(df.loc[i, columns_of_interest])

            row_values_i = df.loc[i, columns_of_interest]
            row_values_seaborn_i = row_values_i.reset_index()
            row_values_seaborn_i.columns = ["Category", "Count"]
            row_values_seaborn_i["Category"] = [
                v.replace(" motifs", "") for v in row_values_seaborn_i["Category"]
            ]
            row_values_seaborn_i["Value"] = 100 * row_values_seaborn_i["Count"] / total

            sns.barplot(x="Category", y="Value", data=row_values_seaborn_i, ax=axes, color="cornflowerblue", errorbar=None)

            axes.set_title(f"{current_catalog_name}\n\n{total:,d} loci", fontsize=20, pad=20)
            axes.set_ylabel("% of loci", fontsize=18)
            axes.set_xlabel(" ", fontsize=18)
            axes.grid(axis="y", linestyle="-", linewidth=0.5, color="lightgray")

            # add second y axis with repeat counts
            # ax2 = axes.twinx()
            # ax2.set_ylabel("Total loci", fontsize=14)
            # ax2.set_ylim(0, 100)

            # remote spines
            axes.spines["top"].set_visible(False)
            axes.spines["right"].set_visible(False)
            #if idx_j > 0:
            #    axes.spines["left"].set_visible(False)

            # add y padding to the plot title
            axes.title.set_position([.5, 1.05])

            axes.tick_params(axis="x", labelsize=14, rotation=0, length=0)  # No x-axis tick marks
            axes.tick_params(axis="y", labelsize=14, length=0)

            if histogram_type in ["motif_sizes", "3bp_motif_gene_regions", "gene_regions"]:
                y_max = 100
                counts_and_patches = zip(row_values_seaborn_i["Count"], axes.patches)
            elif histogram_type == "locus_sizes":
                y_max = 100
                largest_5_values = set(list(sorted(row_values_seaborn_i["Count"]))[-5:])
                counts_and_patches = [
                    (locus_count, p) for locus_count, p in zip(row_values_seaborn_i["Count"], axes.patches)
                    if locus_count == 0 or locus_count in largest_5_values
                ]
            else:
                raise ValueError(f"Invalid histogram_type: {histogram_type}")

            axes.set_ylim(0, y_max)
            axes.set_yticks(range(0, y_max + 1, 10))

            for locus_count, p in counts_and_patches:
                x_pos = p.get_x() + 0.5
                y_pos = p.get_height()
                axes.annotate(f"{int(locus_count):,d}", (x_pos, y_pos),
                              ha="center", va="bottom", fontsize=15, color="gray", xytext=(0, 10),
                              textcoords="offset points", rotation=90)

            # Set fixed ticks before setting labels
            axes.set_xticks(range(len(axes.get_xticklabels())))
            axes.set_xticklabels(axes.get_xticklabels(), rotation=90)
            axes.set_xlabel(axes.get_xlabel(), fontsize=18)

            i += 1
            if i >= len(df):
                break

    # Move tight_layout before saving
    plt.tight_layout()

    output_paths = []
    for suffix in "png", "svg":
        output_path = f"{output_path_prefix}.{suffix}"
        print(f"Saving plot to {output_path}")
        plt.savefig(output_path)  # , dpi=300)
        output_paths.append(output_path)

    return output_paths


if __name__ == "__main__":
    main()
