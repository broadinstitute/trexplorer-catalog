"""This script takes a repeat catalog in BED format and determines how well its loci are represented in
the TRExplorer catalog (or some other main TR catalog like Adotto or Platinum).
"""


import argparse
import collections
import intervaltree
import gzip
import matplotlib.pyplot as plt
import os
import pandas as pd
import pyfaidx
import random
import re
import seaborn as sns
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.find_repeat_unit import find_repeat_unit_without_allowing_interruptions
from str_analysis.utils.find_motif_utils import compute_sequence_purity_stats

MOTIF_MATCH_SCORE_FOR_SAME_MOTIF = 4
MOTIF_MATCH_SCORE_FOR_SAME_MOTIF_LENGTH = 3
MOTIF_MATCH_SCORE_FOR_SHORTER_MOTIF_IN_MAIN_CATALOG = 2
MOTIF_MATCH_SCORE_FOR_LONGER_MOTIF_IN_MAIN_CATALOG = 1
MOTIF_MATCH_SCORE_FOR_ABSENT_FROM_MAIN_CATALOG = 0
MOTIF_MATCH_SCORE_FOR_ABSENT_FROM_NEW_CATALOG = -1
MOTIF_MATCH_SCORE_MAP = {
    MOTIF_MATCH_SCORE_FOR_SAME_MOTIF: "same motif",
    MOTIF_MATCH_SCORE_FOR_SAME_MOTIF_LENGTH: "same motif length",
    MOTIF_MATCH_SCORE_FOR_SHORTER_MOTIF_IN_MAIN_CATALOG: "shorter motif in main catalog",
    MOTIF_MATCH_SCORE_FOR_LONGER_MOTIF_IN_MAIN_CATALOG: "longer motif in main catalog",
    MOTIF_MATCH_SCORE_FOR_ABSENT_FROM_MAIN_CATALOG: "absent from main catalog",
    MOTIF_MATCH_SCORE_FOR_ABSENT_FROM_NEW_CATALOG: "absent from new catalog",
}

OVERLAP_SCORE_FOR_EXACT_MATCH = 7
OVERLAP_SCORE_FOR_DIFF_2_REPEATS = 6
OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_ABOVE_0_66 = 5
OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BETWEEN_0_5_AND_0_66 = 4
OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BETWEEN_0_33_AND_0_5 = 3
OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BETWEEN_0_2_AND_0_33 = 2
OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BELOW_0_2 = 1
OVERLAP_SCORE_FOR_ABSENT_FROM_MAIN_CATALOG = 0
OVERLAP_SCORE_FOR_ABSENT_FROM_NEW_CATALOG = -1
OVERLAP_SCORE_MAP = {
    OVERLAP_SCORE_FOR_EXACT_MATCH: "exact match",
    OVERLAP_SCORE_FOR_DIFF_2_REPEATS: "diff ≤ 2 repeats",
    OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_ABOVE_0_66: "Jaccard similarity > 0.66",
    OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BETWEEN_0_5_AND_0_66: "0.5 < Jaccard similarity ≤ 0.66",
    OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BETWEEN_0_33_AND_0_5: "0.33 < Jaccard similarity ≤ 0.5",
    OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BETWEEN_0_2_AND_0_33: "0.2 < Jaccard similarity ≤ 0.33",
    OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BELOW_0_2: "Jaccard similarity ≤ 0.2",
    OVERLAP_SCORE_FOR_ABSENT_FROM_MAIN_CATALOG: "absent from main catalog",
    OVERLAP_SCORE_FOR_ABSENT_FROM_NEW_CATALOG: "absent from new catalog",
}

NO_MATCH_FOUND = "no (no overlapping definitions)"

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-n", type=int, help="Number of loci to process")
    p.add_argument("-R", "--reference-fasta", help="Reference genome FASTA file", default="~/hg38.fa")
    p.add_argument("--print-stats", type=int, default=1, choices=[0, 1, 2, 3],
                   help="At the end, output some stats at this level of detail. The higher the number, the more detailed the stats")
    p.add_argument("--skip-plots", action="store_true", help="Skip generating plots")
    p.add_argument("--plot-output-dir", default="plots", help="Directory to write plots to")
    p.add_argument("--write-loci-absent-from-new-catalog", action="store_true",
                   help="When generating the output table, include loci that are absent from the new catalog")
    p.add_argument("--write-bed-files-with-subsets", action="store_true",
                   help="Output BED files of loci with different concordance levels to enable manual review")
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
        args.catalog_name = re.sub(".bed(.b?gz)?$", "", os.path.basename(args.catalog_bed_path))

    if not os.path.isfile(args.catalog_bed_path):
        p.error(f"File not found: {args.catalog_bed_path}")

    if not os.path.isdir(args.plot_output_dir):
        print(f"Creating plot output directory: {args.plot_output_dir}")
        os.makedirs(args.plot_output_dir)

    reference_fasta_path = os.path.expanduser(args.reference_fasta)
    if not os.path.isfile(reference_fasta_path):
        p.error(f"File not found: {reference_fasta_path}")
    reference_fasta = pyfaidx.Fasta(reference_fasta_path, as_raw=True, one_based_attributes=False, sequence_always_upper=True)


    new_catalog_names = []
    for i, path in enumerate(args.new_catalog):
        if ":" in path:
            if path.count(":") > 1:
                p.error(f"More than one ':' found in {path}")
            new_catalog_name, path = path.split(":")
            args.new_catalog[i] = path
            new_catalog_names.append(new_catalog_name)
        else:
            new_catalog_names.append(re.sub(".bed(.b?gz)?$", "", os.path.basename(path)))

        if not os.path.isfile(path):
            p.error(f"File not found: {path}")


    OVERLAP_SCORE_MAP[OVERLAP_SCORE_FOR_ABSENT_FROM_MAIN_CATALOG] = f"absent from {args.catalog_name} catalog"
    MOTIF_MATCH_SCORE_MAP[MOTIF_MATCH_SCORE_FOR_ABSENT_FROM_MAIN_CATALOG] = f"absent from {args.catalog_name} catalog"
    MOTIF_MATCH_SCORE_MAP[MOTIF_MATCH_SCORE_FOR_LONGER_MOTIF_IN_MAIN_CATALOG] = f"longer motif in {args.catalog_name} catalog"
    MOTIF_MATCH_SCORE_MAP[MOTIF_MATCH_SCORE_FOR_SHORTER_MOTIF_IN_MAIN_CATALOG] = f"shorter motif in {args.catalog_name} catalog"
    OVERLAP_SCORE_MAP[OVERLAP_SCORE_FOR_ABSENT_FROM_NEW_CATALOG] = f"only in {args.catalog_name} catalog"
    MOTIF_MATCH_SCORE_MAP[MOTIF_MATCH_SCORE_FOR_ABSENT_FROM_NEW_CATALOG] = f"only in {args.catalog_name} catalog"

    main_catalog_loci = load_main_catalog_loci(args)

    for new_catalog_path, new_catalog_name in zip(args.new_catalog, new_catalog_names):
        output_filename_prefix = new_catalog_path.replace(".bed", "").replace(".gz", "")
        output_tsv = f"{output_filename_prefix}.overlap_with_{args.catalog_name}.tsv.gz"
        print(f"Comparing {output_tsv}")
        df = compare_loci(
            main_catalog_loci,
            new_catalog_path,
            main_catalog_name=args.catalog_name,
            new_catalog_name=new_catalog_name,
            reference_fasta=reference_fasta,
            write_loci_absent_from_new_catalog=args.write_loci_absent_from_new_catalog,
        )
        df.sort_values(by=["chrom", "start_0based", "end_1based"], inplace=True)
        df.to_csv(output_tsv, sep="\t", index=False)
        print(f"Wrote {len(df):,d} rows to {output_tsv}")

        if args.write_bed_files_with_subsets:

            # output BED file with loci that have Jaccard similarity < 0.2
            df_filtered = df[df["overlap_score"] == OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BELOW_0_2]
            if len(df_filtered) == 0:
                print(f"No loci were found to have Jaccard similarity < 0.2")
            else:
                output_bed_path = f"loci_from_{new_catalog_name}.jaccard_similarity_below_0_2.bed"
                df_filtered[["chrom", "start_0based", "end_1based", f"{new_catalog_name}_motif"]].to_csv(output_bed_path, sep="\t", index=False, header=False)
                os.system(f"bgzip -f {output_bed_path}")
                os.system(f"tabix -f {output_bed_path}.gz")
                print(f"Wrote {len(df_filtered):,d} rows to {output_bed_path}.gz")

            # output BED file with loci that are absent from the new catalog
            df_filtered = df[df["overlap_score"] == OVERLAP_SCORE_FOR_ABSENT_FROM_NEW_CATALOG]
            if len(df_filtered) == 0:
                print(f"No loci were found to be absent from {new_catalog_name}")
            else:
                output_bed_path = f"loci_from_{args.catalog_name}.absent_from_{new_catalog_name}.bed"
                df_filtered[["chrom", "start_0based", "end_1based", f"{args.catalog_name}_motif"]].to_csv(output_bed_path, sep="\t", index=False, header=False)
                os.system(f"bgzip -f {output_bed_path}")
                os.system(f"tabix -f {output_bed_path}.gz")
                print(f"Wrote {len(df_filtered):,d} rows to {output_bed_path}.gz")

            # output BED file with loci that are absent from the main catalog
            df_filtered = df[df["overlap_score"] == OVERLAP_SCORE_FOR_ABSENT_FROM_MAIN_CATALOG]
            if len(df_filtered) == 0:
                print(f"No loci were found to be absent from {args.catalog_name}")
            else:
                output_bed_path = f"loci_from_{new_catalog_name}.absent_from_{args.catalog_name}.bed"
                df_filtered[["chrom", "start_0based", "end_1based", f"{new_catalog_name}_motif"]].to_csv(output_bed_path, sep="\t", index=False, header=False)
                os.system(f"bgzip -f {output_bed_path}")
                os.system(f"tabix -f {output_bed_path}.gz")
                print(f"Wrote {len(df_filtered):,d} rows to {output_bed_path}.gz")

            # output BED file with Jaccard similarity < 0.2

            # output BED files with loci that have the same boundaries in the two catalogs but different motifs
            df_filtered = df[(df["overlap_score"] == OVERLAP_SCORE_FOR_EXACT_MATCH) & (df["motif_match_score"] <= 1)]
            if len(df_filtered) == 0:
                print("No loci have the same boundaries in the two catalogs but different motifs")
            else:
                print(f"Found {len(df_filtered):,d} out of {len(df):,d} ({len(df_filtered) / len(df):.1%}) loci have the same boundaries in the two catalogs but different motifs")

                output_bed1_path = f"loci_from_{args.catalog_name}.same_boundaries_but_different_motifs.bed"
                df_filtered[["chrom", "start_0based", "end_1based", f"{args.catalog_name}_motif"]].to_csv(output_bed1_path, sep="\t", index=False, header=False)
                os.system(f"bgzip -f {output_bed1_path}")
                os.system(f"tabix -f {output_bed1_path}.gz")
                print(f"Wrote {len(df_filtered):,d} rows to {output_bed1_path}.gz")

                output_bed2_path = f"loci_from_{new_catalog_name}.same_boundaries_but_different_motifs.bed"
                df_filtered[["chrom", "start_0based", "end_1based", f"{new_catalog_name}_motif"]].to_csv(output_bed2_path, sep="\t", index=False, header=False)
                os.system(f"bgzip -f {output_bed2_path}")
                os.system(f"tabix -f {output_bed2_path}.gz")
                print(f"Wrote {len(df_filtered):,d} rows to {output_bed2_path}.gz")

        if args.print_stats > 0:
            print_stats(args, main_catalog_loci, df, new_catalog_name=new_catalog_name)

def print_stats(args, main_catalog_loci, df, new_catalog_name):
    # print some stats
    total_count = collections.Counter()
    main_catalog_motif_size_distribution = collections.Counter()
    main_catalog_reference_repeat_count_distribution = collections.Counter()
    for _, tree in main_catalog_loci.items():
        for interval in tree:
            total_count["TRs"] += 1
            canonical_motif = interval.data["canonical_motif"]  # this motif was simplified before computing the canonical motif
            main_catalog_motif_size_distribution[len(canonical_motif)] += 1
            reference_repeat_count = (interval.end - interval.begin) // len(canonical_motif)
            main_catalog_reference_repeat_count_distribution[reference_repeat_count] += 1
            if len(canonical_motif) <= 6:
                total_count["STRs"] += 1
            else:
                total_count["VNTRs"] += 1


    if args.print_stats >= 1:
        print("Summary:")
        for key, count in sorted(df["match_found?"].apply(lambda x: x.split(" ")[0]).value_counts().items(), reverse=True):
            print(f"{count:10,d}  {key}")

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

    if args.skip_plots:
        return

    new_catalog_has_motifs = sum(df[f"{new_catalog_name}_canonical_motif"].notna()) > 0

    df1 = df[df[f"{new_catalog_name}_canonical_motif"].notna()]
    if len(df1) > 0:
        df1["motif_size"] = df1[f"{new_catalog_name}_canonical_motif"].str.len()
        new_catalog_motif_size_distribution = collections.Counter(df1["motif_size"])
        new_catalog_motif_size_distribution2 = collections.Counter(df1[df1["overlap_score"] == 0]["motif_size"])
        catalog_map = {
            args.catalog_name: main_catalog_motif_size_distribution,
            new_catalog_name: new_catalog_motif_size_distribution,
            f"{new_catalog_name}_loci_absent_from_{args.catalog_name}": new_catalog_motif_size_distribution2,
        }
        for catalog_name, motif_size_distribution in catalog_map.items():
            # bin values between 7 and 24 as "7-24" and 25+ as "25+"
            motif_size_distribution_binned = {f"{motif_size}bp": 0 for motif_size in range(1, 25)}
            motif_size_distribution_binned["25+"] = 0
            for motif_size, count in motif_size_distribution.items():
                if motif_size < 25:
                    motif_size_distribution_binned[f"{motif_size}bp"] += count
                else:
                    motif_size_distribution_binned["25+"] += count
            catalog_map[catalog_name] = motif_size_distribution_binned

        for catalog_name, motif_size_distribution in catalog_map.items():
            catalog_name = catalog_name.replace(" ", "_")
            plt.figure(figsize=(12, 6))
            sns.barplot(x=list(motif_size_distribution.keys()), y=list(motif_size_distribution.values()), color="cornflowerblue")
            # turn the x labels 45 degrees
            plt.xticks(rotation=45)
            plt.xlabel("Motif size")
            plt.ylabel("Count")
            plt.gca().spines["top"].set_visible(False)
            plt.gca().spines["right"].set_visible(False)
            plt.gca().grid(axis="y", linestyle="-", linewidth=0.5, color="lightgray")

            # add , separators to y-axis labels
            plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: "{:,.0f}".format(x)))
            plt.title(f"{catalog_name.replace('_', ' ')} motif size distribution", pad=10)
            output_path = f"{args.plot_output_dir}/{catalog_name}.motif_size_distribution.png"
            plt.savefig(output_path)
            print(f"Wrote motif size distribution to {output_path}")
            plt.close()

        # add a single bar plot that shows multiple bars side by side for each motif size - one for each catalog
        # Use dual y-axes: left axis for first catalog (main), right axis for 2nd and 3rd catalogs
        motif_size_order = [f"{i}bp" for i in range(1, 25)] + ["25+"]
        catalog_names = list(catalog_map.keys())

        fig, ax1 = plt.subplots(figsize=(14, 6))
        ax2 = ax1.twinx()

        # Get data for each catalog
        x_positions = range(len(motif_size_order))
        bar_width = 0.25

        # Plot first catalog (main) on left axis
        first_catalog_name = catalog_names[0]
        first_catalog_counts = [catalog_map[first_catalog_name].get(motif_size, 0) for motif_size in motif_size_order]
        bars1 = ax1.bar([x - bar_width for x in x_positions], first_catalog_counts, width=bar_width,
                        label=first_catalog_name, color="cornflowerblue", alpha=0.8)

        # Plot 2nd and 3rd catalogs on right axis
        if len(catalog_names) >= 2:
            second_catalog_name = catalog_names[1]
            second_catalog_counts = [catalog_map[second_catalog_name].get(motif_size, 0) for motif_size in motif_size_order]
            bars2 = ax2.bar(x_positions, second_catalog_counts, width=bar_width,
                           label=second_catalog_name, color="coral", alpha=0.8)

        if len(catalog_names) >= 3:
            third_catalog_name = catalog_names[2]
            third_catalog_counts = [catalog_map[third_catalog_name].get(motif_size, 0) for motif_size in motif_size_order]
            bars3 = ax2.bar([x + bar_width for x in x_positions], third_catalog_counts, width=bar_width,
                           label=third_catalog_name, color="lightgreen", alpha=0.8)

        # Set x-axis labels
        ax1.set_xticks(x_positions)
        ax1.set_xticklabels(motif_size_order, rotation=45, ha="right")
        ax1.set_xlabel("Motif size")

        # Set y-axis labels and formatting
        ax1.set_ylabel(f"Count ({first_catalog_name})", color="cornflowerblue")
        ax1.tick_params(axis="y", labelcolor="cornflowerblue")
        ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: "{:,.0f}".format(x)))

        if len(catalog_names) >= 2:
            ax2.set_ylabel(f"Count ({catalog_names[1]})\nCount ({catalog_names[2]})", color="coral")
            ax2.tick_params(axis="y", labelcolor="coral")
            ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: "{:,.0f}".format(x)))

        # Styling
        ax1.spines["top"].set_visible(False)
        ax2.spines["top"].set_visible(False)
        ax1.grid(axis="y", linestyle="-", linewidth=0.5, color="lightgray", alpha=0.5)
        ax2.grid(axis="y", linestyle="--", linewidth=0.5, color="lightgray", alpha=0.5)

        # Combine legends
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, frameon=False, loc="best")

        plt.title("Motif size distribution for each catalog", pad=10)
        plt.tight_layout()
        plt.savefig(f"{args.plot_output_dir}/{args.catalog_name}_vs_{new_catalog_name}.motif_size_distribution.png")
        print(f"Wrote motif size distribution to {args.plot_output_dir}/{args.catalog_name}_vs_{new_catalog_name}.motif_size_distribution.png")
        plt.close()

    df2 = df[df[f"{new_catalog_name}_reference_repeat_count"].notna()]
    if len(df2) > 0:
        new_catalog_reference_repeat_count_distribution = collections.Counter(df2[f"{new_catalog_name}_reference_repeat_count"])
        new_catalog_reference_repeat_count_distribution2 = collections.Counter(df2[df2["overlap_score"] == 0][f"{new_catalog_name}_reference_repeat_count"])
        catalog_map = {
            args.catalog_name: main_catalog_reference_repeat_count_distribution,
            new_catalog_name: new_catalog_reference_repeat_count_distribution,
            f"{new_catalog_name}_loci_absent_from_{args.catalog_name}": new_catalog_reference_repeat_count_distribution2,
        }
        for catalog_name, reference_repeat_count_distribution in catalog_map.items():
            # bin values between 1 and 10 as "1-10" and 11+ as "11+"
            reference_repeat_count_distribution_binned = {f"{repeat_count}x": 0 for repeat_count in range(0, 11)}
            reference_repeat_count_distribution_binned["11+"] = 0
            for repeat_count, count in reference_repeat_count_distribution.items():
                if repeat_count < 11:
                    repeat_count = int(repeat_count)
                    reference_repeat_count_distribution_binned[f"{repeat_count}x"] += count
                else:
                    reference_repeat_count_distribution_binned["11+"] += count
            catalog_map[catalog_name] = reference_repeat_count_distribution_binned

        for catalog_name, reference_repeat_count_distribution in catalog_map.items():
            catalog_name = catalog_name.replace(" ", "_")
            plt.figure(figsize=(12, 6))
            sns.barplot(x=list(reference_repeat_count_distribution.keys()), y=list(reference_repeat_count_distribution.values()), color="cornflowerblue")
            # turn the x labels 45 degrees
            plt.xticks(rotation=45)
            plt.xlabel("Reference repeat count")
            plt.ylabel("Count")
            plt.title(f"{catalog_name.replace('_', ' ')} reference repeat count distribution", pad=10)
            plt.gca().spines["top"].set_visible(False)
            plt.gca().spines["right"].set_visible(False)
            plt.gca().grid(axis="y", linestyle="-", linewidth=0.5, color="lightgray")

            # add , separators to y-axis labels
            plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: "{:,.0f}".format(x)))
            output_path = f"{args.plot_output_dir}/{catalog_name}.reference_repeat_count_distribution.png"
            plt.savefig(output_path)
            print(f"Wrote reference repeat count distribution to {output_path}")
            plt.close()


        # add a single bar plot that shows multiple bars side by side for each motif size - one for each catalog
        # Use dual y-axes: left axis for first catalog (main), right axis for 2nd and 3rd catalogs
        repeat_count_order = [f"{i}x" for i in range(0, 11)] + ["11+"]
        catalog_names = list(catalog_map.keys())

        fig, ax1 = plt.subplots(figsize=(14, 6))
        ax2 = ax1.twinx()

        # Get data for each catalog
        x_positions = range(len(repeat_count_order))
        bar_width = 0.25

        # Plot first catalog (main) on left axis
        first_catalog_name = catalog_names[0]
        first_catalog_counts = [catalog_map[first_catalog_name].get(repeat_count, 0) for repeat_count in repeat_count_order]
        bars1 = ax1.bar([x - bar_width for x in x_positions], first_catalog_counts, width=bar_width,
                        label=first_catalog_name, color="cornflowerblue", alpha=0.8)

        # Plot 2nd and 3rd catalogs on right axis
        if len(catalog_names) >= 2:
            second_catalog_name = catalog_names[1]
            second_catalog_counts = [catalog_map[second_catalog_name].get(repeat_count, 0) for repeat_count in repeat_count_order]
            bars2 = ax2.bar(x_positions, second_catalog_counts, width=bar_width,
                           label=second_catalog_name, color="coral", alpha=0.8)

        if len(catalog_names) >= 3:
            third_catalog_name = catalog_names[2]
            third_catalog_counts = [catalog_map[third_catalog_name].get(repeat_count, 0) for repeat_count in repeat_count_order]
            bars3 = ax2.bar([x + bar_width for x in x_positions], third_catalog_counts, width=bar_width,
                           label=third_catalog_name, color="lightgreen", alpha=0.8)

        # Set x-axis labels
        ax1.set_xticks(x_positions)
        ax1.set_xticklabels(repeat_count_order, rotation=45, ha="right")
        ax1.set_xlabel("Reference repeat count")

        # Set y-axis labels and formatting
        ax1.set_ylabel(f"Count ({first_catalog_name})", color="cornflowerblue")
        ax1.tick_params(axis="y", labelcolor="cornflowerblue")
        ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: "{:,.0f}".format(x)))

        if len(catalog_names) >= 2:
            ax2.set_ylabel(f"Count ({catalog_names[1]})\nCount ({catalog_names[2]})", color="coral")
            ax2.tick_params(axis="y", labelcolor="coral")
            ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: "{:,.0f}".format(x)))

        # Styling
        ax1.spines["top"].set_visible(False)
        ax2.spines["top"].set_visible(False)
        ax1.grid(axis="y", linestyle="-", linewidth=0.5, color="lightgray", alpha=0.5)
        ax2.grid(axis="y", linestyle="--", linewidth=0.5, color="lightgray", alpha=0.5)

        # Combine legends
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, frameon=False, loc="best")

        plt.title("Reference repeat count distribution for each catalog", pad=10)
        plt.tight_layout()
        plt.savefig(f"{args.plot_output_dir}/{args.catalog_name}_vs_{new_catalog_name}.reference_repeat_count_distribution.png")
        print(f"Wrote reference repeat count distribution to {args.plot_output_dir}/{args.catalog_name}_vs_{new_catalog_name}.reference_repeat_count_distribution.png")
        plt.close()

    # horizontal bar plot of the values in the 'overlap' column with color by 'motif_match' column
    plt.figure(figsize=(12, 6))
    y_axis_order = list(OVERLAP_SCORE_MAP.values())
    if sum(df["overlap_score"] < 0)/len(df) > 0.85:
        y_axis_order = y_axis_order[:-1]
    elif sum(df["overlap_score"] < 1)/len(df) > 0.85:
        y_axis_order = y_axis_order[:-2]
    elif not args.write_loci_absent_from_new_catalog:
        y_axis_order = y_axis_order[:-1]

    if not new_catalog_has_motifs:
        y_axis_order = y_axis_order[1:]

    # use stacked bar plot to plot the values in the 'overlap' column with color by 'motif_match' column
    plot_data = df.groupby("overlap")["motif_match"].value_counts().reset_index(name="count")
    plot_data_pivot = plot_data.pivot(index="overlap", columns="motif_match", values="count")
    plot_data_pivot = plot_data_pivot.reindex(y_axis_order[::-1]).fillna(0)

    # reorder columns to match legend order
    desired_order = list(MOTIF_MATCH_SCORE_MAP.values())
    plot_data_pivot = plot_data_pivot[[col for col in desired_order[::-1] if col in plot_data_pivot.columns]]
    if new_catalog_has_motifs:
        alpha = 0.5
        motif_match_colors = {
            MOTIF_MATCH_SCORE_MAP[MOTIF_MATCH_SCORE_FOR_SAME_MOTIF]: (0.0, 0.7, 0.0, alpha),                     # green
            MOTIF_MATCH_SCORE_MAP[MOTIF_MATCH_SCORE_FOR_SAME_MOTIF_LENGTH]: (0.95, 0.8, 0.1, alpha),             # yellow
            MOTIF_MATCH_SCORE_MAP[MOTIF_MATCH_SCORE_FOR_SHORTER_MOTIF_IN_MAIN_CATALOG]: (0.95, 0.0, 0.0, alpha),  # light red
            MOTIF_MATCH_SCORE_MAP[MOTIF_MATCH_SCORE_FOR_LONGER_MOTIF_IN_MAIN_CATALOG]: (0.5, 0.0, 0.0, alpha),   # dark red
            MOTIF_MATCH_SCORE_MAP[MOTIF_MATCH_SCORE_FOR_ABSENT_FROM_MAIN_CATALOG]: (0.0, 0.0, 0.7, alpha),       #  blue
            MOTIF_MATCH_SCORE_MAP[MOTIF_MATCH_SCORE_FOR_ABSENT_FROM_NEW_CATALOG]: (0.0, 0.0, 0.7, alpha),        #  blue
        }

        colors = [motif_match_colors.get(col, (0.5, 0.5, 0.5, alpha)) for col in plot_data_pivot.columns]
    else:
        colors = None

    ax = plot_data_pivot.plot(kind="barh", stacked=True, color=colors, ax=plt.gca())

    # make bars thicker
    for patch in ax.patches:
        patch.set_height(patch.get_height() * 1.5)
    plt.title(f"Does {args.catalog_name.replace('_', ' ')} capture TR loci from {new_catalog_name}?", pad=10)
    plt.xlabel("# of TR loci")

    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: "{:,.0f}".format(x)))

    # shift the plot to the right by 10% of the width of the plot
    plt.gca().set_position([0.22, 0.1, 0.75, 0.8])
    plt.gca().grid(axis="x", linestyle="-", linewidth=0.5, color="lightgray")

    # sort colors by motif_match_score
    handles, labels = plt.gca().get_legend_handles_labels()
    handle_label_map = dict(zip(labels, handles))
    ordered_handles = [handle_label_map[label] for label in desired_order if label in handle_label_map]
    ordered_labels = [label for label in desired_order if label in handle_label_map]

    plt.legend(ordered_handles, ordered_labels, title="", frameon=True)
    plt.gca().set_ylabel("")
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    output_path = f"{args.plot_output_dir}/{args.catalog_name}_vs_{new_catalog_name}.overlap_distribution_by_motif_match.png"
    plt.savefig(output_path)
    print(f"Wrote overlap distribution by motif match to {output_path}")
    plt.close()

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
            main_catalog_loci[chrom].add(intervaltree.Interval(start, max(start+1, end), data={
                "motif": motif,
                "canonical_motif": canonical_motif,
            }))

    return main_catalog_loci


def compute_motif_match_score(main_catalog_canonical_motif, other_catalog_canonical_motif):
    if main_catalog_canonical_motif == other_catalog_canonical_motif:
        return MOTIF_MATCH_SCORE_FOR_SAME_MOTIF
    elif len(main_catalog_canonical_motif) == len(other_catalog_canonical_motif):
        return MOTIF_MATCH_SCORE_FOR_SAME_MOTIF_LENGTH
    elif len(main_catalog_canonical_motif) < len(other_catalog_canonical_motif):
        return MOTIF_MATCH_SCORE_FOR_SHORTER_MOTIF_IN_MAIN_CATALOG
    elif len(main_catalog_canonical_motif) > len(other_catalog_canonical_motif):
        return MOTIF_MATCH_SCORE_FOR_LONGER_MOTIF_IN_MAIN_CATALOG
    else:
        raise Exception(f"Logic error")


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

    if main_catalog_interval.begin == start_0based and main_catalog_interval.end == end_1based:
        return OVERLAP_SCORE_FOR_EXACT_MATCH

    union_size = max(main_catalog_interval.end, end_1based) - min(main_catalog_interval.begin, start_0based)
    intersection_size = main_catalog_interval.overlap_size(start_0based, end_1based)
    if abs(intersection_size - union_size) <= 2*min_motif_size:
        return OVERLAP_SCORE_FOR_DIFF_2_REPEATS
    jaccard_similarity = intersection_size / union_size
    if jaccard_similarity > 2/3.0:
        return OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_ABOVE_0_66
    if jaccard_similarity > 0.5:
        return OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BETWEEN_0_5_AND_0_66
    if jaccard_similarity > 1/3.0:
        return OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BETWEEN_0_33_AND_0_5
    if jaccard_similarity > 0.2:
        return OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BETWEEN_0_2_AND_0_33

    return OVERLAP_SCORE_FOR_JACCARD_SIMILARITY_BELOW_0_2


def compute_match_summary(overlap_score, motif_match_score, new_catalog_motif, main_catalog_motif, main_catalog_name=None):
    high_overlap_score = overlap_score >= 5
    medium_or_high_overlap_score = overlap_score >= 2
    high_motif_similarity = (motif_match_score == 3 or (new_catalog_motif is not None and len(new_catalog_motif) > 6 and motif_match_score == 2))

    if main_catalog_motif is not None and new_catalog_motif is not None and not high_motif_similarity:
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


def compute_purity_of_optimal_motif_vs_length_of_random_sequence(min_sequence_length=14, max_sequence_length=1000, sequence_length_step=1, n_trials=1000, verbose=False):
    random.seed(1)
    sequence_length_to_max_purity_and_quality_score = collections.defaultdict(int)
    for length in range(min_sequence_length, max_sequence_length, sequence_length_step):
        max_purity = 0
        max_quality = 0
        for _ in range(n_trials): # N trials for each length
            random_sequence = ''.join(random.choices("ACGT", k=length))
            optimal_motif, motif_purity, quality = find_optimal_motif_length(random_sequence, length//2, verbose=False)
            max_purity = max(max_purity, motif_purity)
            max_quality = max(max_quality, quality)

        if verbose:
            print(f"{length:4d}bp   Max purity: {max_purity:.2f}     Max quality score: {max_quality:.2f}")

        sequence_length_to_max_purity_and_quality_score[length] = (max_purity, max_quality)


    return sequence_length_to_max_purity_and_quality_score

"""Results:
compute_purity_of_optimal_motif_vs_length_of_random_sequence(min_sequence_length=10, max_sequence_length=1001, sequence_length_step=10, n_trials=250, verbose=True)
  10bp   Max purity: 0.83     Max quality score: 0.67
  20bp   Max purity: 0.78     Max quality score: 0.56
  30bp   Max purity: 0.67     Max quality score: 0.39
  40bp   Max purity: 0.73     Max quality score: 0.50
  50bp   Max purity: 0.62     Max quality score: 0.37
  60bp   Max purity: 0.58     Max quality score: 0.32
  70bp   Max purity: 0.59     Max quality score: 0.35
  80bp   Max purity: 0.59     Max quality score: 0.35
  90bp   Max purity: 0.58     Max quality score: 0.33
 100bp   Max purity: 0.51     Max quality score: 0.27
 110bp   Max purity: 0.48     Max quality score: 0.22
 120bp   Max purity: 0.49     Max quality score: 0.24
 130bp   Max purity: 0.47     Max quality score: 0.22
 140bp   Max purity: 0.48     Max quality score: 0.23
 150bp   Max purity: 0.45     Max quality score: 0.20
 160bp   Max purity: 0.47     Max quality score: 0.23
 170bp   Max purity: 0.48     Max quality score: 0.23
 180bp   Max purity: 0.45     Max quality score: 0.20
 190bp   Max purity: 0.44     Max quality score: 0.19
 200bp   Max purity: 0.44     Max quality score: 0.19
 210bp   Max purity: 0.44     Max quality score: 0.20
 220bp   Max purity: 0.42     Max quality score: 0.17
 230bp   Max purity: 0.44     Max quality score: 0.18
 240bp   Max purity: 0.43     Max quality score: 0.19
 250bp   Max purity: 0.41     Max quality score: 0.17
 260bp   Max purity: 0.43     Max quality score: 0.18
 270bp   Max purity: 0.42     Max quality score: 0.16
 280bp   Max purity: 0.41     Max quality score: 0.15
 290bp   Max purity: 0.43     Max quality score: 0.18
 300bp   Max purity: 0.40     Max quality score: 0.15
 310bp   Max purity: 0.44     Max quality score: 0.19
 320bp   Max purity: 0.41     Max quality score: 0.16
 330bp   Max purity: 0.41     Max quality score: 0.16
 340bp   Max purity: 0.39     Max quality score: 0.13
 350bp   Max purity: 0.43     Max quality score: 0.19
 360bp   Max purity: 0.39     Max quality score: 0.14
 370bp   Max purity: 0.38     Max quality score: 0.13
 380bp   Max purity: 0.40     Max quality score: 0.15
 390bp   Max purity: 0.39     Max quality score: 0.14
 400bp   Max purity: 0.42     Max quality score: 0.17
 410bp   Max purity: 0.39     Max quality score: 0.15
 420bp   Max purity: 0.38     Max quality score: 0.13
 430bp   Max purity: 0.39     Max quality score: 0.14
 440bp   Max purity: 0.38     Max quality score: 0.13
 450bp   Max purity: 0.37     Max quality score: 0.12
 460bp   Max purity: 0.40     Max quality score: 0.15
 470bp   Max purity: 0.39     Max quality score: 0.14
 480bp   Max purity: 0.37     Max quality score: 0.12
 490bp   Max purity: 0.40     Max quality score: 0.15
 500bp   Max purity: 0.39     Max quality score: 0.13
 510bp   Max purity: 0.36     Max quality score: 0.12
 520bp   Max purity: 0.37     Max quality score: 0.12
 530bp   Max purity: 0.37     Max quality score: 0.12
 540bp   Max purity: 0.38     Max quality score: 0.12
 550bp   Max purity: 0.38     Max quality score: 0.13
 560bp   Max purity: 0.38     Max quality score: 0.13
 570bp   Max purity: 0.37     Max quality score: 0.12
 580bp   Max purity: 0.38     Max quality score: 0.13
 590bp   Max purity: 0.37     Max quality score: 0.12
 600bp   Max purity: 0.37     Max quality score: 0.11
 610bp   Max purity: 0.37     Max quality score: 0.12
 620bp   Max purity: 0.37     Max quality score: 0.12
 630bp   Max purity: 0.37     Max quality score: 0.12
 640bp   Max purity: 0.37     Max quality score: 0.12
 650bp   Max purity: 0.39     Max quality score: 0.13
 660bp   Max purity: 0.36     Max quality score: 0.11
 670bp   Max purity: 0.37     Max quality score: 0.12
 680bp   Max purity: 0.36     Max quality score: 0.12
 690bp   Max purity: 0.35     Max quality score: 0.10
 700bp   Max purity: 0.37     Max quality score: 0.12
 710bp   Max purity: 0.36     Max quality score: 0.11
 720bp   Max purity: 0.36     Max quality score: 0.11
 730bp   Max purity: 0.36     Max quality score: 0.11
 740bp   Max purity: 0.36     Max quality score: 0.12
 750bp   Max purity: 0.35     Max quality score: 0.11
 760bp   Max purity: 0.36     Max quality score: 0.11
 770bp   Max purity: 0.36     Max quality score: 0.11
 780bp   Max purity: 0.37     Max quality score: 0.12
 790bp   Max purity: 0.37     Max quality score: 0.12
 800bp   Max purity: 0.35     Max quality score: 0.10
 810bp   Max purity: 0.35     Max quality score: 0.10
 820bp   Max purity: 0.37     Max quality score: 0.12
 830bp   Max purity: 0.35     Max quality score: 0.10
 840bp   Max purity: 0.34     Max quality score: 0.09
 850bp   Max purity: 0.37     Max quality score: 0.12
 860bp   Max purity: 0.35     Max quality score: 0.10
 870bp   Max purity: 0.35     Max quality score: 0.10
 880bp   Max purity: 0.35     Max quality score: 0.10
 890bp   Max purity: 0.34     Max quality score: 0.09
 900bp   Max purity: 0.34     Max quality score: 0.09
 910bp   Max purity: 0.35     Max quality score: 0.10
 920bp   Max purity: 0.35     Max quality score: 0.10
 930bp   Max purity: 0.36     Max quality score: 0.11
 940bp   Max purity: 0.35     Max quality score: 0.10
 950bp   Max purity: 0.34     Max quality score: 0.08
 960bp   Max purity: 0.34     Max quality score: 0.09
 970bp   Max purity: 0.35     Max quality score: 0.10
 980bp   Max purity: 0.35     Max quality score: 0.10
 990bp   Max purity: 0.36     Max quality score: 0.10
1000bp   Max purity: 0.34     Max quality score: 0.08
 """




def find_optimal_motif_length_from_reference_fasta(reference_fasta, chrom, start_0based, end_1based, max_motif_length, verbose=False):
    if reference_fasta is None:
        return None, None

    chrom = f"chr{chrom.replace('chr', '')}"  # make sure chrom has "chr" prefix
    reference_sequence = reference_fasta[chrom][start_0based:end_1based]

    max_motif_length = min(max_motif_length, (end_1based - start_0based)//2)
    return find_optimal_motif_length(nucleotide_sequence=reference_sequence, max_motif_length=max_motif_length, verbose=verbose)


def compute_motif_length_quality(motif_length, motif_length_vs_motif_and_purity):
    optimal_motif_length_purities = []
    other_motif_length_purities = []
    for current_motif_length, (current_motif, purity) in motif_length_vs_motif_and_purity.items():
        if current_motif_length >= motif_length and current_motif_length % motif_length == 0:
            optimal_motif_length_purities.append(purity)
        else:
            other_motif_length_purities.append(purity)

    optimal_motif_length_mean_purity = sum(optimal_motif_length_purities) / len(optimal_motif_length_purities)
    if len(other_motif_length_purities) > 0:
        other_motif_lengths_mean_purity = sum(other_motif_length_purities) / len(other_motif_length_purities)
    else:
        other_motif_lengths_mean_purity = 0.31  # this is the base-line quality score upper-bound (the average of the distribution is 0.25 and the upper-bound appears to be ~0.31)

    quality_score = optimal_motif_length_mean_purity - other_motif_lengths_mean_purity

    return max(0, quality_score)


def find_optimal_motif_length(nucleotide_sequence, max_motif_length, verbose=False):
    """Scan different motif lengths from 1 to max_motif_length to find the one that produces the highest
    repeat purity with respect to the reference sequence at the given locus. For each motif length,
    this method finds the most frequent motif of that length within the reference sequence, then constructs
    a synthetic perfect repeat sequence of that motif (making it the same length as the reference sequence),
    then computes purity of that motif length as the fraction of bases in the reference sequence that match the
    previously constructed perfect repeat sequence of that motif.
    """


    #if verbose:
    #    print("--------------------------------")
    #    print(f"Sequence: {nucleotide_sequence}")

    motif_length_vs_motif_and_purity = {}
    for motif_length in range(1, max_motif_length+1):
        _, motif_length_purity, most_common_motif = get_purity_stats_for_locus(
            nucleotide_sequence=nucleotide_sequence, motif="A"*motif_length)
        if motif_length_purity is None:
            raise ValueError(f"Motif length {motif_length}bp purity is None for sequence: {nucleotide_sequence}")

        motif_length_vs_motif_and_purity[len(most_common_motif)] = (most_common_motif, motif_length_purity)
        #if verbose:
        #    print(f"{motif_length:3d}bp: {motif_length_purity:.2f}    {most_common_motif}")

    optimal_motif_length = max(
        motif_length_vs_motif_and_purity, key=lambda motif_length: (motif_length_vs_motif_and_purity[motif_length][1], -motif_length))
    optimal_motif, optimal_purity = motif_length_vs_motif_and_purity[optimal_motif_length]
    optimal_motif_length_quality_score = compute_motif_length_quality(optimal_motif_length, motif_length_vs_motif_and_purity)

    for shorter_motif_length in range(1, optimal_motif_length//2 + 1):
        if optimal_motif_length % shorter_motif_length != 0:
            # only consider shorter motif lengths that are factors of the current optimal motif length
            continue

        shorter_motif_length_quality_score = compute_motif_length_quality(shorter_motif_length, motif_length_vs_motif_and_purity)
        if optimal_motif_length_quality_score - shorter_motif_length_quality_score < 0.03:
            original_optimal_motif = optimal_motif
            optimal_motif, optimal_purity = motif_length_vs_motif_and_purity[shorter_motif_length]
            if verbose:
                print(f"Replacing motif length {original_optimal_motif} ({optimal_motif_length}bp) which had quality {optimal_motif_length_quality_score} "
                      f"with {optimal_motif} ({shorter_motif_length}bp) which has quality {shorter_motif_length_quality_score} ")
                      #f"in sequence: {nucleotide_sequence}")
            optimal_motif_length_quality_score = shorter_motif_length_quality_score
            break

    simplified_optimal_motif, _, _ = find_repeat_unit_without_allowing_interruptions(optimal_motif, allow_partial_repeats=False)
    #if simplified_optimal_motif != optimal_motif:
    #    print(f"WARNING: Simplified optimal motif {simplified_optimal_motif} != optimal motif {optimal_motif} for sequence: {nucleotide_sequence}")  # this happens occasionally due to edge cases

    if verbose:
        print(f"Optimal motif: {len(optimal_motif)}bp   purity: {optimal_purity:.2f}   quality: {optimal_motif_length_quality_score}   (null quality is: ) tried "
              f"{len(motif_length_vs_motif_and_purity)} motif lengths, their average purity was: {sum([x[1] for x in motif_length_vs_motif_and_purity.values()]) / len(motif_length_vs_motif_and_purity):.2f}")

    return simplified_optimal_motif, optimal_purity, optimal_motif_length_quality_score


def get_purity_stats_for_locus(nucleotide_sequence, motif):
    _, motif_purity, _ = compute_sequence_purity_stats(nucleotide_sequence, motif, include_partial_repeats=True)

    # slice the reference sequence into subsequences of length len(motif) and then get the most common motif
    end_index = len(nucleotide_sequence) - len(nucleotide_sequence) % len(motif)
    if end_index == 0:
        return None, None, None

    sliced_motif_list = [nucleotide_sequence[i:i+len(motif)] for i in range(0, end_index, len(motif))]
    if len(sliced_motif_list) == 0:
        return None, None, None

    most_common_motif = collections.Counter(sliced_motif_list).most_common(1)[0][0]
    # remove the first occurrence of the most common motif from the sliced nucleotide sequence list
    sliced_motif_list.remove(most_common_motif)
    if len(sliced_motif_list) == 0:
        return None, None, None

    remaining_nucleotide_sequence = "".join(sliced_motif_list)
    _, motif_length_purity, _ = compute_sequence_purity_stats(remaining_nucleotide_sequence, most_common_motif, include_partial_repeats=True)

    return motif_purity, motif_length_purity, most_common_motif


def get_purity_stats_for_locus_from_reference_fasta(reference_fasta, chrom, start_0based, end_1based, motif):
    if reference_fasta is None:
        return None, None, None

    chrom = f"chr{chrom.replace('chr', '')}"  # make sure chrom has "chr" prefix
    reference_sequence = reference_fasta[chrom][start_0based:end_1based]

    return get_purity_stats_for_locus(nucleotide_sequence=reference_sequence, motif=motif)


def compute_optimal_motif_match_score(
        optimal_motif,
        new_catalog_name,
        new_catalog_motif,
        main_catalog_name,
        main_catalog_motif,
):
    if main_catalog_motif is None:
        return f"absent from {main_catalog_name} catalog"

    if new_catalog_motif is None:
        return f"absent from {new_catalog_name} catalog"

    main_catalog_canonical_motif = compute_canonical_motif(main_catalog_motif)
    new_catalog_canonical_motif = compute_canonical_motif(new_catalog_motif)
    optimal_canonical_motif = compute_canonical_motif(optimal_motif)
    if len(new_catalog_motif) == len(main_catalog_motif):
        common_motif_length = len(new_catalog_motif)
        if main_catalog_canonical_motif == optimal_canonical_motif and new_catalog_canonical_motif == optimal_canonical_motif:
            return "both match optimal motif"
        elif main_catalog_canonical_motif == optimal_canonical_motif:
            return f"both match optimal motif length, {main_catalog_name} matches optimal motif"
        elif new_catalog_canonical_motif == optimal_canonical_motif:
            return f"both match optimal motif length, {new_catalog_name} matches optimal motif"
        else:
            # both differ from the optimal motif seq
            if common_motif_length == len(optimal_motif) and main_catalog_canonical_motif == new_catalog_canonical_motif:
                return f"both match optimal motif length and have same motif, but motif differs from optimal motif"
            elif common_motif_length == len(optimal_motif) and main_catalog_canonical_motif != new_catalog_canonical_motif:
                return f"both match optimal motif length, but all 3 motifs differ"
            else:
                if common_motif_length > len(optimal_motif):
                    return f"same motif length in both, but optimal motif length is shorter"
                else:
                    return f"same motif length in both, but optimal motif length is longer"
    elif len(main_catalog_motif) == len(optimal_motif):
        if main_catalog_canonical_motif == optimal_canonical_motif:
            return f"only {main_catalog_name} has optimal motif"
        else:
            return f"only {main_catalog_name} has optimal motif length"
    elif len(new_catalog_motif) == len(optimal_motif):
        if new_catalog_canonical_motif == optimal_canonical_motif:
            return f"only {new_catalog_name} has optimal motif"
        else:
            return f"only {new_catalog_name} has optimal motif length"
    else:
        if len(optimal_motif) > max(len(main_catalog_motif), len(new_catalog_motif)):
            return f"all three differ in length, optimal motif is longest"
        elif len(optimal_motif) < min(len(main_catalog_motif), len(new_catalog_motif)):
            return f"all three differ in length, optimal motif is shortest"
        else:
            return f"all three differ in length, optimal motif length is in the middle"


def compare_loci(
        main_catalog_loci,
        new_catalog,
        main_catalog_name="catalog1",
        new_catalog_name="catalog2",
        reference_fasta=None,
        write_loci_absent_from_new_catalog=False,
):
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
    new_catalog_duplicate_detector = set()
    new_catalog_duplicate_counter = 0
    new_catalog_total_counter = 0
    main_catalog_loci_that_overlap_new_catalog = set()
    for fields in tqdm.tqdm(fields_iterator, unit=" records", unit_scale=True):
        new_catalog_total_counter += 1
        chrom = fields[0].replace("chr", "")
        start_0based = int(fields[1])
        end_1based = int(fields[2])
        if len(fields) < 4:
            motif = None
            simplified_motif = None
            canonical_motif = None
        else:
            motif = fields[3]
            simplified_motif, _, _ = find_repeat_unit_without_allowing_interruptions(motif, allow_partial_repeats=False)
            canonical_motif = compute_canonical_motif(simplified_motif)

        key = (chrom, start_0based, end_1based, canonical_motif)
        if key in new_catalog_duplicate_detector:
            new_catalog_duplicate_counter += 1
            continue
        new_catalog_duplicate_detector.add(key)

        overlap_score = 0
        motif_match_score = 0
        closest_match_main_catalog_interval = None

        # find the overlapping main catalog locus that overlaps the largest percentage of the new locus and is closest in size
        main_catalog_overlap = main_catalog_loci[chrom].overlap(start_0based, end_1based)
        for main_catalog_interval in main_catalog_overlap:
            main_catalog_canonical_motif = main_catalog_interval.data["canonical_motif"]
            if write_loci_absent_from_new_catalog:
                main_catalog_loci_that_overlap_new_catalog.add((chrom, main_catalog_interval.begin, main_catalog_interval.end, main_catalog_canonical_motif))

            current_overlap_score = compute_overlap_score(
                main_catalog_interval,
                start_0based,
                end_1based,
                min_motif_size=min(len(canonical_motif), len(main_catalog_canonical_motif)) if motif else len(main_catalog_canonical_motif))

            if current_overlap_score < overlap_score:
                continue

            if motif:
                current_motif_match_score = compute_motif_match_score(main_catalog_canonical_motif, canonical_motif)
            else:
                current_motif_match_score = MOTIF_MATCH_SCORE_FOR_SAME_MOTIF

            if (current_overlap_score, current_motif_match_score) < (overlap_score, motif_match_score):
                continue

            overlap_score = current_overlap_score
            motif_match_score = current_motif_match_score
            closest_match_main_catalog_interval = main_catalog_interval

        main_catalog_canonical_motif = closest_match_main_catalog_interval.data["canonical_motif"] if closest_match_main_catalog_interval is not None else None
        was_match_found = compute_match_summary(overlap_score, motif_match_score, canonical_motif, main_catalog_canonical_motif, main_catalog_name)

        reference_repeat_count = (end_1based - start_0based) // len(canonical_motif) if canonical_motif is not None else None

        motif_purity, motif_length_purity, _ = get_purity_stats_for_locus_from_reference_fasta(reference_fasta, chrom, start_0based, end_1based, motif)
        optimal_motif, optimal_motif_purity, optimal_motif_quality_score = find_optimal_motif_length_from_reference_fasta(
            reference_fasta,
            chrom,
            start_0based,
            end_1based,
            max_motif_length=2 * (max(len(motif), len(main_catalog_canonical_motif)) if main_catalog_canonical_motif else len(motif)))

        optimal_motif_match = compute_optimal_motif_match_score(
            optimal_motif,
            new_catalog_name,
            motif,
            main_catalog_name,
            main_catalog_motif = closest_match_main_catalog_interval.data["motif"] if closest_match_main_catalog_interval is not None else None,
        )

        output_row = {
            "chrom": chrom,

            "start_0based": start_0based,
            "end_1based": end_1based,
            "motif": motif,
            "simplified_motif": simplified_motif,
            "canonical_motif": canonical_motif,
            "reference_region": f"{chrom}:{start_0based}-{end_1based}",
            "reference_repeat_count": reference_repeat_count,
            "reference_region_size": end_1based - start_0based,

            f"{new_catalog_name}_start_0based": start_0based,
            f"{new_catalog_name}_end_1based": end_1based,
            f"{new_catalog_name}_motif": motif,
            f"{new_catalog_name}_simplified_motif": simplified_motif,
            f"{new_catalog_name}_canonical_motif": canonical_motif,
            f"{new_catalog_name}_reference_region": f"{chrom}:{start_0based}-{end_1based}",
            f"{new_catalog_name}_reference_repeat_count": reference_repeat_count,
            f"{new_catalog_name}_reference_region_size": end_1based - start_0based,
            f"{new_catalog_name}_purity_of_motif": motif_purity,
            f"{new_catalog_name}_purity_of_motif_length": motif_length_purity,

            "optimal_motif": optimal_motif,
            "optimal_canonical_motif": compute_canonical_motif(optimal_motif),
            "optimal_motif_length": len(optimal_motif),
            "optimal_motif_purity": optimal_motif_purity,
            "optimal_motif_quality_score": optimal_motif_quality_score,
            "optimal_motif_match": optimal_motif_match,

            "overlap_score": overlap_score,
            "motif_match_score": motif_match_score,
            "overlap": OVERLAP_SCORE_MAP[overlap_score],
            "motif_match": MOTIF_MATCH_SCORE_MAP[motif_match_score],
            "match_found?": was_match_found,
        }

        if closest_match_main_catalog_interval is not None:
            closest_match_start_0based = closest_match_main_catalog_interval.begin
            closest_match_end_1based = closest_match_main_catalog_interval.end
            closest_match_motif = closest_match_main_catalog_interval.data["motif"]
            closest_match_canonical_motif = closest_match_main_catalog_interval.data["canonical_motif"]
            closest_match_reference_region = f"{chrom}:{closest_match_start_0based}-{closest_match_end_1based}"
            closest_match_motif_purity, closest_match_motif_length_purity, _ = get_purity_stats_for_locus_from_reference_fasta(
                reference_fasta, chrom, closest_match_start_0based, closest_match_end_1based, closest_match_motif)

            output_row.update({
                f"{main_catalog_name}_start_0based": closest_match_start_0based,
                f"{main_catalog_name}_end_1based": closest_match_end_1based,
                f"{main_catalog_name}_motif": closest_match_motif,
                f"{main_catalog_name}_canonical_motif": closest_match_canonical_motif,
                f"{main_catalog_name}_reference_region": closest_match_reference_region,
                f"{main_catalog_name}_reference_repeat_count": (closest_match_end_1based - closest_match_start_0based) // len(closest_match_canonical_motif),
                f"{main_catalog_name}_reference_region_size": closest_match_end_1based - closest_match_start_0based,
                f"{main_catalog_name}_purity_of_motif": closest_match_motif_purity,
                f"{main_catalog_name}_purity_of_motif_length": closest_match_motif_length_purity,
            })

        output_rows.append(output_row)

    if new_catalog_duplicate_counter:
        print(f"WARNING: Skipped {new_catalog_duplicate_counter:,d} out of {new_catalog_total_counter:,d} loci in {new_catalog_name} because they were exact duplicates of a previous locus based on (chrom, start, end, motif)")

    if write_loci_absent_from_new_catalog:
        for chrom, interval_tree in tqdm.tqdm(main_catalog_loci.items(), unit=" chromosomes", unit_scale=True):
            for main_catalog_interval in interval_tree:
                main_catalog_canonical_motif = main_catalog_interval.data["canonical_motif"]
                if (chrom, main_catalog_interval.begin, main_catalog_interval.end, main_catalog_canonical_motif) in main_catalog_loci_that_overlap_new_catalog:
                    continue

                main_catalog_motif = main_catalog_interval.data["motif"]
                main_catalog_reference_region = f"{chrom}:{main_catalog_interval.begin}-{main_catalog_interval.end}"
                main_catalog_reference_repeat_count = (main_catalog_interval.end - main_catalog_interval.begin) // len(main_catalog_canonical_motif)
                main_catalog_motif_purity, main_catalog_motif_length_purity, _ = get_purity_stats_for_locus_from_reference_fasta(
                    reference_fasta, chrom, main_catalog_interval.begin, main_catalog_interval.end, main_catalog_motif)
                output_rows.append({
                    "chrom": chrom,

                    "start_0based": main_catalog_interval.begin,
                    "end_1based": main_catalog_interval.end,
                    "motif": main_catalog_motif,
                    "simplified_motif": main_catalog_canonical_motif,
                    "canonical_motif": main_catalog_canonical_motif,
                    "reference_region": main_catalog_reference_region,
                    "reference_repeat_count": main_catalog_reference_repeat_count,
                    "reference_region_size": main_catalog_interval.end - main_catalog_interval.begin,

                    "overlap_score": -1,
                    "motif_match_score": -1,
                    "overlap": OVERLAP_SCORE_MAP[-1],
                    "motif_match": MOTIF_MATCH_SCORE_MAP[-1],
                    "match_found?": NO_MATCH_FOUND,
                    f"{main_catalog_name}_start_0based": main_catalog_interval.begin,
                    f"{main_catalog_name}_end_1based": main_catalog_interval.end,
                    f"{main_catalog_name}_motif": main_catalog_motif,
                    f"{main_catalog_name}_canonical_motif": main_catalog_canonical_motif,
                    f"{main_catalog_name}_reference_region": main_catalog_reference_region,
                    f"{main_catalog_name}_reference_repeat_count": main_catalog_reference_repeat_count,
                    f"{main_catalog_name}_reference_region_size": main_catalog_interval.end - main_catalog_interval.begin,
                    f"{main_catalog_name}_purity_of_motif": main_catalog_motif_purity,
                    f"{main_catalog_name}_purity_of_motif_length": main_catalog_motif_length_purity,
                })

    if f is not None:
        f.close()

    return pd.DataFrame(output_rows)



if __name__ == "__main__":
    main()