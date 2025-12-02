"""Generate supplementary figure that provides an overview of the catalog and its core properties:

- Motif Size distribution  (colored by specific motifs?)
- Number of repeats in the reference  (colored by purity?)
- Gene Regions
- Polymorphism rates
"""


import argparse
import collections
import ijson
import gzip
import matplotlib.pyplot as plt
import os
import pandas as pd
import re
import seaborn as sns
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

# column names:
"""
['chrom', 'start_0based', 'end_1based', 'ReferenceRegion', 'LocusId',
 'ReferenceMotif', 'MotifSize', 'CanonicalMotif',
 'NumRepeatsInReference', 'ReferenceRepeatPurity', 'NsInFlanks',
 'TRsInRegion', 'Source', 'FlanksAndLocusMappability',
 'KnownDiseaseAssociatedLocus', 'KnownDiseaseAssociatedMotif',
 'GencodeGeneRegion', 'GencodeGeneName', 'GencodeGeneId',
 'GencodeTranscriptId', 'RefseqGeneRegion', 'RefseqGeneName',
 'RefseqGeneId', 'RefseqTranscriptId', 'ManeGeneRegion', 'ManeGeneName',
 'ManeGeneId', 'ManeTranscriptId', 'TenK10K_AlleleHistogram',
 'TenK10K_ModeAllele', 'TenK10K_Stdev', 'TenK10K_Median',
 'TenK10K_99thPercentile', 'HPRC100_AlleleHistogram',
 'HPRC100_ModeAllele', 'HPRC100_Stdev', 'HPRC100_Median',
 'HPRC100_99thPercentile', 'AoU1027_ModeAllele', 'AoU1027_Stdev',
 'AoU1027_Median', 'AoU1027_99thPercentile', 'AoU1027_MaxAllele',
 'AoU1027_UniqueAlleles', 'AoU1027_OE_Length',
 'AoU1027_OE_LengthPercentile', 'AlleleFrequenciesFromIllumina174k',
 'StdevFromIllumina174k', 'AlleleFrequenciesFromT2TAssemblies',
 'StdevFromT2TAssemblies', 'VariationCluster',
 'VariationClusterSizeDiff']
"""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-n", type=int, help="Number of loci to process")
    p.add_argument("--catalog-name", default="TRExplorer_v1")
    p.add_argument("--catalog-path",
                   #default="../results__2025-09-04/release_draft_2025-09-04/repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz",
                   #default="../results__2025-11-03/release_draft_2025-11-03/repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz",
                   default="../results__2025-09-04/release_draft_2025-09-04/repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz",
                   help="Catalog path")
    args = p.parse_args()

    if not os.path.isfile(args.catalog_path):
        p.error(f"File not found: {args.catalog_path}")


    motif_size_distribution = collections.Counter()
    reference_repeat_count_distribution = collections.Counter()
    gene_region_distribution = collections.Counter()
    for key in ['CDS', "5' UTR", "3' UTR", 'exon', 'intron', 'promoter', 'intergenic']:
        gene_region_distribution[key] = 0
        
    print(f"Parsing {args.catalog_path} to json")
    fopen = gzip.open if args.catalog_path.endswith("gz") else open
    with fopen(args.catalog_path, "rt") as f:
        # use ijson to parse the json
        for i, record in tqdm.tqdm(enumerate(ijson.items(f, "item")), unit=" records", unit_scale=True):
            if args.n is not None and i >= args.n:
                break

            motif_size = len(record["CanonicalMotif"])
            motif_size_distribution[motif_size] += 1
            reference_repeat_count_distribution[record["NumRepeatsInReference"]] += 1
            gene_region_distribution[record["GencodeGeneRegion"]] += 1


    # bin values between 7 and 24 as "7-24" and 25+ as "25+"
    motif_size_distribution_binned = {f"{motif_size}": 0 for motif_size in range(1, 7)}
    motif_size_distribution_binned["7-24"] = 0
    motif_size_distribution_binned["25+"] = 0
    for motif_size, count in motif_size_distribution.items():
        if motif_size < 7:
            motif_size_distribution_binned[f"{motif_size}"] += count
        elif motif_size < 25:
            motif_size_distribution_binned["7-24"] += count
        else:
            motif_size_distribution_binned["25+"] += count
    motif_size_distribution = motif_size_distribution_binned

    catalog_name = args.catalog_name.replace(" ", "_")
    plt.figure(figsize=(12, 8))
    sns.barplot(x=list(motif_size_distribution.keys()), y=list(motif_size_distribution.values()), color="cornflowerblue")
    plt.xticks(rotation=45)
    plt.xlabel("Motif size (bp)", fontsize=16)
    plt.ylabel("# of TRs", fontsize=16)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().set_position([0.15, 0.15, 0.8, 0.8])
    plt.gca().xaxis.labelpad = 15
    plt.gca().yaxis.labelpad = 15
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: "{:,.0f}".format(x)))
    plt.title(f"{catalog_name.replace('_', ' ')} motif size distribution", fontsize=18)
    plt.gca().tick_params(axis='both', which='major', labelsize=14)
    plt.savefig(f"{catalog_name}.motif_size_distribution.png")
    print(f"Wrote motif size distribution to {catalog_name}.motif_size_distribution.png")
    plt.close()

    # bin values between 1 and 10 as "1-10" and 11+ as "11+"
    reference_repeat_count_distribution_binned = {f"{repeat_count}": 0 for repeat_count in range(0, 11)}
    reference_repeat_count_distribution_binned["11+"] = 0
    for repeat_count, count in reference_repeat_count_distribution.items():
        if repeat_count < 11:
            reference_repeat_count_distribution_binned[f"{repeat_count}"] += count
        else:
            reference_repeat_count_distribution_binned["11+"] += count
    reference_repeat_count_distribution = reference_repeat_count_distribution_binned

    catalog_name = catalog_name.replace(" ", "_")
    plt.figure(figsize=(12, 8))
    sns.barplot(x=list(reference_repeat_count_distribution.keys()), y=list(reference_repeat_count_distribution.values()), color="cornflowerblue")
    plt.xticks(rotation=45)
    plt.xlabel("Repeats in hg38", fontsize=16)
    plt.ylabel("# of TRs", fontsize=16)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().set_position([0.15, 0.15, 0.8, 0.8])
    plt.gca().xaxis.labelpad = 15
    plt.gca().yaxis.labelpad = 15
    plt.title(f"{catalog_name.replace('_', ' ')} reference repeat count distribution", fontsize=18)
    plt.gca().tick_params(axis='both', which='major', labelsize=14)
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: "{:,.0f}".format(x)))
    plt.savefig(f"{catalog_name}.reference_repeat_count_distribution.png")
    print(f"Wrote reference repeat count distribution to {catalog_name}.reference_repeat_count_distribution.png")
    plt.close()


    # plot the gene region distribution
    plt.figure(figsize=(12, 8))
    sns.barplot(x=list(gene_region_distribution.keys()), y=list(gene_region_distribution.values()), color="cornflowerblue")
    plt.xticks(rotation=45)
    plt.xlabel("Gene region", fontsize=16)
    plt.ylabel("# of TRs", fontsize=16)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().set_position([0.15, 0.20, 0.8, 0.7])


    plt.gca().xaxis.labelpad = 15
    plt.gca().yaxis.labelpad = 15
    plt.title(f"{catalog_name.replace('_', ' ')} gene region distribution", fontsize=18)
    plt.gca().tick_params(axis='both', which='major', labelsize=14)
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: "{:,.0f}".format(x)))
    plt.savefig(f"{catalog_name}.gene_region_distribution.png")
    print(f"Wrote gene region distribution to {catalog_name}.gene_region_distribution.png")
    plt.close()



if __name__ == "__main__":
    main()