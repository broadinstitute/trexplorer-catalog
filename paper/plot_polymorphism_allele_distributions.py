import collections
import pprint as pp
import glob
import ijson
import json
import matplotlib.pyplot as plt
import os
import pandas as pd
import gzip
import re
import tqdm

from str_analysis.utils.misc_utils import parse_interval

stats_lookup = collections.defaultdict(dict)

annotated_catalog_json_path = "~/code/tandem-repeat-explorer/downloads/TR_catalog.5591917_loci.20251209_110637.json.gz"
with gzip.open(os.path.expanduser(annotated_catalog_json_path), "rt") as f:  
    for i, record in tqdm.tqdm(enumerate(ijson.items(f, "item")), total=5_700_000, unit=" records", unit_scale=True):
        locus_id = record["LocusId"]

        stats = {}
        if "HPRC256_AlleleHistogram" in record:
            stats.update({
                "HPRC256_NumAlleles": int(record["HPRC256_UniqueAlleles"]),
                "HPRC256_Stdev": float(record["HPRC256_Stdev"]),
            })

            if record["HPRC256_AlleleHistogram"].count(",") + 1 != int(record["HPRC256_UniqueAlleles"]):
                print(f"WARNING: For locus {record['LocusId']}, HPRC256_AlleleHistogram has {record['HPRC256_AlleleHistogram'].count(',') + 1} alleles, but HPRC256_UniqueAlleles = {int(record['HPRC256_UniqueAlleles'])}")
            
        if "AoU1027_Stdev" in record:
            stats.update({
                "AoU1027_Stdev": float(record["AoU1027_Stdev"]),
                "AoU1027_NumAlleles": int(record["AoU1027_UniqueAlleles"]),
            })

            s = {int(record[k]) for k in ["AoU1027_MinAllele", "AoU1027_ModeAllele", "AoU1027_MaxAllele"]}
            if len(s) > 1 and record["AoU1027_Stdev"] == 0:
                print(f"WARNING: For locus {record['LocusId']}, {s}, ModeAllele = {record['AoU1027_ModeAllele']}, Median = {record['AoU1027_Median']}, 99thPercentile = {record['AoU1027_99thPercentile']}, MaxAllele == {record['AoU1027_MaxAllele']}, but AoU Stdev == 0")
            if int(record["AoU1027_UniqueAlleles"]) < len(s):
                print(f"WARNING: For locus {record['LocusId']}, {s}, AoU1027_UniqueAlleles = {record['AoU1027_UniqueAlleles']}, but {s}")

        if not stats:
            continue

        chrom, start_0based, end = parse_interval(record["ReferenceRegion"])
        stats_lookup[locus_id] = {
            "LocusId": locus_id,
            "Motif": record["CanonicalMotif"],
            "MotifSize": record["MotifSize"],
            "Source": record["Source"],
            #"VariationCluster": record["VariationCluster"],
            "NumRepeatsInReference": int(record["NumRepeatsInReference"]),
            "NumBasesInReference": end - start_0based,
        }
        stats_lookup[locus_id].update(stats)


output_table = "HPRC256_and_AoU1027_polymorphism_stats.tsv"
df = pd.DataFrame(stats_lookup.values())
df.to_csv(output_table, sep="\t", index=False, header=True)
print(f"Wrote {len(df):,d} rows to {output_table}")


downsampled_HPRC256_table_paths = []
for path in glob.glob(os.path.expanduser("~/code/tandem-repeat-explorer/data-prep/hprc_lps.*.grouped_by_locus_and_motif.with_biallelic_histogram.tsv.gz")):
    """
    Example file name: hprc_lps.90_samples.grouped_by_locus_and_motif.with_biallelic_histogram.tsv.gz
    Example row:
    $1                locus_id : 1-61871-61880-T
    $2                   motif : T
    $3   allele_size_histogram : 9x:144,10x:30
    $4     biallelic_histogram : 9/9:60,9/10:24,10/10:3
    $5             mode_allele : 9
    $6                    mean : 9.172
    $7                   stdev : 0.378
    $8                  median : 9
    $9         99th_percentile : 10
    """  
    if os.path.basename(path) == "hprc_lps.2025_05.grouped_by_locus_and_motif.with_biallelic_histogram.tsv.gz":
        number_of_samples = 256
    else:
        number_of_samples = int(re.search("hprc_lps[.](.*)_samples[.]grouped_by_locus_and_motif", path).group(1))

    downsampled_HPRC256_table_paths.append((number_of_samples, path))

downsampled_stats_rows = []
for number_of_samples, path in sorted(downsampled_HPRC256_table_paths):
    df = pd.read_table(path)
    print(f"Read {len(df):,d} rows from {path}")
    for locus_id, motif, hist in zip(df.locus_id, df.motif, df.allele_size_histogram):
        downsampled_stats_rows.append({
            "LocusId": locus_id,
            "MotifSize": len(motif),
            "Motif": motif,
            "NumberOfSamples": number_of_samples,
            "NumberOfAlleles": hist.count(",") + 1,
        })

df = pd.DataFrame(downsampled_stats_rows)
df.to_csv("HPRC256_downsampled_allele_counts.tsv.gz", sep="\t", index=False, header=True)
print(f"Wrote {len(df):,d} rows to HPRC256_downsampled_allele_counts.tsv.gz")

# AoU sanity checks: 
#  Median <= MaxAllele
#  Median <= 99thPercentile
#  99thPercentile <= MaxAllele
#  Mode <= MaxAllele
#  UniqueAlleles >= len(set({Mode, Median, 99thPercentile, MaxAllele}))


"""
{'AlleleFrequenciesFromIllumina174k': '3x:298,5x:4693,6x:7',
 'AlleleFrequenciesFromT2TAssemblies': '2x:1,3x:155',
 'AoU1027_99thPercentile': '4',
 'AoU1027_MaxAllele': '4',
 'AoU1027_Median': '4',
 'AoU1027_ModeAllele': '1',
 'AoU1027_OE_Length': Decimal('0.3875704817720314'),
 'AoU1027_OE_LengthPercentile': Decimal('0.2488102853902814'),
 'AoU1027_Stdev': Decimal('0.02841416741671275'),
 'AoU1027_UniqueAlleles': '7',
 'CanonicalMotif': 'AAC',
 'FlanksAndLocusMappability': Decimal('0.02'),
 'GencodeGeneId': 'ENSG00000225880',
 'GencodeGeneName': 'LINC00115',
 'GencodeGeneRegion': 'intron',
 'GencodeTranscriptId': 'ENST00000506640',
 'HPRC256_99thPercentile': '4',
 'HPRC256_AlleleHistogram': '4x:200',
 'HPRC256_Median': '4',
 'HPRC256_ModeAllele': '4',
 'HPRC256_Stdev': 0,
 'KnownDiseaseAssociatedMotif': 'AGC',
 'LocusId': '1-747743-747757-AAC',
 'LocusStructure': '(AAC)*',
 'ManeGeneId': 'ENSG00000284662',
 'ManeGeneName': 'OR4F16',
 'ManeGeneRegion': 'intergenic',
 'ManeTranscriptId': 'ENST00000332831',
 'MotifSize': '3',
 'NsInFlanks': '0',
 'NumRepeatsInReference': '4',
 'ReferenceMotif': 'AAC',
 'ReferenceRegion': 'chr1:747743-747757',
 'ReferenceRepeatPurity': 1,
 'RefseqGeneId': 'LOC100288069',
 'RefseqGeneName': 'LOC100288069',
 'RefseqGeneRegion': 'intron',
 'RefseqTranscriptId': 'NR_168328',
 'Source': 'PerfectRepeatsInReference',
 'StdevFromIllumina174k': Decimal('0.475'),
 'StdevFromT2TAssemblies': Decimal('0.08'),
 'TRsInRegion': '1',
 'TenK10K_99thPercentile': '10',
 'TenK10K_AlleleHistogram': '1x:7,2x:1,3x:3,4x:16,5x:3647,6x:50,7x:60,8x:7,9x:5,10x:9,11x:7,12x:4,13x:4,14x:7,15x:2,16x:1,18x:1,21x:1,28x:2,29x:1,41x:1',
 'TenK10K_Median': '5',
 'TenK10K_ModeAllele': '5',
 'TenK10K_Stdev': Decimal('1.217'),
 'VariationCluster': '1:735622-735641',
 'VariationClusterSizeDiff': '10'}
 """

from glob import glob
import os
import re


output_rows = []
for path in sorted(glob(os.path.expanduser("~/code/tandem-repeat-explorer/data-prep/hprc-lps_2025-12-06/hprc_lps.2025_12.per_locus_and_motif.only_AFR.[0-9]0*_samples.tsv.gz"))):
    number_of_samples = int(re.search("[.]([0-9]+)_samples.tsv.gz", path).group(1))
    print(f"{number_of_samples:,} samples: {path}")
    # Read table and count the number of rows with stdev > 0.2
    df = pd.read_table(path)
    for allele_count_threshold in range(2, 6):
        count = len(df[df.unique_alleles >= allele_count_threshold])
        print(f"{os.path.basename(path)}: {count:,} rows with >= {allele_count_threshold} unique_alleles")
        output_rows.append({
            "filename": os.path.basename(path),
            "num_samples": number_of_samples,
            "num_loci": count,
            "allele_count_threshold": allele_count_threshold,
        })

df = pd.DataFrame(output_rows)
df.to_csv("downsampled_numbers_of_polymorphic_loci.tsv", sep="\t", index=False, header=True)



"""
Example row:

locus_id                 1-19175-19184-TCC
motif                                  TCC
allele_size_histogram                3x:60
biallelic_histogram                 3/3:30
min_allele                               3
mode_allele                              3
mean                                   3.0
stdev                                  0.0
median                                   3
99th_percentile                          3
max_allele                               3
unique_alleles                           1
num_called_alleles                      60

"""