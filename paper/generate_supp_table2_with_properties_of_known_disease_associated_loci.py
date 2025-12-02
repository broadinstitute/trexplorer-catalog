"""This script generates Supp. Table 2 in the paper"""

import pandas as pd

df = pd.read_json("paper/TR_catalog.63_loci.20251010_173909.json.gz")

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

df["Motif Size"] = df["CanonicalMotif"].str.len()
df["VariationClusterSizeDiff"] = df["VariationClusterSizeDiff"].apply(
    lambda x : "" if pd.isna(x) else f"+{int(x)} bp")

df["TRsInRegion"] -= 1

# replace "CDS" with "coding"
df["GencodeGeneRegion"] =  df["GencodeGeneRegion"].replace("CDS", "coding")

df.sort_values([
    "MotifSize", "CanonicalMotif", "GencodeGeneRegion", "KnownDiseaseAssociatedLocus",
], inplace=True)


rename_dict = {
  "KnownDiseaseAssociatedLocus": "Name",
  "ReferenceRegion": "Coordinates (hg38, 0-based)",
  "CanonicalMotif": "Motif (normalized)",
  "GencodeGeneRegion": "Gene Region",
  "HPRC100_Stdev": "Polymorphism (HPRC)",
  "FlanksAndLocusMappability": "Mappability",
  "VariationClusterSizeDiff": "Variation Cluster",
  "TRsInRegion": "Other TRs Nearby",
}

df = df[list(rename_dict.keys())]
df.rename(columns=rename_dict, inplace=True)

print(df.columns)
print(df)
df.to_csv("paper/table_of_properties_of_known_disease_associated_loci.tsv", sep="\t", index=False)

print(f"Wrote {len(df):,d} loci to paper/table_of_properties_of_known_disease_associated_loci.tsv")
#%%