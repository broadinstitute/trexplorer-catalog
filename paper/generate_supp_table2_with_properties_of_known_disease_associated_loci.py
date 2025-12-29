"""This script generates Supp. Table 2 in the paper"""

import pandas as pd

df = pd.read_json("TR_catalog.63_loci.20251229_005446.json.gz")

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
 'TenK10K_99thPercentile', 'HPRC256_AlleleHistogram',
 'HPRC256_ModeAllele', 'HPRC256_Stdev', 'HPRC256_Median',
 'HPRC256_99thPercentile', 'AoU1027_ModeAllele', 'AoU1027_Stdev',
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
    "MotifSize", "CanonicalMotif", "KnownDiseaseAssociatedLocus", "GencodeGeneRegion", "KnownDiseaseAssociatedLocus",
], inplace=True)


rename_dict = {
  "KnownDiseaseAssociatedLocus": "Name",
  "ReferenceRegion": "Coordinates (hg38, 0-based)",
  "CanonicalMotif": "Motif (normalized)",
  "MotifSize": "Motif Size",
  "GencodeGeneRegion": "Gene Region",
  "HPRC256_Stdev": "Polymorphism (HPRC256 stdev)",
  "AoU1027_Stdev": "Polymorphism (AoU1027 stdev)",
  "FlanksAndLocusMappability": "Mappability",
  "VariationClusterSizeDiff": "Variation Cluster",
  "TRsInRegion": "Other TRs Nearby",
}

df = df[list(rename_dict.keys())]
df.rename(columns=rename_dict, inplace=True)

print(df.columns)
print(df)
output_path = "table_of_properties_of_known_disease_associated_loci.tsv"
df.to_csv(output_path, sep="\t", index=False)

print(f"Wrote {len(df):,d} loci to {output_path}")
#%%