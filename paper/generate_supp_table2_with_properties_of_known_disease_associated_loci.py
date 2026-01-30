"""This script generates Supp. Table 2 in the paper"""

import pandas as pd

KNOWN_DISEASE_ASSOCIATED_LOCUS_IDS = [
    "1-1435798-1435818-GGCGCGGAGC",
    "1-149390802-149390841-GGC",
    "1-57367043-57367118-AAAAT",
    "1-94418421-94418442-GCC",
    "10-79826383-79826404-GGC",
    "12-111598949-111599018-GCT",
    "12-123533720-123533750-GGC",
    "12-50505001-50505022-GGC",
    "12-6936716-6936773-CAG",
    "13-102161574-102161724-AAG",
    "13-70139383-70139428-CTG",
    "13-99985448-99985493-GCN",
    "14-23321472-23321490-GCG",
    "14-92071010-92071040-CTG",
    "15-22786677-22786701-GCG",
    "16-17470907-17470922-GCC",
    "16-24613439-24613529-TTTTA",
    "16-66490398-66490453-TAAAA",
    "16-67842863-67842950-CAG",
    "16-72787694-72787757-GCC",
    "16-87604287-87604329-CTG",
    "17-80147059-80147139-CCTCGCTGTGCCGCTGCCGA",
    "18-55586155-55586227-CAG",
    "19-13207858-13207897-CTG",
    "19-14496041-14496074-CCG",
    "19-18786034-18786049-GTC",
    "19-45770204-45770264-CAG",
    "2-176093058-176093103-GCG",
    "2-190880872-190880920-GCA",
    "2-96197066-96197121-AAAAT",
    "20-2652733-2652757-GGCCTG",
    "20-4699397-4699493-CCTCATGGTGGTGGCTGGGGGCAG",
    "21-43776443-43776479-CGCGGGGCGGGG",
    "22-19766762-19766807-GCN",
    "22-45795354-45795424-ATTCT",
    "3-129172576-129172656-CAGG",
    "3-138946020-138946062-NGC",
    "3-183712187-183712222-TTTTA",
    "3-63912684-63912714-GCA",
    "4-159342526-159342616-TTTTA",
    "4-3074876-3074933-CAG",
    "4-39348424-39348479-AAAAG",
    "4-41745972-41746032-GCN",
    "5-10356346-10356411-ATTTT",
    "5-146878727-146878757-GCT",
    "6-16327633-16327723-TGC",
    "6-170561906-170562017-GCA",
    "6-45422750-45422792-GCN",
    "7-27199678-27199732-NGC",
    "7-27199825-27199861-NGC",
    "7-27199924-27199966-NGC",
    "8-104588970-104588997-CGC",
    "8-118366812-118366918-AAAAT",
    "9-130681606-130681639-GCC",
    "9-27573528-27573546-GGCCCC",
    "9-69037286-69037304-GAA",
    "X-137566826-137566856-GCC",
    "X-140504316-140504361-NGC",
    "X-147912050-147912110-CGG",
    "X-148500631-148500691-GCC",
    "X-25013529-25013565-NGC",
    "X-25013649-25013697-NGC",
    "X-67545316-67545385-GCA",
]

df = pd.read_json("TR_catalog.63_loci.20251229_005446.json.gz")
df = df[df["LocusId"].isin(KNOWN_DISEASE_ASSOCIATED_LOCUS_IDS)]

missing_loci = set(KNOWN_DISEASE_ASSOCIATED_LOCUS_IDS) - set(df["LocusId"])
if missing_loci:
    raise ValueError(f"Missing {len(missing_loci)} loci from input: {missing_loci}")

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
  "AoU1027_OE_LengthPercentile": "Constraint (O/E Length Percentile)",
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