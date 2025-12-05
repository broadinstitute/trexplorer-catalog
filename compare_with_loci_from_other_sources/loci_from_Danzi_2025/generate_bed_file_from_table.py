import argparse
import collections
import intervaltree
import gzip
import os
import pandas as pd

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

os.chdir(os.path.dirname(__file__))


parser = argparse.ArgumentParser()
parser.add_argument("--adotto-catalog-bed-path", default="adotto_tr_catalog_v1.2.bed.gz")
args = parser.parse_args()

# Parse [Danzi et al. 2025] Supp. Table 3 and filter it to the subset of intervals that are outliers in OE_len or combinedLPSStdev
input_supp_table3_path = "media-3.xlsx"
print(f"Reading {input_supp_table3_path}")
df = pd.read_excel(input_supp_table3_path)
print(f"Read {len(df):,d} intervals from {input_supp_table3_path}")
df = df[["TRID", "combinedLPSStdev", "expectedCombinedLPSStdev", "OE_len", "CPS", "expectedCPS", "OE_motif"]]
df["combinedLPSStdev_percentile"] = df["combinedLPSStdev"].rank(pct=True)  # compute percentile of combinedLPSStdev
df[["chrom", "start_0based", "end_1based"]] = df["TRID"].str.split("_", expand=True)
df["start_0based"] = df["start_0based"].astype(int)
df["end_1based"] = df["end_1based"].astype(int)

# Thresholds are based on Figure 5C in [Danzi et al. 2025]
# https://www.biorxiv.org/content/10.1101/2025.01.06.631535v2.full
OE_len_threshold1 = 0.33
OE_len_threshold2 = 2
combinedLPSStdev_min_threshold = 0.99  # minimum percentile of combinedLPSStdev to keep
total = len(df)
df = df[(df["OE_len"] <= OE_len_threshold1) 
    | (df["OE_len"] >= OE_len_threshold2) 
    | (df["combinedLPSStdev_percentile"] >= combinedLPSStdev_min_threshold)].copy()

print(f"Kept {len(df):,d} out of {total:,d} ({(len(df) / total):.1%}) intervals after filtering by OE_len <= {OE_len_threshold1} or >= {OE_len_threshold2} or combinedLPSStdev_percentile >= {combinedLPSStdev_min_threshold}")
print(f"Of these, {sum(df['OE_len'] <= OE_len_threshold1):,d} have OE_len <= {OE_len_threshold1}, {sum(df['OE_len'] >= OE_len_threshold2):,d} have OE_len >= {OE_len_threshold2}, {sum(df['combinedLPSStdev_percentile'] >= combinedLPSStdev_min_threshold):,d} have combinedLPSStdev_percentile >= {combinedLPSStdev_min_threshold}")

# Parse additional table provided by Matt Danzi that contains the motifs
df2_path = "longestPureSegmentQuantiles.txt.gz"
df2 = pd.read_table(df2_path)
print(f"Read {len(df2):,d} rows from {df2_path}")

df2.sort_values(by=["N_motif"], inplace=True, ascending=False)
before = len(df2)
df2 = df2.groupby("TRID").first().reset_index()
print(f"Kept {len(df2):,d} out of {before:,d} ({(len(df2) / before):.1%}) rows with unique TRIDs in {df2_path}")
trid_to_motif_lookup = {
    trid: motif for trid, motif in zip(df2["TRID"], df2["longestPureSegmentMotif"])  # Example row: "chr10_100003328_100003396", "AAAC"
}

df["motif"] = df["TRID"].map(trid_to_motif_lookup)

output_bed_path = "Danzi_2025_OE_or_LPSStdev_outliers.wide_boundaries.bed"
with open(output_bed_path, "wt") as f:
    for _, row in df.iterrows():
        f.write(f"{row['chrom']}\t{row['start_0based']}\t{row['end_1based']}\t{row['motif']}\n")
os.system(f"bgzip -f {output_bed_path}")
os.system(f"tabix -f {output_bed_path}.gz")
print(f"Wrote {len(df):,d} intervals to {output_bed_path}.gz")


print(f"Added motifs to {sum(df['motif'].notna()):,d} out of {len(df):,d} ({sum(df['motif'].notna()) / len(df):.1%}) outlier intervals from {input_supp_table3_path}")
if sum(df['motif'].isna()) > 0:
    raise(f"ERROR: {df2_path} did not contain a motif for {len(sum(df['motif'].isna())):,d} TRIDs")

# Load Adotto narrow locus definitions
adotto_catalog_interval_trees = collections.defaultdict(intervaltree.IntervalTree)
with gzip.open(args.adotto_catalog_bed_path, "rt") as f:
    for line in f:
        fields = line.strip().split("\t")
        chrom = fields[0].replace("chr", "")
        start_0based = int(fields[1])
        end_1based = int(fields[2])
        motif = fields[3]
        canonical_motif = compute_canonical_motif(motif)
        adotto_catalog_interval_trees[chrom].add(intervaltree.Interval(start_0based, end_1based, data={
            "motif": motif,
            "canonical_motif": canonical_motif,
        }))

# Find loci in the Adotto catalog that overlap the outlier intervals and have the same motif length
adotto_catalog_loci_that_overlap_outlier_intervals = []
adotto_catalog_loci_that_overlap_outlier_intervals_and_have_same_length_motifs = []
adotto_catalog_loci_that_overlap_outlier_intervals_and_have_same_motifs = []
total_outlier_intervals_without_overlapping_adotto_loci = 0
for index, row in df.iterrows():
    chrom = row["chrom"].replace("chr", "")
    start_0based = row["start_0based"]
    end_1based = row["end_1based"]
    motif = row["motif"]
    canonical_motif = compute_canonical_motif(motif)
    overlapping_loci = adotto_catalog_interval_trees[chrom].overlap(start_0based, end_1based)
    if len(overlapping_loci) == 0:
        total_outlier_intervals_without_overlapping_adotto_loci += 1
    else:
        for overlapping_interval in overlapping_loci:
            adotto_catalog_loci_that_overlap_outlier_intervals.append((
                chrom,
                overlapping_interval.begin,
                overlapping_interval.end,
                overlapping_interval.data["motif"],
            ))

            if canonical_motif == overlapping_interval.data["canonical_motif"]:
                adotto_catalog_loci_that_overlap_outlier_intervals_and_have_same_motifs.append((
                    chrom,
                    overlapping_interval.begin,
                    overlapping_interval.end,
                    motif,
                ))

            if len(motif) == len(overlapping_interval.data["motif"]):
                adotto_catalog_loci_that_overlap_outlier_intervals_and_have_same_length_motifs.append((
                    chrom,
                    overlapping_interval.begin,
                    overlapping_interval.end,
                    overlapping_interval.data["motif"],
                ))


if total_outlier_intervals_without_overlapping_adotto_loci > 0:
    print(f"WARNING: Found {total_outlier_intervals_without_overlapping_adotto_loci:,d} out of {len(df):,d} ({total_outlier_intervals_without_overlapping_adotto_loci / len(df):.1%}) outlier intervals that did not overlap any locus in the Adotto catalog.")

print("--------------------------------")
print("Matching canonical motifs:")
print(f"Found {len(adotto_catalog_loci_that_overlap_outlier_intervals_and_have_same_motifs):,d} out of {len(df):,d} outlier intervals "
      f"overlapped a locus in the Adotto catalog that had the same motif "
      f"(average: {len(adotto_catalog_loci_that_overlap_outlier_intervals_and_have_same_motifs)/len(df):.2} matching loci per outlier interval).")
print(f"Found that the outlier intervals overlapped a total of {len(adotto_catalog_loci_that_overlap_outlier_intervals):,d} loci in the Adotto catalog, and "
      f"{len(adotto_catalog_loci_that_overlap_outlier_intervals_and_have_same_motifs):,d} of these had the same motif as the outlier interval.")

print("--------------------------------")
print("Matching motif lengths:")
print(f"Found {len(adotto_catalog_loci_that_overlap_outlier_intervals_and_have_same_length_motifs):,d} out of {len(df):,d} outlier intervals "
      f"overlapped a locus in the Adotto catalog that had a same-length motif "
      f"(average: {len(adotto_catalog_loci_that_overlap_outlier_intervals_and_have_same_length_motifs)/len(df):.2} matching loci per outlier interval).")
print(f"Found that the outlier intervals overlapped a total of {len(adotto_catalog_loci_that_overlap_outlier_intervals):,d} loci in the Adotto catalog, and "
      f"{len(adotto_catalog_loci_that_overlap_outlier_intervals_and_have_same_length_motifs):,d} of these had the same-length motif as the outlier interval.")

for locus_list, output_bed_path in [
    (adotto_catalog_loci_that_overlap_outlier_intervals, "Danzi_2025_OE_or_LPSStdev_outliers.narrow_adotto_boundaries.all_overlapping_loci.bed"),
    (adotto_catalog_loci_that_overlap_outlier_intervals_and_have_same_length_motifs, "Danzi_2025_OE_or_LPSStdev_outliers.narrow_adotto_boundaries.only_matching_motif_lengths.bed"),
    (adotto_catalog_loci_that_overlap_outlier_intervals_and_have_same_motifs, "Danzi_2025_OE_or_LPSStdev_outliers.narrow_adotto_boundaries.only_matching_motifs.bed"),
]:
    locus_list.sort(key=lambda x: (x[0], x[1], x[2]))
    with open(output_bed_path, "wt") as f:
        for chrom, start_0based, end_1based, motif in locus_list:
            f.write(f"{chrom}\t{start_0based}\t{end_1based}\t{motif}\n")
    os.system(f"bgzip -f {output_bed_path}")
    os.system(f"tabix -f {output_bed_path}.gz")
    print(f"Wrote {len(locus_list):,d} loci to {output_bed_path}.gz")