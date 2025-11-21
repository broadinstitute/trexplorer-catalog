import os
import pandas as pd

# source: https://github.com/ChaissonLab/vamos
version = "v2.1"

header = [
	"chrom",
	"start_1based",
	"end_1based",
	"motifs",
	"version",
	"STR_or_VNTR",
	"motif_size",
	"score1",
	"score2",
]

#url = f"https://zenodo.org/records/11625069/files/vamos.motif.hg38.{version}.e0.1.tsv.gz?download=1"
#header += ["segdup", "is_coding"]

url = f"https://zenodo.org/records/11625069/files/vamos.motif.hg38.{version}.orig.tsv.gz?download=1"
output_path_prefix = f"vamos_catalog.ori.{version}"

print(f"Downloading {url}")
df = pd.read_table(url, compression="gzip", names=header)
df["start_0based"] = df["start_1based"] - 1
df["locus_size"] = df["end_1based"] - df["start_0based"]
df["motif1"] = df["motifs"].str.split(",").str[0]
#df["motif_count"] = df["motifs"].str.count(",") + 1
#df["is_coding"] = df["is_coding"] == "coding"
#df["segdup"] = df["segdup"] == "segDup"


# save as TSV
df.to_csv(f"{output_path_prefix}.tsv.gz", sep="\t", index=False, header=True)
print(f"Wrote {len(df):,d} records to {output_path_prefix}.tsv.gz")

# save as BED
df = df[["chrom", "start_0based", "end_1based", "motif1", "motif_size"]]
df.to_csv(f"{output_path_prefix}.bed", sep="\t", index=False, header=False)
os.system(f"bgzip -f {output_path_prefix}.bed")
os.system(f"tabix -f {output_path_prefix}.bed.gz")
print(f"Wrote {len(df):,d} loci to {output_path_prefix}.bed.gz")

# save as TRGT catalog format
with open(f"{output_path_prefix}.TRGT.bed", "wt") as f:
	for _, row in df.iterrows():
		chrom = row["chrom"]
		start_0based = row["start_0based"]
		end_1based = row["end_1based"]
		motif = row["motif1"]
		locus_id = f"{chrom.replace('chr', '')}-{start_0based}-{end_1based}-{motif}"
		f.write("\t".join(map(str,
			[chrom, start_0based, end_1based, f"ID={locus_id};MOTIFS={motif};STRUC=({motif})n"]
		)) + "\n")
	# example: chr4  3074876  3074966  ID=HTT,MOTIFS=CAG,CCG;STRUC=(CAG)nCAACAG(CCG)n

os.system(f"bgzip -f {output_path_prefix}.TRGT.bed")
os.system(f"tabix -f {output_path_prefix}.TRGT.bed.gz")
print(f"Wrote {len(df):,d} records to {output_path_prefix}.TRGT.bed.gz")


notes = """Notes: 
len(df) == 1,175,953
df["version"] is always "v-2.0"
df["STR_or_VNTR"] is:
	STR     805485
	VNTR    370468
	
motif sizes:
	1      175625
	2      194556
	3	  	76515
	4      175497
	5       86683
	6       96609
	7       56180
	8        9294
	9	     5361
        ...
  491         1  
  492		  2
  494         1
  496         1
"""

def run(cmd):
	print(cmd)
	os.system(cmd)

# annotate and filter
cmd = f"python3 -m str_analysis.annotate_and_filter_str_catalog --verbose --show-progress --output-bed "
cmd += "--reference-fasta ~/hg38.fa "
cmd += "--min-motif-size 3 "
cmd += "--min-repeats-in-reference 1 "
cmd += "--skip-gene-annotations "
cmd += "--skip-mappability-annotations "
cmd += "--skip-disease-loci-annotations "
cmd += "--discard-loci-with-non-ACGT-bases-in-reference "
cmd += "--discard-loci-with-non-ACGTN-bases-in-motif "
cmd += "--set-locus-id "
cmd += f"{output_path_prefix}.bed.gz "
run(cmd)

# as a sanity check, see how often the motif matches the catalog
trexplorer_catalog_v1_path = f"results__2024-10-01/release_draft_2024-10-01/repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz"
print(f"Comparing to TRExplorer v1 catalog: {trexplorer_catalog_v1_path}")
cmd = f"python3  compare_with_loci_from_other_papers/compare_loci_with_catalog.py "
cmd += f"--catalog-bed-path TRExplorer_v1:{trexplorer_catalog_v1_path} "
#cmd += f"--write-loci-absent-from-new-catalog "
cmd += f"Vamos_v2_1:{output_path_prefix}.annotated_and_filtered.bed.gz"
run(cmd)

df = pd.read_table("vamos_catalog.ori.v2.1.overlap_with_TRExplorer_v1.tsv.gz")
sum(df["overlap_score"] == "same motif")
