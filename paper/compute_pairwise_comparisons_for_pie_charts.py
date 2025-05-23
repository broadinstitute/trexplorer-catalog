"""Hail Batch pipeline for computing pairwise overlap between catalogs"""

import collections
import json
import hail as hl
import hailtop.fs as hfs
import logging
import os
import pandas as pd
from pprint import pformat
import re
import tqdm

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/tandem-repeat-catalogs@sha256:97703ddfdaf6c61f73def482f9afd6ef5b56ec2c8618cadf4f42ca788ecb95d5"

REFERENCE_FASTA = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"

OUTPUT_BASE_DIR = "gs://bw2-delete-after-60-days/tandem-repeat-catalog/compare/"
DOWNLOAD_TO_DIR = "pairwise_comparisons"

CATALOGS = [
	("KnownDiseaseAssociatedLoci", "https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json"),
	("Illumina174kPolymorphicTRs", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/illumina_variant_catalog.sorted.bed.gz"),
	("UCSC_SimpleRepeatTrack", "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz"),
	("VamosCatalog_v2.1", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/vamos_catalog.v2.1.bed.gz"),
	("GangSTR_v17", "https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver17.bed.gz"),
	("HipSTR_Catalog", "https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/hg38.hipstr_reference.bed.gz"),
	("PolymorphicTRsInT2TAssemblies", "https://storage.googleapis.com/str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers_keeping_loci_that_have_overlapping_variants/combined/merged_expansion_hunter_catalog.78_samples.json.gz"),
	("Adotto_v1.2", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/adotto_tr_catalog_v1.2.bed.gz"),
	("PerfectRepeatsInReference", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp.bed.gz"),
	("TRExplorer_Catalog_v1", "https://github.com/broadinstitute/tandem-repeat-catalog/releases/download/v1.0/repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.json.gz"),
	("PopSTR_Catalog", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/popstr_catalog_v2.bed.gz"),
	("PlatinumTRs_v1.0", "https://zenodo.org/records/13178746/files/human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.bed.gz"),
	("Chiu_et_al", "https://zenodo.org/records/11522276/files/hg38.v1.bed.gz"),
	#("MukamelVNTRs", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/mukamel_VNTR_catalog.bed.gz"),
]

CATALOG_NAMES = {
	"KnownDiseaseAssociatedLoci": "Known disease-associated loci",
	"Illumina174kPolymorphicTRs": "Illumina 174k polymorphic TRs",
	"UCSC_SimpleRepeatTrack": "UCSC simple repeat track",
	"VamosCatalog_v2.1": "Vamos v2.1",
	"GangSTR_v17": "GangSTR v17",
	"HipSTR_Catalog": "HipSTR catalog",
	"Adotto_v1.2": "Adotto v1.2",
	"PopSTR_Catalog": "PopSTR catalog",
	"PlatinumTRs_v1.0": "PlatinumTRs v1.0",
	"Chiu_et_al": "Chiu et al",
	"PerfectRepeatsInReference": "All perfect repeats in hg38",
	"PolymorphicTRsInT2TAssemblies": "Polymorphic TRs in 78 T2T assemblies",
	"TRExplorer_Catalog_v1": "TRExplorer Catalog v1",
}

def trim_catalog(bp, catalog_label, catalog_url, machine_size=1):
	s = bp.new_step(
		name=f"Trim: {catalog_label}",
		step_number=1,
		arg_suffix="trim",
		image=DOCKER_IMAGE,
		storage="20Gi",
		cpu=machine_size,
		memory="highmem" if machine_size > 1 else "standard",
		output_dir=OUTPUT_BASE_DIR,
	)
	s.command("set -ex")

	local_reference_fasta = s.input(REFERENCE_FASTA)
	s.command(f"wget -q {catalog_url}")

	catalog_filename = os.path.basename(catalog_url)
	if catalog_label == "PlatinumTRs_v1.0":
		output_path = catalog_filename.replace(".bed.gz", "") + ".catalog.json.gz"

		s.command(f"python3 -u -m str_analysis.convert_trgt_catalog_to_expansion_hunter_catalog "
				  f"-r {local_reference_fasta} {catalog_filename} -o {output_path}")
		catalog_filename = output_path
	elif catalog_label == "GangSTR_v17":
		output_path = catalog_filename.replace(".bed.gz", "") + ".catalog.json.gz"
		s.command(f"python3 -u -m str_analysis.convert_gangstr_spec_to_expansion_hunter_catalog "
				  f"--verbose {catalog_filename} -o {output_path}")
		catalog_filename = output_path		
	elif catalog_label == "HipSTR_Catalog":
		output_path = catalog_filename.replace(".bed.gz", "") + ".catalog.bed.gz"
		s.command(f"python3 -u /scripts/convert_hipstr_catalog_to_regular_bed_file.py "
				  f"{catalog_filename} -o {output_path}")
		catalog_filename = output_path
	elif catalog_label == "UCSC_SimpleRepeatTrack":
		output_path = catalog_filename.replace(".txt.gz", "") + "_track_from_UCSC.bed.gz"
		s.command(f"python3 -u /scripts/convert_ucsc_simple_repeat_track_to_bed.py "
				  f"{catalog_filename} -o {output_path}")
		catalog_filename = output_path
	elif catalog_label == "Chiu_et_al":
		output_path = "hg38.Chiu_et_al.v1.bed.gz"
		s.command(f"gunzip -c {catalog_filename} | tail -n +3 | cut -f 1-4 | awk 'BEGIN {{OFS=\"\\t\"}} {{ print( $1, $2 - 1, $3, $4 ) }}' | gzip -c - > {output_path}")
		catalog_filename = output_path

	output_filename = f"trimmed.{catalog_label}.json.gz"
	s.command(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog \
		--reference-fasta {local_reference_fasta} \
		--skip-gene-annotations \
		--skip-disease-loci-annotations \
		--skip-mappability-annotations \
		--verbose \
		--trim-loci \
		--output-stats \
		--output-path {output_filename} \
		{catalog_filename}
	""")

	s.output(output_filename, download_to_dir=DOWNLOAD_TO_DIR)

	return s, os.path.join(OUTPUT_BASE_DIR, output_filename)


def compute_catalog_stats(bp, catalog_label, catalog_path, machine_size=1):
	s = bp.new_step(
		name=f"Stats: {catalog_label}",
		step_number=2,
		arg_suffix="stats",
		image=DOCKER_IMAGE,
		cpu=machine_size,
		memory="highmem" if machine_size > 1 else "standard",
		output_dir=OUTPUT_BASE_DIR,
	)
	s.command("set -ex")

	local_catalog_path = s.input(catalog_path)

	s.command(f"python3 -u -m str_analysis.compute_catalog_stats {local_catalog_path}")

	output_filename = re.sub("(.json|.bed)(.gz)?$", "", os.path.basename(catalog_path)) + ".catalog_stats.tsv"
	s.output(output_filename, download_to_dir=DOWNLOAD_TO_DIR)

	return s, os.path.join(OUTPUT_BASE_DIR, output_filename)


def process_pair(bp, catalog1, path1, catalog2, path2, machine_size=1):

	s = bp.new_step(
		name=f"Compare: {catalog1} vs. {catalog2}",
		step_number=3,
		arg_suffix="compare",
		image=DOCKER_IMAGE,
		storage="20Gi",
		cpu=machine_size,
		memory="highmem" if machine_size > 1 else "standard",
		output_dir=OUTPUT_BASE_DIR,
	)

	s.command("set -ex")

	catalog1_local_path = s.input(path1)
	if path1 == path2:
		catalog2_local_path = catalog1_local_path
	else:
		catalog2_local_path = s.input(path2)

	s.command(f"""python3 -u -m str_analysis.merge_loci \
		--discard-extra-fields-from-input-catalogs \
		--write-bed-files-with-unique-loci \
		--overlap-fraction 0.05 \
		--write-outer-join-table \
		--output-prefix merged___{catalog1}___vs___{catalog2} \
		{catalog1}:{catalog1_local_path} \
		{catalog2}:{catalog2_local_path}
	""")

	s.command(f"ls -ltrh")

	output_filename = f"merged___{catalog1}___vs___{catalog2}.outer_join_overlap_table.tsv.gz"
	unique_loci_output_filename = f"merged___{catalog1}___vs___{catalog2}.trimmed.{catalog2}.unique_loci.bed.gz"
	s.output(output_filename, download_to_dir=DOWNLOAD_TO_DIR)
	s.output(unique_loci_output_filename, download_to_dir=DOWNLOAD_TO_DIR)

	return s, os.path.join(OUTPUT_BASE_DIR, output_filename)


def main():
	bp = pipeline("pairwise catalog comparison", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

	parser = bp.get_config_arg_parser()
	parser.add_argument("--by-motif-size", action="store_true", help="Stratify the results by motif size")
	args, _ = parser.parse_known_args()

	# trim catalogs
	step1_map = {}
	step1_output_paths = {}
	local_catalog_stats_tables = {}
	for i, (catalog_label, catalog_url) in enumerate(CATALOGS):
		step1, step1_output_paths[catalog_label] = trim_catalog(
			bp, catalog_label, catalog_url, machine_size=8 if i >= 5 else 4 if i >= 2 else 1)
		step1_map[catalog_label] = step1

		step2, _ = compute_catalog_stats(
			bp, catalog_label, step1_output_paths[catalog_label], machine_size=2 if i >= 2 else 1)
		step2.depends_on(step1)

		local_output_paths = step2.get_output_paths_to_download_when_done()
		assert len(local_output_paths) == 1
		local_catalog_stats_tables[catalog_label] = local_output_paths[0][1]

	# do pair-wise comparisons
	local_join_table_paths = {}
	for i in range(len(CATALOGS)):
		catalog1, _ = CATALOGS[i]
		for j in range(len(CATALOGS)):
			if i == j:
				continue
			catalog2, _ = CATALOGS[j]

			print(f"Processing {catalog1} vs {catalog2}")
			step2, output_path = process_pair(
				bp, catalog1, step1_output_paths[catalog1], catalog2, step1_output_paths[catalog2],
				machine_size=16 if i >= 7 or j >= 7 else 8 if i >= 5 or j >= 5 else 4 if i >= 3 or j >= 3 else 1)
			step2.depends_on(step1_map[catalog1])
			step2.depends_on(step1_map[catalog2])

			local_output_paths = step2.get_output_paths_to_download_when_done()
			assert len(local_output_paths) > 1
			local_join_table_paths[(catalog1, catalog2)] = local_output_paths[0][1]

	bp.run()

	output_table_path = "pairwise_catalog_comparison_results.tsv"
	with open(output_table_path, "wt") as f:
		f.write("\t".join([
			"catalog_size",
			"catalog_size_by_motif_size",
			"motif_size",
			"motif_label",
			"category",
			"catalog1",
			"catalog2",
			"count"]) + "\n")

		catalog_size = {}
		catalog_size_by_motif_size = {}
		for i, (catalog_label, local_catalog_stats_table) in enumerate(local_catalog_stats_tables.items()):
			df = pd.read_table(local_catalog_stats_table)
			catalog_size[catalog_label] = df.iloc[0].total

			if catalog_label != "KnownDiseaseAssociatedLoci":
				filename = f"merged___{catalog_label}___vs___KnownDiseaseAssociatedLoci.outer_join_overlap_table.tsv.gz"
			else:
				filename = f"merged___{catalog_label}___vs___Illumina174kPolymorphicTRs.outer_join_overlap_table.tsv.gz"
			path = os.path.join("pairwise_comparisons", filename)
			local_join_table_paths[(catalog_label, catalog_label)] = path
			df_full = pd.read_table(path)
			df_full["motif_size"] = df_full["LocusId"].apply(lambda x: len(x.split("-")[3]))
			df_full = df_full[df_full[catalog_label] == "Yes"]

			for min_motif_size, max_motif_size in [
				(1,1), (2,2), (3,3), (4,4), (5,5), (6,6), (7, 24), (25, None),
			]:
				filter_condition = (df_full["motif_size"] >= min_motif_size)
				if max_motif_size is not None:
					filter_condition &= (df_full["motif_size"] <= max_motif_size)
				df = df_full[filter_condition]
				catalog_size_by_motif_size[(catalog_label, min_motif_size)] = len(df)

		for (catalog1, catalog2), local_outer_join_tsv_path in tqdm.tqdm(local_join_table_paths.items(), unit=" comparison"):
			print(f"Comparison of {catalog1} vs {catalog2}: {local_outer_join_tsv_path}")

			df_full = pd.read_table(local_outer_join_tsv_path)
			df_full["motif_size"] = df_full["LocusId"].apply(lambda x: len(x.split("-")[3]))
			if catalog1 == catalog2:
				df_full = df_full[df_full[catalog1] == "Yes"]
				df_full = df_full[["LocusId", catalog1, "motif_size"]]

			if args.by_motif_size:
				motif_size_ranges = [
					(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (7, 24), (25, None),
				]
			else:
				motif_size_ranges = [(None, None)]

			for min_motif_size, max_motif_size in motif_size_ranges:
				if min_motif_size == max_motif_size:
					motif_size_label = min_motif_size
				elif max_motif_size is None:
					motif_size_label = f"{min_motif_size}+"
				else:
					motif_size_label = f"{min_motif_size}-{max_motif_size}"

				filter_condition = None
				if min_motif_size is not None:
					filter_condition = (df_full["motif_size"] >= min_motif_size)
				if max_motif_size is not None:
					filter_condition &= (df_full["motif_size"] <= max_motif_size)

				df = df_full[filter_condition] if filter_condition is not None else df_full

				for count_type in "intersection", "unique1", "unique2", "widerInCatalog1", "widerInCatalog2":
					if count_type == "intersection":
						count = sum(df[catalog1].isin({"Yes", "YesButShifted"}) & df[catalog2].isin({"Yes", "YesButShifted"}))
					elif count_type == "unique1":
						count = sum(~df[catalog1].isna() & df[catalog2].isna())
					elif count_type == "unique2":
						count = sum(df[catalog1].isna() & ~df[catalog2].isna())
					elif count_type == "widerInCatalog1":
						count = sum((df[catalog1] == "YesButWider") | (df[catalog2] == "YesButNarrower"))
					elif count_type == "widerInCatalog2":
						count = sum((df[catalog2] == "YesButWider") | (df[catalog1] == "YesButNarrower"))
					else:
						raise ValueError(f"Unexepcted count_type: {count_type}")

					f.write("\t".join(map(str, [
						catalog_size[catalog1],
						catalog_size_by_motif_size.get((catalog1, min_motif_size)),
						min_motif_size,
						motif_size_label,
						count_type,
						CATALOG_NAMES[catalog1],
						CATALOG_NAMES[catalog2],
						count,
					])) + "\n")

				assert all(df[catalog1].isin({"Yes", "YesButShifted", "YesButNarrower", "YesButWider", float("nan")}))
				assert sum(df[catalog1].isna() & df[catalog2].isna()) == 0
				assert sum(df[catalog1].isna() & (df[catalog2] != "Yes")) == 0
				assert sum(df[catalog2].isna() & (df[catalog1] != "Yes")) == 0

	print(f"Wrote results to {output_table_path}")


if __name__ == "__main__":
	main()
