"""Add variation cluster annotations to catalog"""

import argparse
import collections
import gzip
import ijson
import os
import simplejson as json
import sys
import tqdm

from str_analysis.utils.misc_utils import parse_interval

MINIMUM_CHANGE_TO_BOUNDARIES_THRESHOLD = 6

def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--known-pathogenic-loci-json-path", required=True, help="Path of ExpansionHunter catalog "
						"containing known pathogenic loci. This is used to retrieve the original locus boundaries for "
						"these loci since their IDs don't contain these coordinates the way that IDs of other loci do.")
	parser.add_argument("--verbose", action="store_true")
	parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
	parser.add_argument("--output-catalog-json-path", help="Path of the output catalog JSON file that includes variation cluster annotations")
	parser.add_argument("--generate-plot", action="store_true", help="Generate a plot of the size differences between variation clusters and simple repeats")
	parser.add_argument("variation_clusters_bed_path", help="Path of the variation clusters BED file")
	parser.add_argument("catalog_json_path", help="Path of the JSON catalog to annotate")
	args = parser.parse_args()

	for path in args.variation_clusters_bed_path, args.catalog_json_path, args.known_pathogenic_loci_json_path:
		if not os.path.isfile(path):
			parser.error(f"File not found: {path}")

	if not args.output_catalog_json_path:
		args.output_catalog_json_path = args.catalog_json_path.replace(".json", ".with_variation_clusters.json")

	print(f"Parsing {args.known_pathogenic_loci_json_path}")
	fopen = gzip.open if args.known_pathogenic_loci_json_path.endswith("gz") else open
	with fopen(args.known_pathogenic_loci_json_path, "rt") as f:
		known_pathogenic_loci = json.load(f)
		known_pathogenic_reference_regions_lookup = {}
		for locus in known_pathogenic_loci:
			if isinstance(locus["ReferenceRegion"], list):
				assert isinstance(locus["VariantId"], list)
				assert len(locus["ReferenceRegion"]) == len(locus["VariantId"])
				for variant_id, reference_region in zip(locus["VariantId"], locus["ReferenceRegion"]):
					known_pathogenic_reference_regions_lookup[variant_id] = reference_region
			else:
				known_pathogenic_reference_regions_lookup[locus["LocusId"]] = locus["ReferenceRegion"]

	locus_id_to_variation_cluster_interval = {}
	locus_id_to_variation_cluster_size_difference_from_simple_repeat_boundaries = {}
	size_diff_histogram = collections.Counter()
	input_variation_clusters_counter = 0
	input_locus_ids_counter = 0
	almost_no_change_to_boundaries = 0
	output_variation_clusters_counter = 0
	examples = set()
	if args.verbose:
		print(f"Parsing {args.variation_clusters_bed_path}")

	fopen = gzip.open if args.variation_clusters_bed_path.endswith("gz") else open
	with fopen(args.variation_clusters_bed_path, "rt") as f:
		if args.show_progress_bar:
			f = tqdm.tqdm(f, unit=" records", unit_scale=True)

		for line in f:
			input_variation_clusters_counter += 1
			fields = line.strip("\n").split("\t")
			chrom = fields[0]
			start_0based = int(fields[1])
			end_1based = int(fields[2])
			info_fields = fields[3]

			info_fields_dict = {}
			for key_value in info_fields.split(";"):
				key_value = key_value.split("=")
				if len(key_value) != 2:
					print(f"WARNING: skipping invalid key-value pair '{key_value}' in line {fields}")
					continue
				key, value = key_value
				info_fields_dict[key] = value

			if "VC" not in info_fields_dict["STRUC"]:
				continue

			variation_cluster_differs_from_simple_repeat = False
			region = f"{chrom.replace('chr', '')}:{start_0based}-{end_1based}"
			for locus_id in info_fields_dict["ID"].split(","):
				input_locus_ids_counter += 1
				if locus_id in known_pathogenic_reference_regions_lookup:
					region2 = known_pathogenic_reference_regions_lookup[locus_id]
					original_chrom, original_start_0based, original_end_1based = parse_interval(region2)
				elif locus_id.count("-") == 3:
					original_chrom, original_start_0based, original_end_1based, _ = locus_id.split("-")
					region2 = f"{original_chrom}:{original_start_0based}-{original_end_1based}"
				else:
					raise ValueError(f"Unexpected locus_id '{locus_id}'")

				original_start_0based = int(original_start_0based)
				original_end_1based = int(original_end_1based)

				if abs(end_1based - original_end_1based) < MINIMUM_CHANGE_TO_BOUNDARIES_THRESHOLD and abs(original_start_0based - start_0based) < MINIMUM_CHANGE_TO_BOUNDARIES_THRESHOLD:
					almost_no_change_to_boundaries += 1
					if len(examples) < 5:
						examples.add(f"VC:{region} and locus:{locus_id}")
					print(f"{region} doesn't change {locus_id} by {MINIMUM_CHANGE_TO_BOUNDARIES_THRESHOLD}bp or more")
				else:
					size_diff = abs(end_1based - original_end_1based) + abs(original_start_0based - start_0based)
					variation_cluster_differs_from_simple_repeat = True
					locus_id_to_variation_cluster_interval[locus_id] = region
					locus_id_to_variation_cluster_size_difference_from_simple_repeat_boundaries[locus_id] = size_diff
					size_diff_histogram[size_diff] += 1

			if variation_cluster_differs_from_simple_repeat:
				output_variation_clusters_counter += 1

	locus_ids_in_variation_cluster_above_threshold = len(locus_id_to_variation_cluster_interval)
	if args.verbose:
		print(f"Parsed {input_variation_clusters_counter:,d} variation clusters that contained {input_locus_ids_counter:,d} simple TR ids")
		if almost_no_change_to_boundaries:
			print(f"Found {almost_no_change_to_boundaries:,d} out of {input_variation_clusters_counter:,d} "
				  f"({almost_no_change_to_boundaries/input_variation_clusters_counter:.1%}) "
				  f"variation clusters that did not change the original locus boundaries "
				  f"by {MINIMUM_CHANGE_TO_BOUNDARIES_THRESHOLD} bases or more")
			print(f"These contained {input_locus_ids_counter - locus_ids_in_variation_cluster_above_threshold:,d} out of {input_locus_ids_counter:,d} "
				  f"({(input_locus_ids_counter - locus_ids_in_variation_cluster_above_threshold)/input_locus_ids_counter:.1%}) locus IDs. "
				  f"Examples: ", ", ".join(examples))
		print(f"Found {output_variation_clusters_counter:,d} out of {input_variation_clusters_counter:,d} "
			  f"({output_variation_clusters_counter/input_variation_clusters_counter:.1%}) variation clusters "
			  f"differed from simple TRs by at least {MINIMUM_CHANGE_TO_BOUNDARIES_THRESHOLD}bp")

	print(f"Annotating {args.catalog_json_path} with variation cluster annotations")
	fopen = gzip.open if args.catalog_json_path.endswith("gz") else open
	with fopen(args.catalog_json_path, "rt") as f:
		f2open = gzip.open if args.output_catalog_json_path.endswith("gz") else open
		with f2open(args.output_catalog_json_path, "wt") as f2:
			input_locus_counter = locus_without_variation_cluster_counter = output_locus_counter = 0
			iterator = ijson.items(f, "item")
			if args.show_progress_bar:
				iterator = tqdm.tqdm(iterator, unit=" records", unit_scale=True)
			f2.write("[")
			for i, record in enumerate(iterator):
				locus_id = record["LocusId"]
				input_locus_counter += 1
				if locus_id in locus_id_to_variation_cluster_interval:
					#record["VariationClusterId"] = locus_id_to_variation_cluster_id[locus_id]
					record["VariationCluster"] = locus_id_to_variation_cluster_interval[locus_id]
					record["VariationClusterSizeDiff"] = locus_id_to_variation_cluster_size_difference_from_simple_repeat_boundaries[locus_id]
					#del locus_id_to_variation_cluster_id[locus_id]
					del locus_id_to_variation_cluster_interval[locus_id]
				else:
					locus_without_variation_cluster_counter += 1
					#if args.verbose:
					#	print(f"WARNING: locus_id {locus_id} not found in variation cluster catalog")

				output_locus_counter += 1
				if i > 0:
					f2.write(", ")
				f2.write(json.dumps(record, f2, use_decimal=True, indent=4))
			f2.write("]")

	if args.generate_plot:
		print(f"{locus_without_variation_cluster_counter:,d} out of {input_locus_counter:,d} "
			  f"({locus_without_variation_cluster_counter/input_locus_counter:.1%}) loci from {args.catalog_json_path} are not in variation clusters")
		print(f"Wrote {output_locus_counter:,d} out of {input_locus_counter:,d} ({output_locus_counter/input_locus_counter:.1%}) records "
			  f"to {args.output_catalog_json_path}")

		# plot using seaborn
		print(f"Generating VC size diff histograms")
		import seaborn as sns
		import matplotlib.pyplot as plt
		plt.figure(figsize=(12, 6))
		sns.barplot(x=list(size_diff_histogram.keys()), y=list(size_diff_histogram.values()))
		plt.xlabel("Size difference")
		plt.ylabel("Count")
		plt.title("Size difference between variation clusters and original loci")
		output_prefix = args.output_catalog_json_path.replace(".json", "").replace(".gz", "") + ".size_diff_histogram"
		plt.savefig(f"{output_prefix}.png")
		# also create a plot with log y-scale
		plt.yscale("log")
		plt.savefig(f"{output_prefix}.log.png")
		print(f"Wrote VC size diff histograms to {output_prefix}.png and {output_prefix}.log.png")

if __name__ == "__main__":
	main()

