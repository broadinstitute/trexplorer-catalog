"""NOTE: deprecated. This script was only used for the v1.0 catalog.

This script takes a BED file of variation clusters and a JSON file of all tandem repeats and writes out a new BED file
with all input variation clusters as well as any tandem repeats from teh input catalog that aren't embedded in variation clusters
(ie. are isolated repeats).
"""


import argparse
import collections
import gzip
import os
import simplejson as json
import re
import tqdm

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator
from str_analysis.convert_expansion_hunter_catalog_to_trgt_catalog import convert_expansion_hunter_record_to_trgt_rows

MINIMUM_CHANGE_TO_BOUNDARIES_THRESHOLD = 6

def run(cmd):
	print(cmd)
	os.system(cmd)


def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
	parser.add_argument("--known-pathogenic-loci-json-path", required=True, help="Path of ExpansionHunter catalog "
						"containing known pathogenic loci. This is used to retrieve the original locus boundaries for "
						"these loci since their IDs don't contain these coordinates the way that IDs of other loci do.")
	parser.add_argument("-o", "--output-bed-path", help="Path of output BED file.")
	parser.add_argument("--verbose", action="store_true")
	parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
	parser.add_argument("input_variation_clusters_bed_path", help="Path of the input variation clusters BED file")
	parser.add_argument("input_repeat_catalog", help="Catalog of all tandem repeats in JSON or BED format")
	args = parser.parse_args()

	if ".bed" not in args.input_variation_clusters_bed_path:
		parser.error("--input-variation-clusters-bed-path must have a '.bed' suffix")

	if not args.output_bed_path:
		args.output_bed_path = re.sub(".bed(.b?gz)?$", "", args.input_variation_clusters_bed_path)
		args.output_bed_path += ".and_isolated_TRs.bed"
	elif args.output_bed_path.endswith(".bed.gz"):
		args.output_bed_path = args.output_bed_path[:-3]
	elif not args.output_bed_path.endswith(".bed"):
		parser.error("--output-bed-path must have a '.bed' suffix")

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

	# get locus ids of TRs in variation clusters
	locus_ids_in_variation_clusters = set()
	counter = collections.Counter()
	output_bed_file = open(args.output_bed_path, "wt")
	fopen = gzip.open if args.input_variation_clusters_bed_path.endswith("gz") else open
	with fopen(args.input_variation_clusters_bed_path, "rt") as f:
		if args.show_progress_bar:
			f = tqdm.tqdm(f, unit=" records", unit_scale=True)

		for line in f:
			counter['output_total'] += 1
			output_bed_file.write(line)

			counter["variation_clusters"] += 1
			fields = line.strip("\n").split("\t")
			#chrom = fields[0]
			#start_0based = int(fields[1])
			#end_1based = int(fields[2])
			info_fields = fields[3]

			info_fields_dict = {}
			for key_value in info_fields.split(";"):
				key_value = key_value.split("=")
				if len(key_value) != 2:
					print(f"WARNING: skipping invalid key-value pair '{key_value}' in line {fields}")
					continue
				key, value = key_value
				info_fields_dict[key] = value

			#region = f"{chrom.replace('chr', '')}:{start_0based}-{end_1based}"
			for locus_id in info_fields_dict["ID"].split(","):
				if locus_id in known_pathogenic_reference_regions_lookup or locus_id.count("-") == 3:
					if locus_id in locus_ids_in_variation_clusters:
						raise ValueError(f"locus_id '{locus_id}' occurs more than once")
					locus_ids_in_variation_clusters.add(locus_id)
				else:
					raise ValueError(f"Unexpected locus_id '{locus_id}'")

	print(f"Parsed {counter['variation_clusters']:,d} variation clusters from {args.input_variation_clusters_bed_path}")

	all_catalog_ids = set()
	for record_i, record in enumerate(get_variant_catalog_iterator(args.input_repeat_catalog, show_progress_bar=args.show_progress_bar)):
		counter["TRs_from_catalog"] += 1
		all_catalog_ids.add(record["LocusId"])

		if record["LocusId"] in locus_ids_in_variation_clusters:
			counter["TRs_in_variation_clusters"] += 1
			continue
		for output_row in convert_expansion_hunter_record_to_trgt_rows(record_i, record):
			output_bed_file.write("\t".join(map(str, output_row)) + "\n")
			counter['output_total'] += 1
			counter['isolated_TRs'] += 1

	print(f"Parsed {counter['TRs_from_catalog']:,d} TRs from {args.input_repeat_catalog}")
	output_bed_file.close()

	unexpected_locus_ids_in_variation_cluster_catalog = locus_ids_in_variation_clusters - all_catalog_ids
	if unexpected_locus_ids_in_variation_cluster_catalog:
		raise ValueError(f"{len(unexpected_locus_ids_in_variation_cluster_catalog)} locus IDs in the variation "
						 f"cluster catalog were not found in the input repeat catalog. Examples: "
						 f"{list(unexpected_locus_ids_in_variation_cluster_catalog)[:100]}")

	print(f"{counter['TRs_in_variation_clusters']:,d} out of {counter['TRs_from_catalog']:,d} "
		  f"({counter['TRs_in_variation_clusters']/counter['TRs_from_catalog']*100:.2f}%) TRs were in variation clusters")
	run(f"bedtools sort -i {args.output_bed_path} | bgzip > {args.output_bed_path}.sorted")
	run(f"mv {args.output_bed_path}.sorted {args.output_bed_path}.gz")
	os.remove(args.output_bed_path)

	print(f"Added {counter['isolated_TRs']:,d} isolated TRs to {args.output_bed_path}.gz")
	print(f"Wrote {counter['output_total']:,d} rows to {args.output_bed_path}.gz")

	
if __name__ == "__main__":
	main()

