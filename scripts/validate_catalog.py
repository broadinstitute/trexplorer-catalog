"""Run basic sanity checks on the final catalogs prior to release"""

import argparse
import gzip
import ijson
import json


from str_analysis.utils.file_utils import open_file
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator

EXPECTED_KEYS_IN_ANNOTATED_CATALOG = {
	"ReferenceRegion": str,
	"LocusStructure": str,
	"VariantType": str,
	"LocusId": str,
	"Source": str,
	"CanonicalMotif": str,
	"GencodeGeneRegion": str,
	"RefseqGeneRegion": str,
	"ManeGeneRegion": str,
	#"GencodeGeneName": str,    <== these are optional,
	#"GencodeGeneId": str,
	#"GencodeTranscriptId": str,
	"NumRepeatsInReference": int,
	"NsInFlanks": int,
	"ReferenceRepeatPurity": float,
	"LeftFlankMappability": float,
	"FlanksAndLocusMappability": float,
	"RightFlankMappability": float,
}

# Define acceptable types for each expected type (allows for type compatibility)
TYPE_COMPATIBILITY = {
	float: (int, float),  # int is a narrower type than float
}

PATH_OF_CATALOG_V1_JSON = "https://github.com/broadinstitute/tandem-repeat-catalog/releases/download/v1.0/repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.json.gz"


def parse_known_pathogenic_loci(known_pathogenic_loci_json_path):
	catalog = {}
	fopen = gzip.open if known_pathogenic_loci_json_path.endswith("gz") else open
	with fopen(known_pathogenic_loci_json_path, "rt") as f:
		known_pathogenic_loci = json.load(f)
		for record in known_pathogenic_loci:
			catalog[record["LocusId"]] = record

	return catalog

def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--verbose", action="store_true")
	parser.add_argument("--check-for-presence-of-annotations", action="store_true", help="Check for the presence of "
						"annotations in the catalog")
	parser.add_argument("--check-for-presence-of-all-known-loci", action="store_true", help="Check that all known disease "
						"loci are present in the catalog")
	parser.add_argument("--check-for-presence-of-all-loci-from-v1", action="store_true", help="Check that all loci from the v1 "
						"catalog are also present in the latest catalog")
	parser.add_argument("--known-pathogenic-loci-json-path", required=True, help="Path of ExpansionHunter catalog "
						"containing known pathogenic loci. This is used to retrieve the original locus boundaries for "
						"these loci since their IDs don't contain these coordinates the way that IDs of other loci do.")
	parser.add_argument("simple_repeat_catalog_path")
	args = parser.parse_args()

	known_disease_loci_catalog = parse_known_pathogenic_loci(args.known_pathogenic_loci_json_path)

	known_locus_id_to_reference_region = {}
	for locus_id, record in known_disease_loci_catalog.items():
		if isinstance(record["ReferenceRegion"], list):
			for reference_region, variant_id in zip(record["ReferenceRegion"], record["VariantId"]):
				known_locus_id_to_reference_region[variant_id] = reference_region
		else:
			known_locus_id_to_reference_region[locus_id] = record["ReferenceRegion"]

	error_counter = 0


	# confirm that the locus ids and ReferenceRegions are present in the catalog
	locus_ids = set()
	reference_regions = set()

	v1_catalog_locus_ids = set()
	if args.check_for_presence_of_all_loci_from_v1:
		v1_catalog = open_file(PATH_OF_CATALOG_V1_JSON)
		for record in ijson.items(v1_catalog, "item"):
			v1_catalog_locus_ids.add(record["LocusId"])

	current_catalog_locus_ids = set()
	input_file_iterator = get_variant_catalog_iterator(args.simple_repeat_catalog_path)
	for i, record in enumerate(input_file_iterator):
		if not record["ReferenceRegion"].startswith("chr"):
			print(f"ERROR: ReferenceRegion {record['ReferenceRegion']} does not start with 'chr'")
			error_counter += 1

		if record["ReferenceRegion"].replace(":", "-").replace("chr", "") not in record["LocusId"] and record["LocusId"] not in known_locus_id_to_reference_region:
			print(f"ERROR: ReferenceRegion {record['ReferenceRegion']} does not match LocusId {record['LocusId']}")
			error_counter += 1

		if args.check_for_presence_of_all_loci_from_v1:
			current_catalog_locus_ids.add(record["LocusId"])

		locus_ids.add(record["LocusId"])
		reference_regions.add(record["ReferenceRegion"])

		if args.check_for_presence_of_annotations:
			for key, value_type in EXPECTED_KEYS_IN_ANNOTATED_CATALOG.items():

				if key not in record or record[key] is None:
					print(f"ERROR: Missing key {key} in record {i}: {record}")
					error_counter += 1
				else:
					# Get acceptable types for this value_type, defaulting to just the value_type itself
					acceptable_types = TYPE_COMPATIBILITY.get(value_type, (value_type,))
					if not isinstance(record[key], acceptable_types):
						print(f"ERROR: Expected {key} to be of type {value_type} but got {type(record[key])}")
						error_counter += 1

			for prefix in "Gencode", "Refseq", "Mane":
				if record[f"{prefix}GeneRegion"] != "intergenic":
					for key in [f"{prefix}GeneId", f"{prefix}TranscriptId"]:  # f"{prefix}GeneName",
						if key not in record or record[key] is None:
							print(f"ERROR: Missing {key} in record {i}: {record}")
							error_counter += 1

			#for allele_frequency_key, found_in_key in [
			#	("AlleleFrequenciesFromIllumina174k", "FoundInIllumina174kPolymorphicTRs"),
			#	("AlleleFrequenciesFromT2TAssemblies", "FoundInPolymorphicTRsInT2TAssemblies"),
			#]:
			#	if allele_frequency_key in record or record.get(found_in_key) == "Yes":
			#		if not record.get(allele_frequency_key):
			#			print(f"ERROR: {allele_frequency_key} is missing in record {i}: {record}")
			#			error_counter += 1
			#		if not record.get(found_in_key):
			#			print(f"ERROR: {found_in_key} is missing in record {i}: {record}")
			#			error_counter += 1


	# check that all known disease loci are in the simple repeat catalog
	if args.check_for_presence_of_all_known_loci:
		for locus_id, reference_region in known_locus_id_to_reference_region.items():
			if reference_region not in reference_regions:
				print(f"ERROR: {locus_id} ReferenceRegion {reference_region} not found in the simple repeat catalog")
				error_counter += 1

	if args.check_for_presence_of_all_loci_from_v1:
		missing_loci = v1_catalog_locus_ids - current_catalog_locus_ids
		if missing_loci:
			print(f"ERROR: {len(missing_loci):,d} v1 catalog LocusIds are missing from the current catalog: ",
				  ', '.join(list(sorted(missing_loci))[:10]), "...")
			error_counter += len(missing_loci)


	if error_counter > 0:
		raise ValueError(f"Found {error_counter} errors in {args.simple_repeat_catalog_path}")

	print(f"Done. Catalog passed validation.")

if __name__ == "__main__":
	main()