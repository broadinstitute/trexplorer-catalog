"""Add longest pure segment (LPS) annotations to catalog:

LPSLengthStdevFromHPRC100 is the standard deviation of the length of the longest pure segment (LPS) detected at
	a tandem repeat locus in 100 high-coverage long-read (LR) samples from the HPRC  (Range: 0 to 1226.0)
LPSMotifFractionFromHPRC100 is the fraction of alleles at a tandem repeat locus in which the given motif
	composed the longest pure segment (LPS) among 100 high-coverage long-read (LR) samples from the HPRC.
"""

import argparse
import collections
import gzip
import ijson
import os
import pandas as pd
import simplejson as json
import sys
import tqdm

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure
"""
Expected columns in lps table:

$1                       TRID : 1-49834-49853-AAAC
$2    longestPureSegmentMotif : AAAC
$3                    N_motif : 2018.0
$4              0thPercentile : 4.0
$5              1stPercentile : 8.0
$6              5thPercentile : 12.0
$7             10thPercentile : 12.0
$8             15thPercentile : 12.0
$9             20thPercentile : 12.0
$10            25thPercentile : 12.0
$11            30thPercentile : 12.0
$12            35thPercentile : 16.0
$13            40thPercentile : 16.0
$14            45thPercentile : 16.0
$15            50thPercentile : 16.0
$16            55thPercentile : 16.0
$17            60thPercentile : 16.0
$18            65thPercentile : 16.0
$19            70thPercentile : 16.0
$20            75thPercentile : 16.0
$21            80thPercentile : 16.0
$22            85thPercentile : 16.0
$23            90thPercentile : 16.0
$24            95thPercentile : 16.0
$25            99thPercentile : 16.0
$26          99.9thPercentile : 16.0
$27           100thPercentile : 16.0
$28                       MAD : 0.0
$29                      Mean : 14.545094152626362
$30                     Stdev : 0.5168253213498479
$31                      Mode : 4.0
$32                     TRID2 : 1-49834-49853-AAAC
$33                     motif : AAAC
$34           canonical_motif : AAAC
$35                numAlleles : 27
$36          numCalledAlleles : 2020
$37          combinedLPSStdev : 2.068555633665116
$38  expectedCombinedLPSStdev : 1.226411599717534
$39                    OE_len : 1.6349040027446389
$40         OE_len_percentile : 0.9381970214187625
"""

def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--known-pathogenic-loci-json-path", required=True, help="Path of ExpansionHunter catalog "
						"containing known pathogenic loci. This is used to retrieve the original locus boundaries for "
						"these loci since their IDs don't contain these coordinates the way that IDs of other loci do.")
	parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
	parser.add_argument("--output-catalog-json-path",
						help="Path of the output catalog JSON file that includes variation cluster annotations")
	parser.add_argument("aou1027_table", help="Path of the LPS data table", default="HPRC_100_LongestPureSegmentQuantiles.txt.gz")
	parser.add_argument("catalog_json_path", help="Path of the JSON catalog to annotate")
	args = parser.parse_args()

	for path in args.aou1027_table, args.catalog_json_path, args.known_pathogenic_loci_json_path:
		if not os.path.isfile(path):
			parser.error(f"{path} file not found")

	if not args.output_catalog_json_path:
		args.output_catalog_json_path = args.catalog_json_path.replace(".json", ".with_LPS_annotations.json")

	print(f"Parsing {args.known_pathogenic_loci_json_path}")
	fopen = gzip.open if args.known_pathogenic_loci_json_path.endswith("gz") else open
	with fopen(args.known_pathogenic_loci_json_path, "rt") as f:
		known_pathogenic_loci = json.load(f)
		known_pathogenic_reference_regions_lookup = {}
		for locus in known_pathogenic_loci:
			motifs = parse_motifs_from_locus_structure(locus["LocusStructure"])
			if isinstance(locus["ReferenceRegion"], list):
				assert isinstance(locus["VariantId"], list)
				assert len(locus["ReferenceRegion"]) == len(locus["VariantId"])
				assert len(locus["ReferenceRegion"]) == len(motifs)
				for variant_id, reference_region, motif in zip(locus["VariantId"], locus["ReferenceRegion"], motifs):
					known_pathogenic_reference_regions_lookup[variant_id] = (reference_region, motif)
			else:
				known_pathogenic_reference_regions_lookup[locus["LocusId"]] = (locus["ReferenceRegion"], motifs[0])

	print(f"Parsed {len(known_pathogenic_reference_regions_lookup)} known pathogenic loci")

	print(f"Parsing {args.aou1027_table}")
	df = pd.read_table(args.aou1027_table)
	missing_columns = {"TRID", "longestPureSegmentMotif", "N_motif", "Stdev"} - set(df.columns)
	if missing_columns:
		parser.error(f"{args.aou1027_table} is missing expected columns: {missing_columns}")

	before = len(df)
	df = df[~df["longestPureSegmentMotif"].isna() & ~df["Stdev"].isna() & ~df["N_motif"].isna()]
	print(f"Filtered out {before - len(df):,d} out of {before:,d} ({(before - len(df)) / before:.1%}) records with missing values")

	annotation_lookup = {}
	row_iterator = df.iterrows()
	if args.show_progress_bar:
		row_iterator = tqdm.tqdm(row_iterator, total=len(df), unit=" records", unit_scale=True)
		
	for _, row in row_iterator:
		locus_id = row["TRID2"]
		motif_size = len(row['longestPureSegmentMotif'])
		# convert stdev in bp to stdev in repeat units
		annotation_lookup[locus_id] = {
			"StdevFromAoU1027": round(row["Stdev"] / motif_size, 3),
			"MedianAlleleFromAoU1027": row["50thPercentile"] // motif_size,
			"MaxAlleleFromAoU1027": row["100thPercentile"] // motif_size,
			"NumUniqueAllelesFromAoU1027": int(row["numAlleles"]),
		}

	input_locus_counter = annotated_locus_counter = 0
	print(f"Adding AoU1027 annotations to {args.catalog_json_path}")
	fopen = gzip.open if args.catalog_json_path.endswith("gz") else open
	with fopen(args.catalog_json_path, "rt") as f:
		f2open = gzip.open if args.output_catalog_json_path.endswith("gz") else open
		with f2open(args.output_catalog_json_path, "wt") as f2:
			iterator = ijson.items(f, "item")
			if args.show_progress_bar:
				iterator = tqdm.tqdm(iterator, unit=" records", unit_scale=True)
			f2.write("[")
			for i, record in enumerate(iterator):
				locus_id = record["LocusId"]
				input_locus_counter += 1
				if locus_id in annotation_lookup:
					record.update(annotation_lookup[locus_id])
					annotated_locus_counter += 1
				if i > 0:
					f2.write(", ")
				f2.write(json.dumps(record, f2, use_decimal=True, indent=4))
			f2.write("]")

	print(f"Annotated {annotated_locus_counter:,d} out of {input_locus_counter:,d} loci")
	print(f"Wrote annotated catalog to {args.output_catalog_json_path}")

if __name__ == "__main__":
	main()

