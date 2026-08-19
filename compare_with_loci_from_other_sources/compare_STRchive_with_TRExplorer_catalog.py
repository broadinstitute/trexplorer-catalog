# This script checks which of the disease-associated TR loci in STRchive (https://harrietdashnow.com/STRchive)
# are present in or absent from the TRExplorer v2 genome-wide catalog.
#
# The TRExplorer catalog contains millions of loci, so rather than reading all of it, this script queries a local
# tabix-indexed copy of the catalog BED file for only the intervals that overlap each STRchive locus.

import argparse
import json
from pathlib import Path
import pysam
import requests

from str_analysis.utils.misc_utils import reverse_complement

STRCHIVE_LOCI_URL = "https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/STRchive-loci.json"

DEFAULT_TREXPLORER_CATALOG_PATH = ("~/code/tandem-repeat-catalogs/results__2026-02-01/release_draft_2026-02-01/"
								  "TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz")

# motifs longer than this are considered VNTRs, for which only the motif lengths need to match
MAX_STR_MOTIF_SIZE = 6


def get_strchive_loci(strchive_loci_path=None):
	"""Load the STRchive locus metadata, either from a local json path or from the STRchive github repo.

	Args:
		strchive_loci_path (str): Optional local path of the STRchive-loci.json file.

	Return:
		list: dicts with the fields needed for the comparison, sorted by locus id.
	"""
	if strchive_loci_path:
		print(f"Loading STRchive loci from {strchive_loci_path}")
		with open(Path(strchive_loci_path).expanduser()) as f:
			strchive_data = json.load(f)
	else:
		print(f"Downloading STRchive loci from {STRCHIVE_LOCI_URL}")
		strchive_data = requests.get(STRCHIVE_LOCI_URL, headers={"Cache-Control": "no-cache"}).json()

	strchive_loci = []
	for d in strchive_data:
		motifs = d["reference_motif_reference_orientation"]
		if len(motifs) != 1:
			raise ValueError(f"Expected exactly 1 reference motif for STRchive locus {d['id']}, found: {motifs}")

		strchive_loci.append({
			"LocusId": d["id"],
			"Gene": d["gene"],
			"Chrom": d["chrom"],
			"Start0Based": int(d["start_hg38"]),
			"End1Based": int(d["stop_hg38"]),
			"Motif": motifs[0],
		})

	return sorted(strchive_loci, key=lambda locus: locus["LocusId"])


def sequences_match(sequence1, sequence2):
	"""Return True if two equal-length sequences match base by base, with N matching any base"""
	return all(base1 == base2 or base1 == "N" or base2 == "N" for base1, base2 in zip(sequence1, sequence2))


def motifs_match(motif1, motif2):
	"""Return True if two motifs should be considered to describe the same repeat.

	For VNTRs (motifs longer than MAX_STR_MOTIF_SIZE) only the motif lengths need to match since the same VNTR is
	often represented using different motifs. For STRs, the motifs must be the same after taking rotations and
	reverse complement into account, with N in either motif matching any base (STRchive represents some motifs using
	N, such as the GCN motif of the PABPN1 locus).

	Args:
		motif1 (str): repeat unit.
		motif2 (str): repeat unit.

	Return:
		bool: whether the motifs match.
	"""
	motif1 = motif1.upper()
	motif2 = motif2.upper()
	if len(motif1) != len(motif2):
		return False

	if len(motif1) > MAX_STR_MOTIF_SIZE:
		return True

	for motif in motif2, reverse_complement(motif2):
		doubled_motif = motif + motif
		if any(sequences_match(motif1, doubled_motif[i:i+len(motif)]) for i in range(len(motif))):
			return True

	return False


def get_overlapping_trexplorer_loci(trexplorer_catalog, locus):
	"""Query the tabix-indexed TRExplorer catalog for loci that overlap the given STRchive locus.

	Args:
		trexplorer_catalog (pysam.TabixFile): the TRExplorer catalog BED file.
		locus (dict): STRchive locus record from get_strchive_loci(..).

	Return:
		list: dicts with the chrom, start_0based, end_1based and motif of each overlapping TRExplorer locus.
	"""
	try:
		rows = list(trexplorer_catalog.fetch(locus["Chrom"], locus["Start0Based"], locus["End1Based"]))
	except ValueError:
		print(f"WARNING: {locus['Chrom']} is not present in the TRExplorer catalog")
		return []

	overlapping_loci = []
	for row in rows:
		fields = row.split("\t")
		overlapping_loci.append({
			"Chrom": fields[0],
			"Start0Based": int(fields[1]),
			"End1Based": int(fields[2]),
			"Motif": fields[3],
		})

	return overlapping_loci


def boundary_distance(locus1, locus2):
	"""Return the larger of the start and end coordinate differences between two loci (in base pairs)"""
	return max(abs(locus1["Start0Based"] - locus2["Start0Based"]), abs(locus1["End1Based"] - locus2["End1Based"]))


def output(output_file, line):
	"""Print line to screen and write it to the file"""
	print(line)
	output_file.write(line + "\n")


def format_locus(locus):
	"""Return a string like 'chr1:94418421-94418444 (  8.0 x CCG    23bp)' for the given locus record"""
	span = locus["End1Based"] - locus["Start0Based"]
	return "%s:%d-%d (%5.1f x %s %6dbp)" % (
		locus["Chrom"], locus["Start0Based"], locus["End1Based"], span/len(locus["Motif"]), locus["Motif"], span)


def compare_catalogs(strchive_loci, trexplorer_catalog, output_path):
	"""Compare STRchive loci with the TRExplorer catalog and write the results to the output path.

	Args:
		strchive_loci (list): STRchive locus records from get_strchive_loci(..).
		trexplorer_catalog (pysam.TabixFile): the TRExplorer catalog BED file.
		output_path (str): path of the output text file.
	"""
	loci_present = []
	loci_absent = []
	loci_absent_without_any_overlap = []
	loci_with_different_boundaries = []

	with open(output_path, "wt") as output_file:
		output(output_file, "Comparing STRchive loci with the TRExplorer v2 catalog")
		output(output_file, "-"*100)

		for locus in strchive_loci:
			overlapping_loci = get_overlapping_trexplorer_loci(trexplorer_catalog, locus)
			matching_loci = [
				overlapping_locus for overlapping_locus in overlapping_loci
				if motifs_match(locus["Motif"], overlapping_locus["Motif"])
			]

			locus_label = f"{locus['LocusId']} ({locus['Gene']})"
			if matching_loci:
				best_match = min(matching_loci, key=lambda matching_locus: boundary_distance(locus, matching_locus))
				loci_present.append(locus)
				output(output_file, "%-30s %-45s is present in TRExplorer: %-45s %s" % (
					locus_label, format_locus(locus), format_locus(best_match),
					f"({len(matching_loci)} matching loci overlap it)" if len(matching_loci) > 1 else ""))

				if boundary_distance(locus, best_match) >= len(locus["Motif"]):
					loci_with_different_boundaries.append(locus)
			else:
				loci_absent.append(locus)
				if not overlapping_loci:
					loci_absent_without_any_overlap.append(locus)
					reason = "no TRExplorer loci overlap it"
				else:
					reason = "%d overlapping TRExplorer loci, but none with a matching motif: %s" % (
						len(overlapping_loci), ", ".join(sorted({d["Motif"] for d in overlapping_loci})))

				output(output_file, "%-30s %-45s is ABSENT from TRExplorer: %s" % (locus_label, format_locus(locus), reason))

		output(output_file, "-"*100)
		output(output_file, f"{len(loci_absent):,d} out of {len(strchive_loci):,d} STRchive loci are absent from the "
							f"TRExplorer catalog:")
		for locus in loci_absent:
			output(output_file, "%5s%-30s %s" % (" ", f"{locus['LocusId']} ({locus['Gene']})", format_locus(locus)))

		output(output_file, "-"*100)
		output(output_file, "Summary:")
		for label, loci in [
			("STRchive loci", strchive_loci),
			("present in the TRExplorer catalog", loci_present),
			("absent from the TRExplorer catalog", loci_absent),
			("absent and not overlapped by any TRExplorer locus", loci_absent_without_any_overlap),
			("present, but with boundaries that differ by 1 or more repeats", loci_with_different_boundaries),
		]:
			str_locus_count = sum(1 for locus in loci if len(locus["Motif"]) <= MAX_STR_MOTIF_SIZE)
			output(output_file, "%5s%4d out of %4d (%5.1f%%) %-60s  %3d STRs, %3d VNTRs" % (
				" ", len(loci), len(strchive_loci), 100*len(loci)/len(strchive_loci) if strchive_loci else 0,
				label, str_locus_count, len(loci) - str_locus_count))

	print(f"Wrote results to {output_path}")


def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--trexplorer-catalog-path", default=DEFAULT_TREXPLORER_CATALOG_PATH,
						help="Local path of the tabix-indexed TRExplorer v2 catalog BED file.")
	parser.add_argument("--strchive-loci-path",
						help="Optional local path of the STRchive-loci.json file (defaults to downloading it from "
							 "the STRchive github repo).")
	parser.add_argument("-o", "--output-path", default="TRExplorer_v2_STRchive_comparison.txt",
						help="Path of the output text file.")
	args = parser.parse_args()

	trexplorer_catalog_path = Path(args.trexplorer_catalog_path).expanduser()
	if not trexplorer_catalog_path.is_file():
		parser.error(f"TRExplorer catalog file not found: {trexplorer_catalog_path}")
	if not Path(f"{trexplorer_catalog_path}.tbi").is_file():
		parser.error(f"TRExplorer catalog tabix index not found: {trexplorer_catalog_path}.tbi. It can be created by "
					 f"running: tabix -p bed {trexplorer_catalog_path}")

	print(f"Using the TRExplorer catalog at {trexplorer_catalog_path}")
	compare_catalogs(get_strchive_loci(args.strchive_loci_path), pysam.TabixFile(str(trexplorer_catalog_path)),
					 args.output_path)


if __name__ == "__main__":
	main()
