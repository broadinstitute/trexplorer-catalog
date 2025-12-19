"""Convert the simple repeat track from UCSC to BED catalog format for comparison.

Source:
https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema
"""

import argparse
import gzip
import json
import os
import tqdm
import pandas as pd

# from: https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=rep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema
COLUMNS = [
	"bin",	 # Indexing field to speed chromosome range queries.
	"chrom",	 # Reference sequence chromosome or scaffold	
	"chromStart",	 # Start position in chromosome
	"chromEnd",	 # End position in chromosome
	"name",	 # Simple Repeats tag name
	"period",	 # Length of repeat unit
	"copyNum",	 # Mean number of copies of repeat
	"consensusSize",	 # Length of consensus sequence	
	"perMatch",	 # Percentage Match
	"perIndel",	 # Percentage Indel
	"score",
	"A",	 # Percent of A's in repeat unit			
	"C",	 # Percent of C's in repeat unit
	"G",	 # Percent of G's in repeat unit
	"T",	 # Percent of T's in repeat unit	
	"entropy",	 # Entropy
	"sequence",	 # Sequence of repeat unit element
]

def main():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-o", "--output-bed", help="Output BED file path")
	parser.add_argument("--simple-repeat-track-url", default="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz")
	args = parser.parse_args()

	print(f"Downloading and converting {args.simple_repeat_track_url} to BED format")
	if not args.output_bed:
		args.output_bed = os.path.basename(args.simple_repeat_track_url).replace(".txt.gz", ".bed.gz")
	if args.output_bed.endswith("gz"):
		args.output_bed = args.output_bed.replace(".gz", "")

	total = supercontig_loci_counter = output_counter = 0
	df = pd.read_table(args.simple_repeat_track_url, names=COLUMNS)

	with open(args.output_bed, "wt") as out:
		for _, row in tqdm.tqdm(df.iterrows(), unit=" rows", unit_scale=True):
			total += 1
			
			chrom = row["chrom"]
			if "_" in chrom:
				supercontig_loci_counter += 1
				continue
			start_0based = int(row["chromStart"])
			end = int(row["chromEnd"])
			motif = row["sequence"]
			out.write(f"{chrom}\t{start_0based}\t{end}\t{motif}\n")
			output_counter += 1

	os.system(f"bgzip -f {args.output_bed}")
	if supercontig_loci_counter:
		print(f"Skipped {supercontig_loci_counter:,d} out of {total:,d} loci "
			  f"({supercontig_loci_counter/total:.2%}) because they are on supercontigs")
	print(f"Wrote {output_counter:,d} rows to {args.output_bed}.gz")

if __name__ == "__main__":
	main()
