import json
import os
from pprint import pprint

from str_analysis.utils.misc_utils import parse_interval

os.chdir(os.path.dirname(__file__))

unique_coordinates = set()
known_disease_associated_loci = []

examples = {}
with open("./variant_catalog_without_offtargets.GRCh38.json", "r") as f:
    gnomad_catalog_data = json.load(f)
    for record in gnomad_catalog_data:
        chrom, start_0based, end_1based = parse_interval(record["MainReferenceRegion"])
        motif = record["RepeatUnit"]
        key = (chrom, start_0based, end_1based)
        if key not in unique_coordinates:
            unique_coordinates.add(key)
            known_disease_associated_loci.append((chrom, start_0based, end_1based, motif))
        if record["LocusId"] == "HTT":
            for key, value in record.items():
                examples[key] = value

print(f"Added {len(known_disease_associated_loci):,d} loci from the gnomAD variant catalog")

examples = {}
total = 0
with open("./STRchive_loci.json", "r") as f:
    strchive_data = json.load(f)
    for record in strchive_data:
        chrom = record["chrom"]
        start_0based = record["start_hg38"]
        end_1based = record["stop_hg38"]
        motif = record["reference_motif_reference_orientation"][0]
        key = (chrom, start_0based, end_1based)
        if key not in unique_coordinates:
            unique_coordinates.add(key)
            known_disease_associated_loci.append((chrom, start_0based, end_1based, motif))
        total += 1
        for key, value in record.items():
            examples[key] = value

    #pprint(examples)


print(f"Added {total:,d} more loci from STRchive of which {len(known_disease_associated_loci) - total:,d} have unique coordinates")

output_path = "known_disease_associated_loci_v2.bed"
with open(output_path, "w") as f:
    for chrom, start_0based, end_1based, motif in known_disease_associated_loci:
        f.write(f"{chrom}\t{start_0based}\t{end_1based}\t{motif}\n")

print(f"Wrote {len(known_disease_associated_loci):,d} loci to {output_path}")
