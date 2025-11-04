"""This script determines if there are any new disease-associated loci that need to be added to the catalog by
comparing the existing catalog with STRchive and the gnomAD catalog of known TR loci"""

#%%

import argparse
import collections
from datetime import datetime
import gzip
import intervaltree
import json
import os
import requests
import tqdm

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif


def get_json_from_url(url):
	"""Download and parse json from a url"""
	return requests.get(url, headers={'Cache-Control': 'no-cache'}).json()


def get_strchive_entries_not_in_gnomad(gnomad_catalog):
    """Download the locus metadata json from the STRchive github repo which underlies https://harrietdashnow.com/STRchive

    Return:
        dict: locus_id -> locus_info
    """

    gnomad_loci = {d["LocusId"] for d in gnomad_catalog}

    print("Downloading loci from STRchive")
    strchive_data = get_json_from_url(
        #"https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/STRchive-loci.json"
        "https://raw.githubusercontent.com/dashnowlab/STRchive/1954f32525846ae91363d3f308b840cbe772f01b/data/STRchive-loci.json"
    )

    strchive_id_to_gnomad_map = {
        "pre-MIR7-2_CHNG3": "PRE-MIR7-2",
        "HOXA13_1": "HOXA13",
        "HOXA13_2": "HOXA13",
        "HOXA13_3": "HOXA13",
        "C9orf72": "C9ORF72",
        "ARX_1": "ARX",
        "ARX_2": "ARX",
    }

    
    strchive_loci_not_in_gnomad = []
    previously_seen_gnomad_genes = set()
    for d in strchive_data:
        if d["id"] in strchive_id_to_gnomad_map:
            d["gnomad"] = strchive_id_to_gnomad_map[d["id"]]
        elif d["gene"] in gnomad_loci:
            d["gnomad"] = d["gene"]
        elif len(d["gnomad"]) > 0:
            d["gnomad"] = d["gnomad"][0]
        else:
            d["gnomad"] = None
            print(f"WARNING: STRchive locus is absent from gnomAD: {d['id']}")
            strchive_loci_not_in_gnomad.append(d)
            continue

        if d["gnomad"]:
            if d["gnomad"] in previously_seen_gnomad_genes:
                print(f"WARNING: Duplicate gnomad field: {d['gnomad']}")

            previously_seen_gnomad_genes.add(d["gnomad"])

    results = {
        d["gene"]: d for d in strchive_data if d.get("gnomad") is None
    }


    for d in results.values():
        assert len(d["reference_motif_reference_orientation"]) == 1, d
        d["reference_motif_reference_orientation"] = d["reference_motif_reference_orientation"][0]
        d["CanonicalMotif"] = compute_canonical_motif(d["reference_motif_reference_orientation"])

        d["Diseases"] = [{
            "Name": d["disease"],
            "Symbol": d["disease_id"],
            "Inheritance": ",".join(d["inheritance"]),
            "Onset": d["age_onset"],
            "NormalMax": d["benign_max"],
            "IntermediateRange": f"{d['intermediate_min']}-{d['intermediate_max']}",
            "PathogenicMin": d["pathogenic_min"],
        }]

    return results



#%%

# date timestamp (YYYY_MM_DD)
date_timestamp = datetime.now().strftime('%Y_%m_%d')
default_output_path = f"known_disease_associated_loci_to_add.{date_timestamp}.json"

#parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#parser.add_argument("-o", "--output-path", default=default_output_path, help="output JSON path")
#args = parser.parse_args()

# collect known loci from gnomAD and STRchive
gnomad_catalog_data = get_json_from_url(
    "https://raw.githubusercontent.com/broadinstitute/str-analysis/refs/heads/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json")

strchive_data = get_strchive_entries_not_in_gnomad(gnomad_catalog_data)
del strchive_data["POLG"]  # since it has insufficient evidence

known_locus_chroms = set()
known_locus_intervals = {}
print(f"Adding {len(gnomad_catalog_data)} loci from gnomAD")
for locus in gnomad_catalog_data:
    chrom, start_0based, end_1based = parse_interval(locus["MainReferenceRegion"])
    chrom = chrom.replace("chr", "")
    known_locus_chroms.add(chrom)
    known_locus_intervals[(chrom, start_0based, end_1based, len(locus["RepeatUnit"]))] = locus["LocusId"]

print(f"Adding {len(strchive_data)} more loci from STRchive")
for locus_id, locus in strchive_data.items():
    chrom = locus["chrom"].replace("chr", "")
    known_locus_chroms.add(chrom)
    known_locus_intervals[(chrom, int(locus["start_hg38"]), int(locus["stop_hg38"]), len(locus["reference_motif_reference_orientation"]))] = locus_id

print(f"Found {len(known_locus_chroms)} chromosomes have known loci")

#%%



# figure out which loci are already well-represented in the catalog
trexplorer_locus_interval_trees = collections.defaultdict(intervaltree.IntervalTree)

os.chdir(os.path.expanduser("~/code/tandem-repeat-catalogs"))
with gzip.open("results__2025-11-03/1_to_1000bp_motifs/repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz", "rt") as f:
    for line in tqdm.tqdm(f, unit=" records", unit_scale=True):
        fields = line.strip("\n").split("\t")
        chrom = fields[0].replace("chr", "")

        start_0based = int(fields[1])
        end_1based = int(fields[2])
        motif = fields[3]

        key = (chrom, start_0based, end_1based, len(motif))
        if key in known_locus_intervals:
            print(f"Found TRExplorer locus {key} in known loci")
            del known_locus_intervals[key]
            continue
    
        trexplorer_locus_interval_trees[chrom].add(intervaltree.Interval(start_0based, end_1based, data=motif))

#%%

missing_loci = []
for (chrom, start_0based, end_1based, motif_length), locus_id in known_locus_intervals.items():
    print('--------------------------------')
    print(f"Known locus: {locus_id}: {chrom}:{start_0based}-{end_1based} ({end_1based - start_0based}bp) {motif_length}bp")
    for overlapping_interval in trexplorer_locus_interval_trees[chrom].overlap(start_0based, end_1based):
        if len(overlapping_interval.data) != motif_length:
            continue

        trexplorer_motif = overlapping_interval.data
        size_difference = abs((overlapping_interval.end - overlapping_interval.begin) - (end_1based - start_0based))

        if size_difference//motif_length < 1:
            print(f"Found matching TRExplorer locus for {locus_id} with size {overlapping_interval.end - overlapping_interval.begin}bp and size difference of {size_difference}bp which is {size_difference/motif_length:.2f} x {trexplorer_motif} repeats")
            break

        print(f"Found overlapping TRExplorer locus for {locus_id} with size {overlapping_interval.end - overlapping_interval.begin}bp and size difference of {size_difference}bp which is {size_difference/motif_length:.2f} x {trexplorer_motif} repeats")
    else:
        print(f"No matching TRExplorer locus found for {locus_id}")

        missing_loci.append({
            "locus_id": locus_id,
            "chrom": chrom,
            "start_0based": start_0based,
            "end_1based": end_1based,
            "motif_length": motif_length,
        })

# Only the updated EP400, RUNX2 defintions, as well MUC1 60bp motif VNTR(?) need to be added to the catalog:
# - EP400: 12:132062524-132062611 (87bp) 3bp motif
# - RUNX2: 6:45422750-45422801 (51bp) 3bp motif
# - MUC1: 1:155188505-155192239 (3734bp) 60bp motif

# TODO review definitions for these loci
# TODO add other VNTRs


