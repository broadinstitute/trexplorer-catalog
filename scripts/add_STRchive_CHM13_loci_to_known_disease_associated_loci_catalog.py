"""Add STRchive's CHM13 definitions of the known disease-associated loci that didn't lift over from hg38.

The CHM13 catalog pipeline lifts the str-analysis known disease-associated loci catalog over from hg38, and uses the
result to set each locus's KnownDiseaseAssociatedLocus field. A few loci don't lift over (for example, because their
repeat sequence differs too much between the two references), and STRchive's CHM13 catalog is then the only source of
their definitions. For each such locus, this script adds a record that has the hg38 record's LocusId, gene and disease
fields, and the CHM13 coordinates and motif of the STRchive locus that STRchive maps to it (in its "gnomad" field), so
that the STRchive definition gets the same KnownDiseaseAssociatedLocus label as in the hg38 catalog.
"""

import argparse
import simplejson as json


def make_record_from_STRchive_coordinates(hg38_record, strchive_locus):
    """Returns a known disease-associated locus record with the hg38 record's fields and STRchive's CHM13 definition.

    STRchive's start_t2t and stop_t2t are a 0-based start and a 1-based end, like the ReferenceRegions of the
    str-analysis catalog (eg. DM2_CNBP has start_hg38 129172576 and stop_hg38 129172656, and the str-analysis CNBP
    locus has ReferenceRegion chr3:129172576-129172656).

    Args:
        hg38_record (dict): record of the str-analysis known disease-associated loci catalog on hg38
        strchive_locus (dict): record of STRchive-loci.json

    Returns:
        dict: the new record
    """
    reference_region = f"{strchive_locus['chrom']}:{strchive_locus['start_t2t']}-{strchive_locus['stop_t2t']}"
    motif = strchive_locus["reference_motif_reference_orientation"][0]
    return {
        "LocusId": hg38_record["LocusId"],
        "LocusStructure": f"({motif})*",
        "ReferenceRegion": reference_region,
        "MainReferenceRegion": reference_region,
        "VariantType": "Repeat",
        "RepeatUnit": motif,
        "Gene": hg38_record.get("Gene"),
        "GeneRegion": hg38_record.get("GeneRegion"),
        "GeneId": hg38_record.get("GeneId"),
        "Diseases": hg38_record["Diseases"],
        # the LocusId of the hg38 TRExplorer catalog's definition of this locus, which the CHM13 catalog's
        # Hg38LocusId falls back to for this locus, since its hg38 definition didn't lift over
        "TRExplorerV2": hg38_record.get("TRExplorerV2"),
        "STRchiveId": strchive_locus["id"],
    }


def make_records_for_loci_that_did_not_lift_over(hg38_records, lifted_records, strchive_loci):
    """Returns records with STRchive's CHM13 definitions of the hg38 disease-associated loci that didn't lift over.

    Alternative definitions (LocusIds with "__", like RUNX2__6-45422750-45422801-GCN) are skipped, since they describe
    a locus that their primary definition already represents.

    Args:
        hg38_records (list): records of the str-analysis known disease-associated loci catalog on hg38
        lifted_records (list): records of the same catalog lifted over to CHM13
        strchive_loci (list): records of STRchive-loci.json

    Returns:
        2-tuple: (list of new records, list of LocusIds that STRchive has no CHM13 definition for)
    """
    lifted_locus_ids = {r["LocusId"] for r in lifted_records}
    strchive_loci_by_locus_id = {}
    for strchive_locus in strchive_loci:
        if strchive_locus.get("start_t2t") and strchive_locus.get("stop_t2t"):
            for locus_id in strchive_locus.get("gnomad") or []:
                strchive_loci_by_locus_id.setdefault(locus_id, strchive_locus)

    new_records = []
    locus_ids_without_STRchive_definition = []
    for hg38_record in hg38_records:
        locus_id = hg38_record["LocusId"]
        if not hg38_record.get("Diseases") or "__" in locus_id or locus_id in lifted_locus_ids:
            continue
        if locus_id in strchive_loci_by_locus_id:
            new_records.append(make_record_from_STRchive_coordinates(hg38_record, strchive_loci_by_locus_id[locus_id]))
        else:
            locus_ids_without_STRchive_definition.append(locus_id)

    return new_records, locus_ids_without_STRchive_definition


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--strchive-loci-json", required=True, help="Path of STRchive's data/STRchive-loci.json")
    parser.add_argument("--hg38-known-disease-associated-loci", required=True, help="Path of the str-analysis known "
                        "disease-associated loci catalog on hg38 (variant_catalog_without_offtargets.GRCh38.json)")
    parser.add_argument("-o", "--output-path", required=True, help="Path of the output catalog JSON file")
    parser.add_argument("lifted_known_disease_associated_loci_path", help="Path of the known disease-associated loci "
                        "catalog lifted over to CHM13 by lift_over_catalog.py")
    args = parser.parse_args()

    with open(args.hg38_known_disease_associated_loci) as f:
        hg38_records = json.load(f)
    with open(args.lifted_known_disease_associated_loci_path) as f:
        lifted_records = json.load(f)
    with open(args.strchive_loci_json) as f:
        strchive_loci = json.load(f)

    new_records, locus_ids_without_STRchive_definition = make_records_for_loci_that_did_not_lift_over(
        hg38_records, lifted_records, strchive_loci)
    for record in new_records:
        print(f"Adding STRchive locus {record['STRchiveId']} as {record['LocusId']} @ {record['ReferenceRegion']}")
    for locus_id in locus_ids_without_STRchive_definition:
        print(f"WARNING: {locus_id} didn't lift over, and STRchive has no CHM13 definition for it")

    with open(args.output_path, "wt") as f:
        json.dump(lifted_records + new_records, f, indent=4)
    print(f"Wrote {len(lifted_records) + len(new_records):,d} loci to {args.output_path}")


if __name__ == "__main__":
    main()
