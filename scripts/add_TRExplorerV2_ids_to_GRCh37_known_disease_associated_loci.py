"""Copy the TRExplorerV2 field of the GRCh38 known disease-associated loci catalog to the GRCh37 one.

The TRExplorerV2 field of the str-analysis known disease-associated loci catalog holds the LocusId of the TRExplorer hg38
catalog's definition of each locus. The hg19 pipeline uses it to find the hg38 LocusId of a known disease-associated
locus that matches no lifted hg38 locus, but the GRCh37 version of the catalog has no such field. So this copies each
GRCh38 record's TRExplorerV2 value to the GRCh37 record that defines the same locus: the one with the same LocusId, or,
for an alternative definition whose LocusId ends in its coordinates and motif after "__" (eg.
PABPN1__14-23790681-23790712-GCG on GRCh37 and PABPN1__14-23321472-23321503-GCG on GRCh38), the one whose LocusId is the
same apart from its coordinates.
"""

import argparse
import simplejson as json

from lift_over_catalog import LOCUS_ID_WITH_COORDINATES_REGEX


def get_locus_id_without_coordinates(locus_id):
    """Returns the LocusId without the coordinates of an alternative definition's LocusId.

    Args:
        locus_id (str): a LocusId of the known disease-associated loci catalog, eg. PABPN1__14-23790681-23790712-GCG

    Returns:
        str: eg. PABPN1__GCG, or the LocusId itself if it doesn't end in coordinates (eg. HTT)
    """
    match = LOCUS_ID_WITH_COORDINATES_REGEX.match(locus_id)
    return f"{match.group(1)}{match.group(2)}" if match else locus_id


def add_TRExplorerV2_ids(grch37_records, grch38_records):
    """Sets the TRExplorerV2 field of each GRCh37 record that has a GRCh38 counterpart with one.

    Args:
        grch37_records (list): records of the GRCh37 known disease-associated loci catalog, which this modifies
        grch38_records (list): records of the GRCh38 known disease-associated loci catalog

    Returns:
        int: the number of GRCh37 records that got a TRExplorerV2 value
    """
    TRExplorerV2_ids = {}
    for record in grch38_records:
        if not record.get("TRExplorerV2"):
            continue
        key = get_locus_id_without_coordinates(record["LocusId"])
        if key in TRExplorerV2_ids:
            raise ValueError(f"More than one GRCh38 record has the LocusId {key}, apart from its coordinates")
        TRExplorerV2_ids[key] = record["TRExplorerV2"]

    num_records_with_TRExplorerV2_id = 0
    for record in grch37_records:
        key = get_locus_id_without_coordinates(record["LocusId"])
        if key in TRExplorerV2_ids:
            record["TRExplorerV2"] = TRExplorerV2_ids[key]
            num_records_with_TRExplorerV2_id += 1

    return num_records_with_TRExplorerV2_id


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--grch38-known-disease-associated-loci", required=True,
                        help="Path of the GRCh38 known disease-associated loci catalog, whose records have a "
                             "TRExplorerV2 field")
    parser.add_argument("-o", "--output-catalog-json-path", required=True, help="Path of the output catalog JSON file")
    parser.add_argument("grch37_known_disease_associated_loci", help="Path of the GRCh37 known disease-associated loci "
                        "catalog")
    args = parser.parse_args()

    with open(args.grch37_known_disease_associated_loci) as f:
        grch37_records = json.load(f)
    with open(args.grch38_known_disease_associated_loci) as f:
        grch38_records = json.load(f)

    num_records_with_TRExplorerV2_id = add_TRExplorerV2_ids(grch37_records, grch38_records)

    with open(args.output_catalog_json_path, "wt") as f:
        json.dump(grch37_records, f, indent=4)

    print(f"Wrote {len(grch37_records):,d} records to {args.output_catalog_json_path}, of which "
          f"{num_records_with_TRExplorerV2_id:,d} got a TRExplorerV2 value from the GRCh38 catalog")


if __name__ == "__main__":
    main()
