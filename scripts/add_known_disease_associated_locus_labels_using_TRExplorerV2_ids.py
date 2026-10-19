"""Label the catalog loci that the exact-boundary disease locus match missed, using the TRExplorerV2 field.

str_analysis.annotate_and_filter_str_catalog sets KnownDiseaseAssociatedLocus only when a catalog locus has exactly
the same boundaries as a known disease-associated locus. A few known disease-associated loci are defined slightly
differently from the catalog's locus at the same position. For example, PRE-MIR7-2 is chr15:88569433-88569449 (4 pure
TTTG repeats), while the catalog locus 15-88569433-88569452-TTTG also includes the trailing partial repeat. The known
disease-associated loci catalog's TRExplorerV2 field gives the LocusId of the corresponding TRExplorer catalog locus,
so this script labels each catalog locus that still has no KnownDiseaseAssociatedLocus and whose LocusId equals the
TRExplorerV2 value of a known disease-associated locus. Labels from the exact-boundary match are never changed, and a
known disease-associated locus that the exact-boundary match already labeled a catalog locus with isn't used again, so
that no two catalog loci get the same label.

For a catalog lifted over from the hg38 TRExplorer catalog, --catalog-locus-id-field names the field that holds each
locus's hg38 LocusId, which is what the TRExplorerV2 values refer to.
"""

import argparse
import gzip
import ijson
import simplejson as json


def get_labels_already_in_catalog(catalog_json_path):
    """Returns the KnownDiseaseAssociatedLocus values of the catalog's loci.

    Args:
        catalog_json_path (str): path of the JSON catalog

    Returns:
        set: the labels
    """
    fopen = gzip.open if catalog_json_path.endswith("gz") else open
    with fopen(catalog_json_path, "rb") as f:
        return {label for label in ijson.items(f, "item.KnownDiseaseAssociatedLocus") if label}


def build_known_disease_associated_locus_ids_by_TRExplorerV2_id(known_disease_associated_loci,
                                                                 labels_already_in_catalog=()):
    """Maps each TRExplorerV2 id to the LocusId of the known disease-associated locus that has it.

    Entries without diseases are skipped, as in str_analysis.annotate_and_filter_str_catalog, and so are the entries
    whose LocusId is already a label in the catalog. When several of the remaining entries share a TRExplorerV2 id (eg.
    SAMD12 and SAMD12__8-118366812-118366918-AAAAT), the first one in the file is used.

    Args:
        known_disease_associated_loci (list): records of the known disease-associated loci catalog
        labels_already_in_catalog (set): KnownDiseaseAssociatedLocus values that catalog loci already have

    Returns:
        dict: TRExplorerV2 id => LocusId of the known disease-associated locus
    """
    known_disease_associated_locus_ids_by_TRExplorerV2_id = {}
    for known_locus in known_disease_associated_loci:
        if not known_locus.get("Diseases") or not known_locus.get("TRExplorerV2"):
            continue
        if known_locus["LocusId"] in labels_already_in_catalog:
            continue
        known_disease_associated_locus_ids_by_TRExplorerV2_id.setdefault(
            known_locus["TRExplorerV2"], known_locus["LocusId"])

    return known_disease_associated_locus_ids_by_TRExplorerV2_id


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Label the catalog loci that the exact-boundary disease locus match missed, using the "
                    "TRExplorerV2 field of the known disease-associated loci catalog.")
    parser.add_argument("--known-disease-associated-loci", required=True,
                        help="Path of the known disease-associated loci JSON catalog "
                             "(eg. str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json)")
    parser.add_argument("--catalog-locus-id-field", default="LocusId",
                        help="Field of each catalog record to match against the TRExplorerV2 values, eg. the field "
                             "that holds the hg38 LocusId of a catalog lifted over from hg38")
    parser.add_argument("-o", "--output-catalog-json-path", required=True,
                        help="Path of the output catalog JSON file")
    parser.add_argument("catalog_json_path", help="Path of the JSON catalog that "
                        "str_analysis.annotate_and_filter_str_catalog annotated with KnownDiseaseAssociatedLocus")
    args = parser.parse_args()

    with open(args.known_disease_associated_loci, "rt") as f:
        known_disease_associated_loci = json.load(f)

    if not build_known_disease_associated_locus_ids_by_TRExplorerV2_id(known_disease_associated_loci):
        parser.error(f"None of the disease-associated loci in {args.known_disease_associated_loci} have a "
                     f"TRExplorerV2 field")

    known_disease_associated_locus_ids_by_TRExplorerV2_id = build_known_disease_associated_locus_ids_by_TRExplorerV2_id(
        known_disease_associated_loci, get_labels_already_in_catalog(args.catalog_json_path))

    input_locus_counter = 0
    known_disease_associated_locus_ids_seen_as_labels = set()
    fopen = gzip.open if args.catalog_json_path.endswith("gz") else open
    with fopen(args.catalog_json_path, "rb") as f:
        f2open = gzip.open if args.output_catalog_json_path.endswith("gz") else open
        with f2open(args.output_catalog_json_path, "wt") as f2:
            f2.write("[")
            for i, record in enumerate(ijson.items(f, "item")):
                input_locus_counter += 1
                catalog_locus_id = record.get(args.catalog_locus_id_field)
                if not record.get("KnownDiseaseAssociatedLocus") and (
                        catalog_locus_id in known_disease_associated_locus_ids_by_TRExplorerV2_id):
                    record["KnownDiseaseAssociatedLocus"] = known_disease_associated_locus_ids_by_TRExplorerV2_id[
                        catalog_locus_id]
                    print(f"Labeled {record['LocusId']} as {record['KnownDiseaseAssociatedLocus']} based on its "
                          f"TRExplorerV2 id {catalog_locus_id}")

                if record.get("KnownDiseaseAssociatedLocus"):
                    known_disease_associated_locus_ids_seen_as_labels.add(record["KnownDiseaseAssociatedLocus"])

                if i > 0:
                    f2.write(", ")
                f2.write(json.dumps(record, use_decimal=True, indent=4))
            f2.write("]")

    for known_locus in known_disease_associated_loci:
        if known_locus.get("Diseases") and known_locus["LocusId"] not in known_disease_associated_locus_ids_seen_as_labels:
            print(f"No catalog locus is labeled {known_locus['LocusId']} (TRExplorerV2: "
                  f"{known_locus.get('TRExplorerV2')})")

    print(f"Wrote {input_locus_counter:,d} records to {args.output_catalog_json_path}")


if __name__ == "__main__":
    main()
