"""Append the loci of another catalog that the input catalog doesn't already define, then sort.

The hg38 pipeline uses this for its extended loci catalog, which holds wider definitions for loci whose boundaries
stopped short of the repeat they sit in (see extended_loci/generate_extended_loci_catalog.sh). The CHM13 pipeline uses
it for the known disease-associated loci, some of which have an alternative definition that overlaps the primary one
(like EP400 and EP400__12-132062524-132062611-CAG). These loci are appended rather than merged because two definitions
of the same locus necessarily overlap, and str_analysis.merge_loci would discard one of the two.
Sorting matters: the annotation steps that follow group overlapping loci and reject a catalog whose coordinates ever
decrease, and an appended definition can start before the locus it overlaps.
"""

import argparse

from str_analysis.extend_str_catalog_boundaries import sort_catalog_records
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator
from str_analysis.utils.export_json import export_json

# The fields that str_analysis.merge_loci --discard-extra-fields-from-input-catalogs keeps
CORE_FIELDS = ("LocusId", "ReferenceRegion", "LocusStructure", "VariantType")


def select_loci_not_already_in_catalog(catalog_records, records_to_append, source_name, discard_extra_fields=False):
    """Returns the records to append whose ReferenceRegion and LocusStructure aren't already in the catalog.

    Each returned record is labeled with source_name as its Source and FoundIn field. Without this they reach the
    sources-stats table as "Unknown", since it reads record["Source"], and the catalog would not say where they came
    from. The FoundIn field is set to "Yes", the value str_analysis.merge_loci --add-found-in-fields uses, since the
    source can also be one of the merge's sources.

    Args:
        catalog_records (list): records of the input catalog
        records_to_append (list): records of the catalog of loci to append
        source_name (str): name of the source that the appended records come from
        discard_extra_fields (bool): keep only the CORE_FIELDS of the records to append

    Returns:
        list: records to append to the catalog
    """
    already_in_catalog = {
        (r["ReferenceRegion"], r["LocusStructure"]) for r in catalog_records
    }
    new_records = [
        {k: v for k, v in r.items() if not discard_extra_fields or k in CORE_FIELDS} for r in records_to_append
        if (r["ReferenceRegion"], r["LocusStructure"]) not in already_in_catalog
    ]
    for record in new_records:
        record["Source"] = source_name
        record[f"FoundIn{source_name}"] = "Yes"

    return new_records


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--catalog-of-loci-to-append", required=True,
                        help="Path of the catalog (BED or JSON) whose loci to append")
    parser.add_argument("--source-name", required=True, help="Source name to give the appended loci, eg. "
                        "TRExplorerV2.1:ExtendedLocusBoundaries")
    parser.add_argument("--discard-extra-fields", action="store_true", help="Keep only " + ", ".join(CORE_FIELDS) +
                        " of the appended loci, like str_analysis.merge_loci --discard-extra-fields-from-input-catalogs")
    parser.add_argument("-o", "--output-catalog-json-path", required=True, help="Path of the output catalog JSON file")
    parser.add_argument("catalog_json_path", help="Path of the input catalog JSON file")
    args = parser.parse_args()

    records = list(get_variant_catalog_iterator(args.catalog_json_path))
    print(f"Read {len(records):,d} loci from {args.catalog_json_path}")

    records_to_append = list(get_variant_catalog_iterator(args.catalog_of_loci_to_append))
    print(f"Read {len(records_to_append):,d} loci to append from {args.catalog_of_loci_to_append}")

    new_records = select_loci_not_already_in_catalog(
        records, records_to_append, args.source_name, args.discard_extra_fields)
    if len(new_records) < len(records_to_append):
        print(f"Skipped {len(records_to_append) - len(new_records):,d} loci already defined in the catalog")

    export_json(sort_catalog_records(records + new_records), args.output_catalog_json_path)
    print(f"Appended {len(new_records):,d} loci and wrote {len(records) + len(new_records):,d} loci")


if __name__ == "__main__":
    main()
