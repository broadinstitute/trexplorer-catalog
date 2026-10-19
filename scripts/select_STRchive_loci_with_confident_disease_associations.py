"""Keep only the STRchive loci whose disease association is confident enough for the catalog.

STRchive (https://strchive.org) classifies the evidence for each locus's disease association as Definitive, Strong,
Moderate, Limited, Disputed, Refuted, No Known Relationship or Provisional. This keeps the loci whose evidence is one of
the --evidence-levels (Definitive, Strong, Moderate or Limited by default), plus the loci that STRchive's website tags as
New. The website tags a locus as New when it was first reported (the year at the start of its "year" field) at most 2
years before the current year, so newly reported loci are kept even while their evidence is still Provisional. Since
that depends on the date the website is viewed, this script instead takes the earliest year that counts as New
(--min-year-first-reported-for-new-loci), which keeps the output reproducible.

The input can be STRchive's locus metadata (data/STRchive-loci.json, whose records have an "id" field) or one of its
catalogs (like data/catalogs/STRchive-disease-loci.T2T-chm13.expansionhunter.json, whose LocusIds are those ids).
"""

import argparse
import re
import simplejson as json

DEFAULT_EVIDENCE_LEVELS = ("Definitive", "Strong", "Moderate", "Limited")


def get_year_first_reported(strchive_locus):
    """Returns the year at the start of an STRchive locus's "year" field, or None if the field doesn't start with one.

    Like the parseInt() call that STRchive's website uses for its New tag, this reads only leading digits. So
    "2026 [@pmid:42482100]" gives 2026, while "Causation identified in 2015 [@pmid:26220009]; ..." gives None.

    Args:
        strchive_locus (dict): record from STRchive-loci.json

    Returns:
        int: the year, or None
    """
    match = re.match(r"\s*(\d+)", strchive_locus.get("year") or "")
    return int(match.group(1)) if match else None


def select_STRchive_locus_ids_with_confident_disease_associations(
        strchive_loci, evidence_levels, min_year_first_reported_for_new_loci):
    """Returns the ids of the STRchive loci to keep.

    Args:
        strchive_loci (list): records from STRchive-loci.json
        evidence_levels (iterable): keep the loci whose evidence includes one of these levels
        min_year_first_reported_for_new_loci (int): also keep the loci first reported in or after this year

    Returns:
        set: ids of the loci to keep
    """
    selected_ids = set()
    for strchive_locus in strchive_loci:
        year_first_reported = get_year_first_reported(strchive_locus)
        if set(strchive_locus.get("evidence") or []) & set(evidence_levels) or (
                year_first_reported is not None and year_first_reported >= min_year_first_reported_for_new_loci):
            selected_ids.add(strchive_locus["id"])

    return selected_ids


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=__doc__.split("\n")[0])
    parser.add_argument("--strchive-loci-json", required=True, help="Path of STRchive's data/STRchive-loci.json")
    parser.add_argument("--evidence-levels", nargs="+", default=list(DEFAULT_EVIDENCE_LEVELS),
                        help="Keep the loci whose evidence includes one of these levels")
    parser.add_argument("--min-year-first-reported-for-new-loci", type=int, required=True,
                        help="Also keep the loci first reported in or after this year, which reproduces STRchive's "
                             "New tag when set to 2 years before the year of the STRchive version")
    parser.add_argument("-o", "--output-json-path", required=True, help="Path of the output JSON file")
    parser.add_argument("input_json_path", help="Path of STRchive-loci.json or of an STRchive JSON catalog whose "
                        "LocusIds are STRchive ids")
    args = parser.parse_args()

    with open(args.strchive_loci_json) as f:
        strchive_loci = json.load(f)

    selected_ids = select_STRchive_locus_ids_with_confident_disease_associations(
        strchive_loci, args.evidence_levels, args.min_year_first_reported_for_new_loci)
    for strchive_locus in strchive_loci:
        if strchive_locus["id"] not in selected_ids:
            print(f"Dropping {strchive_locus['id']} (evidence: {', '.join(strchive_locus.get('evidence') or [])}, "
                  f"first reported: {get_year_first_reported(strchive_locus)})")

    with open(args.input_json_path) as f:
        records = json.load(f)

    output_records = [r for r in records if (r["LocusId"] if "LocusId" in r else r["id"]) in selected_ids]
    with open(args.output_json_path, "wt") as f:
        json.dump(output_records, f, indent=4)

    print(f"Wrote {len(output_records)} out of {len(records)} records to {args.output_json_path}")


if __name__ == "__main__":
    main()
