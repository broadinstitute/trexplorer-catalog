import unittest

import append_loci_not_already_in_catalog as script


def make_record(reference_region, locus_structure="(CAG)*", **extra_fields):
    return {"LocusId": reference_region, "ReferenceRegion": reference_region, "LocusStructure": locus_structure,
            "VariantType": "Repeat", **extra_fields}


class TestSelectLociNotAlreadyInCatalog(unittest.TestCase):

    def test_skips_loci_already_in_catalog_and_labels_the_others_with_the_source(self):
        catalog_records = [make_record("chr12:132111385-132111451")]
        records_to_append = [
            make_record("chr12:132111385-132111451"),            # already in the catalog
            make_record("chr12:132111385-132111451", "(GCA)*"),  # same region, different LocusStructure
            make_record("chr12:132111361-132111451", Diseases=[{"Symbol": "SCA"}]),
        ]

        new_records = script.select_loci_not_already_in_catalog(catalog_records, records_to_append, "Disease")

        self.assertEqual([(r["ReferenceRegion"], r["LocusStructure"]) for r in new_records], [
            ("chr12:132111385-132111451", "(GCA)*"), ("chr12:132111361-132111451", "(CAG)*")])
        for record in new_records:
            self.assertEqual(record["Source"], "Disease")
            self.assertEqual(record["FoundInDisease"], "Yes")
        self.assertEqual(new_records[1]["Diseases"], [{"Symbol": "SCA"}])

    def test_discard_extra_fields(self):
        new_records = script.select_loci_not_already_in_catalog(
            [], [make_record("chr12:132111361-132111451", Diseases=[{"Symbol": "SCA"}], PathogenicMotifs=["CAG"])],
            "Disease", discard_extra_fields=True)

        self.assertEqual(new_records, [{
            "LocusId": "chr12:132111361-132111451", "ReferenceRegion": "chr12:132111361-132111451",
            "LocusStructure": "(CAG)*", "VariantType": "Repeat", "Source": "Disease", "FoundInDisease": "Yes"}])


if __name__ == "__main__":
    unittest.main()
