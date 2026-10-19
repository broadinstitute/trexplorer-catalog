import json
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

import select_STRchive_loci_with_confident_disease_associations as script

STRCHIVE_LOCI = [
    {"id": "SCA1_ATXN1", "evidence": ["Definitive"], "year": "1993 [@pmid:8358429]"},
    {"id": "SCA_EP400", "evidence": ["Provisional"], "year": "2026"},
    {"id": "MODY8_CEL", "evidence": ["Provisional"], "year": "2005"},
    {"id": "DMD_DMD", "evidence": ["Refuted"], "year": "2016 [@pmid:27417533]"},
    {"id": "CPEO_POLG", "evidence": ["Disputed"], "year": None},
    {"id": "CCD_RUNX2", "evidence": ["Definitive"], "year": "Causation identified in 2015 [@pmid:26220009]"},
    {"id": "OPDM_LRP12", "evidence": ["Limited"], "year": "2024"},
    {"id": "FRA12A_DIP2B", "evidence": ["Disputed"], "year": "2024"},
]


class TestSelectSTRchiveLociWithConfidentDiseaseAssociations(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_get_year_first_reported(self):
        self.assertEqual(script.get_year_first_reported({"year": "2026 [@pmid:42482100]"}), 2026)
        self.assertEqual(script.get_year_first_reported({"year": "1991 [@pmid:2062380]; the first"}), 1991)
        self.assertIsNone(script.get_year_first_reported({"year": "Causation identified in 2015 [@pmid:26220009]"}))
        self.assertIsNone(script.get_year_first_reported({"year": None}))
        self.assertIsNone(script.get_year_first_reported({}))

    def test_select_STRchive_locus_ids_with_confident_disease_associations(self):
        # loci first reported in or after 2024 are kept whatever their evidence, like the website's New tag
        self.assertEqual(
            script.select_STRchive_locus_ids_with_confident_disease_associations(
                STRCHIVE_LOCI, script.DEFAULT_EVIDENCE_LEVELS, 2024),
            {"SCA1_ATXN1", "SCA_EP400", "CCD_RUNX2", "OPDM_LRP12", "FRA12A_DIP2B"})
        self.assertEqual(
            script.select_STRchive_locus_ids_with_confident_disease_associations(
                STRCHIVE_LOCI, ["Definitive"], 2027),
            {"SCA1_ATXN1", "CCD_RUNX2"})

    def test_main_filters_catalog_by_LocusId_and_metadata_by_id(self):
        strchive_loci_path = os.path.join(self.temp_dir, "STRchive-loci.json")
        with open(strchive_loci_path, "wt") as f:
            json.dump(STRCHIVE_LOCI, f)
        catalog_path = os.path.join(self.temp_dir, "catalog.json")
        with open(catalog_path, "wt") as f:
            json.dump([{"LocusId": "SCA1_ATXN1", "ReferenceRegion": "chr6:16327633-16327723"},
                       {"LocusId": "DMD_DMD", "ReferenceRegion": "chrX:31279056-31279121"}], f)

        for input_path, expected_ids in [
            (catalog_path, ["SCA1_ATXN1"]),
            (strchive_loci_path, ["SCA1_ATXN1", "SCA_EP400", "CCD_RUNX2", "OPDM_LRP12", "FRA12A_DIP2B"]),
        ]:
            output_path = os.path.join(self.temp_dir, "output.json")
            argv = ["script", "--strchive-loci-json", strchive_loci_path, "--min-year-first-reported-for-new-loci",
                    "2024", "-o", output_path, input_path]
            with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print"):
                script.main()

            with open(output_path) as f:
                self.assertEqual([r.get("LocusId", r.get("id")) for r in json.load(f)], expected_ids)


if __name__ == "__main__":
    unittest.main()
