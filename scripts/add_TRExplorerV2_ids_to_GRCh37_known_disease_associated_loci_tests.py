import os
import shutil
import simplejson as json
import sys
import tempfile
import unittest
from unittest import mock

import add_TRExplorerV2_ids_to_GRCh37_known_disease_associated_loci as script


class TestAddTRExplorerV2IdsToGRCh37KnownDiseaseAssociatedLoci(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_get_locus_id_without_coordinates(self):
        self.assertEqual(script.get_locus_id_without_coordinates("PABPN1__14-23790681-23790712-GCG"), "PABPN1__GCG")
        self.assertEqual(script.get_locus_id_without_coordinates("PABPN1__14-23321472-23321503-GCG"), "PABPN1__GCG")
        self.assertEqual(script.get_locus_id_without_coordinates("HTT"), "HTT")
        self.assertEqual(script.get_locus_id_without_coordinates("HOXA13_1"), "HOXA13_1")

    def test_add_TRExplorerV2_ids(self):
        grch37_records = [
            {"LocusId": "HTT"},
            {"LocusId": "PABPN1__14-23790681-23790712-GCG"},
            {"LocusId": "DMD"},          # its GRCh38 record has no TRExplorerV2 value
            {"LocusId": "ONLY_ON_GRCH37"},
        ]
        grch38_records = [
            {"LocusId": "HTT", "TRExplorerV2": "4-3074876-3074933-CAG"},
            {"LocusId": "PABPN1__14-23321472-23321503-GCG", "TRExplorerV2": "14-23321472-23321490-GCG"},
            {"LocusId": "DMD"},
            {"LocusId": "ONLY_ON_GRCH38", "TRExplorerV2": "1-100-130-CAG"},
        ]
        self.assertEqual(script.add_TRExplorerV2_ids(grch37_records, grch38_records), 2)
        self.assertEqual(grch37_records, [
            {"LocusId": "HTT", "TRExplorerV2": "4-3074876-3074933-CAG"},
            {"LocusId": "PABPN1__14-23790681-23790712-GCG", "TRExplorerV2": "14-23321472-23321490-GCG"},
            {"LocusId": "DMD"},
            {"LocusId": "ONLY_ON_GRCH37"},
        ])

    def test_add_TRExplorerV2_ids_rejects_ambiguous_GRCh38_records(self):
        grch38_records = [
            {"LocusId": "PABPN1__14-23321472-23321503-GCG", "TRExplorerV2": "a"},
            {"LocusId": "PABPN1__14-23321500-23321530-GCG", "TRExplorerV2": "b"},
        ]
        with self.assertRaises(ValueError):
            script.add_TRExplorerV2_ids([], grch38_records)

    def test_main(self):
        grch37_path = os.path.join(self.temp_dir, "grch37.json")
        grch38_path = os.path.join(self.temp_dir, "grch38.json")
        output_path = os.path.join(self.temp_dir, "output.json")
        with open(grch37_path, "wt") as f:
            json.dump([{"LocusId": "HTT", "Diseases": [{"Symbol": "HD"}]}], f)
        with open(grch38_path, "wt") as f:
            json.dump([{"LocusId": "HTT", "TRExplorerV2": "4-3074876-3074933-CAG"}], f)

        argv = ["script", "--grch38-known-disease-associated-loci", grch38_path, "-o", output_path, grch37_path]
        with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print"):
            script.main()

        with open(output_path) as f:
            self.assertEqual(json.load(f), [
                {"LocusId": "HTT", "Diseases": [{"Symbol": "HD"}], "TRExplorerV2": "4-3074876-3074933-CAG"}])


if __name__ == "__main__":
    unittest.main()
