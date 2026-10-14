import gzip
import os
import shutil
import simplejson as json
import sys
import tempfile
import unittest
from unittest import mock

import add_known_disease_associated_locus_labels_using_TRExplorerV2_ids as script


def make_known_locus(locus_id, main_reference_region, trexplorer_v2_id, diseases=("some disease",)):
    known_locus = {
        "LocusId": locus_id,
        "MainReferenceRegion": main_reference_region,
        "Diseases": list(diseases),
    }
    if trexplorer_v2_id is not None:
        known_locus["TRExplorerV2"] = trexplorer_v2_id
    return known_locus


class TestBuildKnownDiseaseAssociatedLocusIdsByTRExplorerV2Id(unittest.TestCase):

    def test_maps_TRExplorerV2_ids_to_locus_ids(self):
        lookup = script.build_known_disease_associated_locus_ids_by_TRExplorerV2_id([
            make_known_locus("PRE-MIR7-2", "chr15:88569433-88569449", "15-88569433-88569452-TTTG"),
            make_known_locus("HTT", "chr4:3074876-3074933", "4-3074876-3074933-CAG"),
        ])
        self.assertEqual(lookup, {
            "15-88569433-88569452-TTTG": "PRE-MIR7-2",
            "4-3074876-3074933-CAG": "HTT",
        })

    def test_skips_entries_without_diseases_or_without_TRExplorerV2(self):
        lookup = script.build_known_disease_associated_locus_ids_by_TRExplorerV2_id([
            make_known_locus("NoDiseases", "chr1:100-120", "1-100-120-AT", diseases=()),
            make_known_locus("NoTRExplorerV2", "chr1:200-220", None),
            make_known_locus("EmptyTRExplorerV2", "chr1:300-320", ""),
        ])
        self.assertEqual(lookup, {})

    def test_uses_the_first_entry_when_several_share_a_TRExplorerV2_id(self):
        lookup = script.build_known_disease_associated_locus_ids_by_TRExplorerV2_id([
            make_known_locus("SAMD12", "chr8:118366815-118366880", "8-118366812-118366918-AAAAT"),
            make_known_locus("SAMD12__8-118366812-118366918-AAAAT", "chr8:118366812-118366918",
                             "8-118366812-118366918-AAAAT"),
        ])
        self.assertEqual(lookup, {"8-118366812-118366918-AAAAT": "SAMD12"})


class TestMain(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.known_loci_path = os.path.join(self.temp_dir, "known_loci.json")
        self.catalog_path = os.path.join(self.temp_dir, "catalog.json.gz")
        self.output_path = os.path.join(self.temp_dir, "output.json.gz")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def run_main(self, known_loci, catalog_records):
        with open(self.known_loci_path, "wt") as f:
            json.dump(known_loci, f)
        with gzip.open(self.catalog_path, "wt") as f:
            json.dump(catalog_records, f)

        argv = ["script", "--known-disease-associated-loci", self.known_loci_path,
                "-o", self.output_path, self.catalog_path]
        with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print"):
            script.main()

        with gzip.open(self.output_path, "rt") as f:
            return json.load(f)

    def test_labels_only_loci_that_the_exact_boundary_match_left_unlabeled(self):
        known_loci = [
            make_known_locus("PRE-MIR7-2", "chr15:88569433-88569449", "15-88569433-88569452-TTTG"),
            make_known_locus("SAMD12", "chr8:118366815-118366880", "8-118366812-118366918-AAAAT"),
            make_known_locus("SAMD12__8-118366812-118366918-AAAAT", "chr8:118366812-118366918",
                             "8-118366812-118366918-AAAAT"),
        ]
        catalog_records = [
            {"LocusId": "15-88569433-88569452-TTTG", "KnownDiseaseAssociatedLocus": None, "NumRepeatsInReference": 4.75},
            {"LocusId": "8-118366812-118366918-AAAAT",
             "KnownDiseaseAssociatedLocus": "SAMD12__8-118366812-118366918-AAAAT"},
            {"LocusId": "1-100-120-AT", "KnownDiseaseAssociatedLocus": None},
            {"LocusId": "2-100-130-CAG"},
        ]

        output_records = self.run_main(known_loci, catalog_records)

        self.assertEqual([record["LocusId"] for record in output_records],
                         [record["LocusId"] for record in catalog_records])
        self.assertEqual([record.get("KnownDiseaseAssociatedLocus") for record in output_records], [
            "PRE-MIR7-2",
            "SAMD12__8-118366812-118366918-AAAAT",  # the exact-boundary label is kept
            None,
            None,
        ])
        self.assertEqual(output_records[0]["NumRepeatsInReference"], 4.75)

    def test_exits_when_no_known_locus_has_a_TRExplorerV2_field(self):
        known_loci = [make_known_locus("PRE-MIR7-2", "chr15:88569433-88569449", None)]
        with self.assertRaises(SystemExit), mock.patch("sys.stderr"):
            self.run_main(known_loci, [{"LocusId": "15-88569433-88569452-TTTG"}])


if __name__ == "__main__":
    unittest.main()
