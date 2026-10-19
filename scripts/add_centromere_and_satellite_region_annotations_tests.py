import gzip
import os
import shutil
import simplejson as json
import sys
import tempfile
import unittest
from unittest import mock

import add_centromere_and_satellite_region_annotations as script

CENSAT_BED_LINES = [
    'track name="hg38.formatted.RM" visibility=2 itemRgb="On"',
    "chr1\t100\t200\tct\t0\t.\t100\t200\t224,224,224",
    "chr1\t200\t300\tactive_hor(S1C1/5/19H1L)\t0\t.\t200\t300\t255,0,0",
    "chr1\t300\t400\tGAP,HSat3\t0\t.\t300\t400\t0,0,255",
    "chr2\t50\t60\tbSat(COMP-subunit_LSAU-BSAT_rnd-1_family-1)\t0\t.\t50\t60\t250,153,255",
]


class TestAddCentromereAndSatelliteRegionAnnotations(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.bed_path = os.path.join(self.temp_dir, "censat.bed")
        with open(self.bed_path, "wt") as f:
            f.write("\n".join(CENSAT_BED_LINES) + "\n")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_parse_region_class(self):
        self.assertEqual(script.parse_region_class("active_hor(S1C1/5/19H1L)"), "active_hor")
        self.assertEqual(script.parse_region_class("bSat(COMP-subunit_LSAU-BSAT_rnd-1_family-1)"), "bSat")
        self.assertEqual(script.parse_region_class("GAP,HSat3"), "GAP,HSat3")
        self.assertEqual(script.parse_region_class("ct"), "ct")

    def test_get_class_of_region_with_largest_overlap(self):
        regions_by_chrom = script.read_centromere_and_satellite_regions(self.bed_path)
        self.assertEqual(sum(len(t) for t in regions_by_chrom.values()), 4)

        get_class = script.get_class_of_region_with_largest_overlap
        self.assertEqual(get_class(regions_by_chrom, "chr1", 120, 130), "ct")
        self.assertEqual(get_class(regions_by_chrom, "chr1", 190, 230), "active_hor")  # 30bp vs. 10bp of ct
        self.assertEqual(get_class(regions_by_chrom, "chr1", 280, 320), "active_hor")  # tie goes to the leftmost
        self.assertEqual(get_class(regions_by_chrom, "chr1", 350, 360), "GAP,HSat3")
        self.assertEqual(get_class(regions_by_chrom, "chr2", 55, 70), "bSat")
        self.assertIsNone(get_class(regions_by_chrom, "chr1", 400, 410))  # BED intervals are half-open
        self.assertIsNone(get_class(regions_by_chrom, "chr1", 90, 100))
        self.assertIsNone(get_class(regions_by_chrom, "chr3", 100, 200))

    def test_main(self):
        catalog_path = os.path.join(self.temp_dir, "catalog.json.gz")
        output_path = os.path.join(self.temp_dir, "output.json.gz")
        catalog_records = [
            {"LocusId": "1-190-230-AT", "ReferenceRegion": "chr1:190-230", "LocusStructure": "(AT)*",
             "VariantType": "Repeat", "Source": "SomeSource"},
            {"LocusId": "1-500-530-CAG", "ReferenceRegion": "chr1:500-530", "LocusStructure": "(CAG)*",
             "VariantType": "Repeat"},
        ]
        with gzip.open(catalog_path, "wt") as f:
            json.dump(catalog_records, f)

        argv = ["script", "--centromere-and-satellite-regions-bed", self.bed_path, "-o", output_path, catalog_path]
        with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print"):
            script.main()

        with gzip.open(output_path, "rt") as f:
            output_records = json.load(f)

        self.assertEqual(output_records, [
            {**catalog_records[0], "CentromereAndSatelliteRegion": "active_hor"},
            {**catalog_records[1], "CentromereAndSatelliteRegion": None},
        ])


if __name__ == "__main__":
    unittest.main()
