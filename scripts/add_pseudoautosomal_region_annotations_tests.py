import gzip
import os
import shutil
import simplejson as json
import sys
import tempfile
import unittest
from unittest import mock

import add_pseudoautosomal_region_annotations as script


class TestAddPseudoautosomalRegionAnnotations(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_get_pseudoautosomal_region_hg38(self):
        hg38_regions = script.PSEUDOAUTOSOMAL_REGIONS_BY_GENOME_VERSION["hg38"]
        get_region = script.get_pseudoautosomal_region
        self.assertEqual(get_region(hg38_regions, "chrX", 10000, 10030), "PAR1")     # starts at PAR1's first base
        self.assertEqual(get_region(hg38_regions, "chrX", 2781449, 2781479), "PAR1")  # ends at PAR1's last base
        self.assertEqual(get_region(hg38_regions, "X", 500000, 500030), "PAR1")       # no chr prefix
        self.assertEqual(get_region(hg38_regions, "chrX", 155701382, 155701400), "PAR2")
        self.assertEqual(get_region(hg38_regions, "chrY", 56887902, 56887920), "PAR2")
        self.assertIsNone(get_region(hg38_regions, "chrX", 9990, 10030))         # crosses PAR1's start
        self.assertIsNone(get_region(hg38_regions, "chrX", 2781470, 2781490))    # crosses PAR1's end
        self.assertIsNone(get_region(hg38_regions, "chrX", 2781479, 2781500))    # the end is exclusive
        self.assertIsNone(get_region(hg38_regions, "chrX", 155701370, 155701400))
        self.assertIsNone(get_region(hg38_regions, "chrX", 100000000, 100000030))
        self.assertIsNone(get_region(hg38_regions, "chr1", 500000, 500030))

    def test_get_pseudoautosomal_region_chm13(self):
        chm13_regions = script.PSEUDOAUTOSOMAL_REGIONS_BY_GENOME_VERSION["chm13"]
        get_region = script.get_pseudoautosomal_region
        self.assertEqual(get_region(chm13_regions, "chrX", 0, 30), "PAR1")
        self.assertEqual(get_region(chm13_regions, "chrX", 2394400, 2394410), "PAR1")
        self.assertIsNone(get_region(chm13_regions, "chrX", 2394400, 2394411))
        self.assertIsNone(get_region(chm13_regions, "chrX", 153925822, 153925839))  # crosses PAR2's start
        self.assertEqual(get_region(chm13_regions, "chrX", 153925834, 153925850), "PAR2")
        self.assertEqual(get_region(chm13_regions, "chrY", 62460000, 62460029), "PAR2")
        # inside hg38's chrX PAR1, but not CHM13's
        self.assertIsNone(get_region(chm13_regions, "chrX", 2700000, 2700030))

    def test_main(self):
        catalog_path = os.path.join(self.temp_dir, "catalog.json.gz")
        output_path = os.path.join(self.temp_dir, "output.json.gz")
        catalog_records = [
            {"LocusId": "X-500000-500030-CAG", "ReferenceRegion": "chrX:500000-500030", "LocusStructure": "(CAG)*",
             "VariantType": "Repeat", "Source": "SomeSource"},
            {"LocusId": "X-155701390-155701410-AT", "ReferenceRegion": "chrX:155701390-155701410",
             "LocusStructure": "(AT)*", "VariantType": "Repeat"},
            {"LocusId": "X-3000000-3000030-CAG", "ReferenceRegion": "chrX:3000000-3000030", "LocusStructure": "(CAG)*",
             "VariantType": "Repeat", "PseudoautosomalRegion": "PAR1"},   # a stale value from an earlier run
            {"LocusId": "1-500000-500030-CAG", "ReferenceRegion": "chr1:500000-500030", "LocusStructure": "(CAG)*",
             "VariantType": "Repeat"},
        ]
        with gzip.open(catalog_path, "wt") as f:
            json.dump(catalog_records, f)

        argv = ["script", "--genome-version", "hg38", "-o", output_path, catalog_path]
        with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print"):
            script.main()

        with gzip.open(output_path, "rt") as f:
            output_records = json.load(f)

        self.assertEqual(output_records, [
            {**catalog_records[0], "PseudoautosomalRegion": "PAR1"},
            {**catalog_records[1], "PseudoautosomalRegion": "PAR2"},
            {**catalog_records[2], "PseudoautosomalRegion": None},
            {**catalog_records[3], "PseudoautosomalRegion": None},
        ])


if __name__ == "__main__":
    unittest.main()
