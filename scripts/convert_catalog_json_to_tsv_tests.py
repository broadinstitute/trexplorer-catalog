import gzip
import json
import os
import tempfile
import unittest
from unittest import mock

import pandas as pd

import convert_catalog_json_to_tsv
from convert_catalog_json_to_tsv import CORE_COLUMNS, DROP_COLUMNS, write_catalog_tsv

# Records whose columns mix types in the ways that change how pandas writes a column, eg. integers that some records
# lack (written as 12.0), integers mixed with strings (written as 12), and booleans that some records lack. The
# OnlyInLastRecord column appears after the first chunk when the chunk size is 2.
RECORDS = [
    {"LocusId": "1-100-110-A", "ReferenceRegion": "chr1:100-110", "VariantType": "Repeat", "IntsSomeMissing": 3,
     "IntsAndFloats": 1, "IntsAndStrings": 7, "AllInts": 1, "AllBools": True, "BoolsSomeNone": True,
     "Lists": ["a", "b"], "Floats": 0.1, "AllNone": None, "StringWithTab": "x\ty"},
    {"LocusId": "1-200-230-CAG", "ReferenceRegion": "chr1:200-230", "VariantType": "Repeat",
     "IntsAndFloats": 2.5, "IntsAndStrings": "PAR1", "AllInts": 2, "AllBools": False, "BoolsSomeNone": None,
     "Lists": [], "Floats": 1e-05, "AllNone": None, "StringWithTab": "plain"},
    {"LocusId": "2-50-60-AT", "ReferenceRegion": "chr2:50-60", "VariantType": "Repeat", "IntsSomeMissing": 12,
     "IntsAndFloats": 3, "IntsAndStrings": 9, "AllInts": 3, "AllBools": True, "BoolsSomeNone": False,
     "Lists": ["c"], "Floats": 0.30000000000000004, "AllNone": None, "StringWithTab": 'has "quotes"',
     "OnlyInLastRecord": 5},
]


def write_tsv_by_loading_the_whole_catalog(catalog_json_path, output_tsv_path, omit_core_columns_missing_from_catalog):
    """The converter before it read records one at a time, which the streaming version must match."""
    with gzip.open(catalog_json_path, "rt") as f:
        df = pd.DataFrame(json.load(f))
    if omit_core_columns_missing_from_catalog:
        core_columns = [c for c in CORE_COLUMNS if c in df.columns]
    else:
        core_columns = CORE_COLUMNS
        for c in set(core_columns) - set(df.columns):
            df[c] = None
    df = df[core_columns + [c for c in df.columns if c not in (core_columns + DROP_COLUMNS)]]
    df.to_csv(output_tsv_path, sep="\t", index=False)


class ConvertCatalogJsonToTsvTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.catalog_json_path = os.path.join(self.temp_dir.name, "catalog.json.gz")
        with gzip.open(self.catalog_json_path, "wt") as f:
            json.dump(RECORDS, f)

    def tearDown(self):
        self.temp_dir.cleanup()

    def assert_output_matches_loading_the_whole_catalog(self, omit_core_columns_missing_from_catalog):
        expected_tsv_path = os.path.join(self.temp_dir.name, "expected.tsv.gz")
        write_tsv_by_loading_the_whole_catalog(
            self.catalog_json_path, expected_tsv_path, omit_core_columns_missing_from_catalog)

        output_tsv_path = os.path.join(self.temp_dir.name, "output.tsv.gz")
        with mock.patch.object(convert_catalog_json_to_tsv, "RECORDS_PER_CHUNK", 2):
            num_rows, _ = write_catalog_tsv(
                self.catalog_json_path, output_tsv_path, omit_core_columns_missing_from_catalog)

        self.assertEqual(num_rows, len(RECORDS))
        with gzip.open(expected_tsv_path, "rt") as f:
            expected_tsv = f.read()
        with gzip.open(output_tsv_path, "rt") as f:
            self.assertEqual(f.read(), expected_tsv)
        return expected_tsv

    def test_output_matches_loading_the_whole_catalog(self):
        expected_tsv = self.assert_output_matches_loading_the_whole_catalog(False)
        # check that the records exercise the cases they're meant to
        header, row1, row2, row3 = expected_tsv.splitlines()
        self.assertIn("IntsSomeMissing", header)
        self.assertIn("12.0", row3)
        self.assertNotIn("VariantType", header)
        self.assertTrue(header.startswith("\t".join(CORE_COLUMNS)))

    def test_output_matches_loading_the_whole_catalog_without_missing_core_columns(self):
        expected_tsv = self.assert_output_matches_loading_the_whole_catalog(True)
        self.assertTrue(expected_tsv.startswith("LocusId\tReferenceRegion\tIntsSomeMissing"))


if __name__ == "__main__":
    unittest.main()
