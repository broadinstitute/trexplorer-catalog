import gzip
import os
import shutil
import simplejson as json
import sys
import tempfile
import unittest
from unittest import mock

import remove_loci_near_chromosome_ends as script


class TestRemoveLociNearChromosomeEnds(unittest.TestCase):

    def test_is_near_chromosome_end(self):
        chromosome_sizes = {"chr1": 10_000}
        self.assertTrue(script.is_near_chromosome_end({"ReferenceRegion": "chr1:0-60"}, chromosome_sizes, 1000))
        self.assertTrue(script.is_near_chromosome_end({"ReferenceRegion": "chr1:999-1060"}, chromosome_sizes, 1000))
        self.assertFalse(script.is_near_chromosome_end({"ReferenceRegion": "chr1:1000-1060"}, chromosome_sizes, 1000))
        self.assertFalse(script.is_near_chromosome_end({"ReferenceRegion": "chr1:8900-9000"}, chromosome_sizes, 1000))
        self.assertTrue(script.is_near_chromosome_end({"ReferenceRegion": "chr1:8900-9001"}, chromosome_sizes, 1000))
        self.assertTrue(script.is_near_chromosome_end(
            {"ReferenceRegion": ["chr1:5000-5010", "chr1:9500-9600"]}, chromosome_sizes, 1000))

    def test_main(self):
        temp_dir = tempfile.mkdtemp()
        try:
            fasta_path = os.path.join(temp_dir, "ref.fa")
            with open(f"{fasta_path}.fai", "wt") as f:
                f.write("chr1\t10000\t6\t60\t61\n")
            catalog_path = os.path.join(temp_dir, "catalog.json.gz")
            with gzip.open(catalog_path, "wt") as f:
                json.dump([
                    {"LocusId": "1-0-60-TTAGGG", "ReferenceRegion": "chr1:0-60", "NumRepeatsInReference": 10.5},
                    {"LocusId": "1-5000-5010-A", "ReferenceRegion": "chr1:5000-5010", "NumRepeatsInReference": 10.5},
                ], f)
            output_path = os.path.join(temp_dir, "output.json.gz")
            argv = ["script", "-R", fasta_path, "-o", output_path, catalog_path]
            with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print"):
                script.main()
            with gzip.open(output_path, "rt") as f:
                self.assertEqual(json.load(f), [
                    {"LocusId": "1-5000-5010-A", "ReferenceRegion": "chr1:5000-5010", "NumRepeatsInReference": 10.5}])
        finally:
            shutil.rmtree(temp_dir)


if __name__ == "__main__":
    unittest.main()
