import gzip
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

import lift_over_bed_regions as script
from lift_over_catalog import parse_chain_file

# Three chains from a 1,000 bp chr1 of the source reference to the target reference:
# - chain 1 maps chr1:100-150 to chr1:500-550 and chr1:160-300 to chr1:570-710 of the target (between them, the source
#   has 10 unaligned bases and the target has 20)
# - chain 2 maps chr1:400-450 to the reverse strand of the target's chr1, at chr1:850-900
# - chain 3 maps chr1:600-650 to the target's chr2, at chr2:0-50
CHAIN_FILE_CONTENTS = """chain 1000 chr1 1000 + 100 300 chr1 1000 + 500 710 1
50\t10\t20
140

chain 500 chr1 1000 + 400 450 chr1 1000 - 100 150 2
50

chain 300 chr1 1000 + 600 650 chr2 1000 + 0 50 3
50

"""

# The chains back from the target to the source reference:
# - chain 1 maps the target's chr1:500-550 and chr1:570-710 back to chr1:100-150 and chr1:160-300, the reverse of
#   chain 1 above
# - chain 2 maps the target's chr1:850-900 to the source's chr1:950-1000 rather than back to chr1:400-450, as for the
#   paralogous copy of a satellite array
REVERSE_CHAIN_FILE_CONTENTS = """chain 1000 chr1 1000 + 500 710 chr1 1000 + 100 300 1
50\t20\t10
140

chain 500 chr1 1000 + 850 900 chr1 1000 + 950 1000 2
50

"""


class TestLiftOverBedRegions(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.chain_file_path = os.path.join(self.temp_dir, "test.over.chain")
        with open(self.chain_file_path, "wt") as f:
            f.write(CHAIN_FILE_CONTENTS)
        self.reverse_chain_file_path = os.path.join(self.temp_dir, "test.reverse.over.chain")
        with open(self.reverse_chain_file_path, "wt") as f:
            f.write(REVERSE_CHAIN_FILE_CONTENTS)
        self.chain_index = parse_chain_file(self.chain_file_path)
        self.reverse_chain_index = parse_chain_file(self.reverse_chain_file_path)

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_lift_region(self):
        def lift(chrom, start_0based, end_1based, max_chain_gap_to_fill=20):
            return script.lift_region(self.chain_index, chrom, start_0based, end_1based, max_chain_gap_to_fill)

        # inside one block
        self.assertEqual(lift("chr1", 110, 130), [("chr1", 510, 530)])
        # across chain 1's gap: one region that also covers the target's 20 bases in the gap
        self.assertEqual(lift("chr1", 120, 200), [("chr1", 520, 610)])
        # the same region when the gap is larger than the limit: the two aligned pieces
        self.assertEqual(lift("chr1", 120, 200, max_chain_gap_to_fill=19), [("chr1", 520, 550), ("chr1", 570, 610)])
        # a region that starts in chain 1's gap
        self.assertEqual(lift("chr1", 152, 170), [("chr1", 570, 580)])
        # a reverse-strand chain
        self.assertEqual(lift("chr1", 400, 420), [("chr1", 880, 900)])
        # a region that two chains align lifts to one region per chain
        self.assertEqual(lift("chr1", 280, 420), [("chr1", 690, 710), ("chr1", 880, 900)])
        # a chain that aligns the region to another chromosome is ignored
        self.assertEqual(lift("chr1", 600, 650), [])
        self.assertEqual(lift("chr1", 420, 650), [("chr1", 850, 880)])
        # a region that no chain aligns, and a chromosome that isn't in the chain file
        self.assertEqual(lift("chr1", 310, 390), [])
        self.assertEqual(lift("chr1", 151, 159), [])
        self.assertEqual(lift("chr2", 100, 200), [])

    def test_lift_bed_regions_keeps_pieces_that_map_back_into_any_of_the_regions(self):
        # chr1:400-450 lifts to the target's chr1:850-900, which maps back to the source's chr1:950-1000. That piece is
        # kept only when the input regions also cover chr1:950-1000.
        lifted_regions, counters = script.lift_bed_regions(
            self.chain_index, self.reverse_chain_index, ["chr1\t400\t450\tHSat3\n"], max_chain_gap_to_fill=20)
        self.assertEqual(lifted_regions, [])
        self.assertEqual(counters["lifted pieces dropped since less than half of their bases map back into the regions"], 1)

        lifted_regions, _ = script.lift_bed_regions(
            self.chain_index, self.reverse_chain_index, ["chr1\t400\t450\tHSat3\n", "chr1\t940\t1000\tHSat2\n"],
            max_chain_gap_to_fill=20)
        self.assertEqual(lifted_regions, [("chr1", 850, 900, "HSat3")])

    def test_main(self):
        bed_path = os.path.join(self.temp_dir, "regions.bed")
        output_path = os.path.join(self.temp_dir, "lifted.bed.gz")
        with open(bed_path, "wt") as f:
            f.write('track name="test"\n')
            f.write("chr1\t280\t420\tHSat3\t0\t.\t280\t420\t0,0,0\n")
            f.write("chr1\t310\t390\tGAP\t0\t.\t310\t390\t0,0,0\n")
            f.write("chr1\t600\t650\tct\t0\t.\t600\t650\t0,0,0\n")
            f.write("chr1\t110\t130\tactive_hor(S1C1/5/19H1L)\t0\t.\t110\t130\t0,0,0\n")

        argv = ["script", "--chain-file", self.chain_file_path, "--reverse-chain-file", self.reverse_chain_file_path,
                "-o", output_path, bed_path]
        with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print"):
            script.main()

        # HSat3's piece at chr1:880-900 maps back to chr1:980-1000, outside the regions, so it's dropped
        with gzip.open(output_path, "rt") as f:
            self.assertEqual(f.read(), "chr1\t510\t530\tactive_hor(S1C1/5/19H1L)\n"
                                       "chr1\t690\t710\tHSat3\n")


if __name__ == "__main__":
    unittest.main()
