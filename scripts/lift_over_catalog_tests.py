import os
import shutil
import tempfile
import unittest

import pysam

import lift_over_catalog as script

LEFT_FLANK = "ACGTTGCAAG"
RIGHT_FLANK = "TTGACCTAGGATCCA"

# Every source contig is 10 flanking bases, CAG x 5, then 15 flanking bases (40bp)
SOURCE_SEQUENCE = LEFT_FLANK + "CAGCAGCAGCAGCAG" + RIGHT_FLANK

# target chr1: the same flanks around a longer repeat (CAG x 7)
TARGET_CHR1 = LEFT_FLANK + "CAGCAGCAGCAGCAGCAGCAG" + RIGHT_FLANK
# target chr2: the reverse complement of SOURCE_SEQUENCE
TARGET_CHR2 = "TGGATCCTAGGTCAA" + "CTGCTGCTGCTGCTG" + "CTTGCAACGT"
# target chr3: the same flanks around a sequence that isn't a CAG repeat
TARGET_CHR3 = LEFT_FLANK + "CAGTTAGGCATCTAG" + RIGHT_FLANK
# target chr4: the same flanks around a repeat with Ns
TARGET_CHR4 = LEFT_FLANK + "CAGCAGNNNNNNCAG" + RIGHT_FLANK
# target chr5: the flanks without the repeat
TARGET_CHR5 = LEFT_FLANK + RIGHT_FLANK
# target chr6: the same flanks around a repeat with an extra leading G, so it starts at a different point of the motif
TARGET_CHR6 = LEFT_FLANK + "GCAGCAGCAGCAGCAG" + RIGHT_FLANK

# Like HTT: a CAG repeat followed by CAACAG, whose first two bases continue the CAG period. Target chr7 has 2 more CAGs.
SOURCE_CHR8 = LEFT_FLANK + "CAGCAGCAG" + "CAACAG" + RIGHT_FLANK
TARGET_CHR7 = LEFT_FLANK + "CAGCAGCAGCAGCAG" + "CAACAG" + RIGHT_FLANK
# Like FXN: a poly-A repeat followed by a GAA repeat. Target chr8 has 2 more GAAs.
SOURCE_CHR9 = LEFT_FLANK + "A" * 12 + "GAA" * 3 + "TTGACCTAG"
TARGET_CHR8 = LEFT_FLANK + "A" * 12 + "GAA" * 5 + "TTGACCTAG"
# Like NOP56: GGCCTG copies followed by CGCCTG copies
TARGET_CHR9 = LEFT_FLANK + "GGCCTG" * 3 + "CGCCTG" * 2 + RIGHT_FLANK

# The two bases before every source repeat (LEFT_FLANK's last two bases, AG) continue the CAG period.
# source chr1 => target chr1 with a 6bp target-side gap between source bases 14 and 15 (inside the repeat)
# source chr2 => reverse strand of target chr2
# source chr3 => target chr1, with a 1bp source-side gap at source base 7, just left of those two bases
# source chr4 => target chr3, source chr5 => target chr4
# source chr6 => target chr5, with the 15bp repeat as a source-side gap
# source chr7 => target chr6, with a 1bp target-side gap right after the left anchor
# source chr8 => target chr7, with the 6bp length difference placed after the C of CAACAG, inside target chr7's repeat
# source chr9 => target chr8, with the 6bp length difference placed inside the poly-A repeat
CHAIN_FILE_CONTENTS = """chain 1000 chr1 40 + 0 40 chr1 46 + 0 46 1
15	0	6
25

chain 1000 chr2 40 + 0 40 chr2 40 - 0 40 2
40

chain 1000 chr3 40 + 0 40 chr1 46 + 0 39 3
7	1	0
32

chain 1000 chr4 40 + 0 40 chr3 40 + 0 40 4
40

chain 1000 chr5 40 + 0 40 chr4 40 + 0 40 5
40

chain 1000 chr6 40 + 0 40 chr5 25 + 0 25 6
10	15	0
15

chain 1000 chr7 40 + 0 40 chr6 41 + 0 41 7
10	0	1
30

chain 1000 chr8 40 + 0 40 chr7 46 + 0 46 8
20	0	6
20

chain 1000 chr9 40 + 0 40 chr8 46 + 0 46 9
16	0	6
24

"""


def write_fasta(path, contigs):
    with open(path, "wt") as f:
        for name, sequence in contigs:
            f.write(f">{name}\n{sequence}\n")
    pysam.faidx(path)


def make_record(reference_region, locus_structure="(CAG)*", **extra_fields):
    return {"LocusId": "TEST", "ReferenceRegion": reference_region, "LocusStructure": locus_structure,
            "VariantType": "Repeat", **extra_fields}


class TestLiftOverCatalog(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        chain_path = os.path.join(self.temp_dir, "test.over.chain")
        with open(chain_path, "wt") as f:
            f.write(CHAIN_FILE_CONTENTS)
        self.chain_index = script.parse_chain_file(chain_path)

        source_fasta_path = os.path.join(self.temp_dir, "source.fa")
        write_fasta(source_fasta_path, [(f"chr{i}", SOURCE_SEQUENCE) for i in range(1, 8)] + [
            ("chr8", SOURCE_CHR8), ("chr9", SOURCE_CHR9)])
        target_fasta_path = os.path.join(self.temp_dir, "target.fa")
        write_fasta(target_fasta_path, [
            ("chr1", TARGET_CHR1), ("chr2", TARGET_CHR2), ("chr3", TARGET_CHR3), ("chr4", TARGET_CHR4),
            ("chr5", TARGET_CHR5), ("chr6", TARGET_CHR6), ("chr7", TARGET_CHR7), ("chr8", TARGET_CHR8),
            ("chr9", TARGET_CHR9)])
        self.source_fasta = pysam.FastaFile(source_fasta_path)
        self.target_fasta = pysam.FastaFile(target_fasta_path)

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def lift(self, record, keep_only_core_fields=False, max_anchor_search_distance=0):
        return script.lift_record(record, self.chain_index, self.source_fasta, self.target_fasta,
                                  max_purity_decrease=0.2, keep_only_core_fields=keep_only_core_fields,
                                  original_locus_id_field="LocusIdBeforeLiftover",
                                  max_anchor_search_distance=max_anchor_search_distance)

    def test_map_position(self):
        self.assertEqual(script.map_position(self.chain_index, "chr1", 9),
                         script.MappedPosition("chr1", 9, "+", "1"))
        self.assertEqual(script.map_position(self.chain_index, "chr1", 14).position_0based, 14)
        self.assertEqual(script.map_position(self.chain_index, "chr1", 15).position_0based, 21)
        self.assertEqual(script.map_position(self.chain_index, "chr2", 0),
                         script.MappedPosition("chr2", 39, "-", "2"))
        self.assertIsNone(script.map_position(self.chain_index, "chr3", 7))
        self.assertEqual(script.map_position(self.chain_index, "chr3", 8).position_0based, 7)
        self.assertIsNone(script.map_position(self.chain_index, "chr1", 40))
        self.assertIsNone(script.map_position(self.chain_index, "chrUn", 10))

    def test_lifts_repeat_whose_length_differs_between_references(self):
        lifted_record, reason = self.lift(make_record("chr1:10-25", SomeAnnotation=1))
        self.assertIsNone(reason)
        self.assertEqual(lifted_record["ReferenceRegion"], "chr1:10-31")
        self.assertEqual(self.target_fasta.fetch("chr1", 10, 31), "CAG" * 7)
        self.assertEqual(lifted_record["LocusStructure"], "(CAG)*")
        self.assertEqual(lifted_record["LocusIdBeforeLiftover"], "TEST")
        self.assertEqual(lifted_record["SomeAnnotation"], 1)

    def test_keep_only_core_fields(self):
        lifted_record, _ = self.lift(make_record("chr1:10-25", SomeAnnotation=1), keep_only_core_fields=True)
        self.assertNotIn("SomeAnnotation", lifted_record)
        self.assertEqual(set(lifted_record), {
            "LocusId", "ReferenceRegion", "LocusStructure", "VariantType", "LocusIdBeforeLiftover"})

    def test_reverse_strand_reverse_complements_and_rotates_the_motif(self):
        # The source repeat CAGCAGCAGCAGCA ends in a partial repeat. On the reverse strand it's TGCTGCTGCTGCTG,
        # which starts with a partial repeat, so the reverse complemented motif CTG is rotated to TGC.
        lifted_record, reason = self.lift(make_record("chr2:10-24"))
        self.assertIsNone(reason)
        self.assertEqual(lifted_record["ReferenceRegion"], "chr2:16-30")
        self.assertEqual(self.target_fasta.fetch("chr2", 16, 30), "TGCTGCTGCTGCTG")
        self.assertEqual(lifted_record["LocusStructure"], "(TGC)*")

    def test_rotates_motif_when_the_repeat_starts_at_a_different_point_of_the_motif(self):
        # GCAGCAGCAGCAGCAG has a purity near 0 against CAG repeated from its first base, but is a pure GCA repeat
        lifted_record, reason = self.lift(make_record("chr7:10-25"))
        self.assertIsNone(reason)
        self.assertEqual(lifted_record["ReferenceRegion"], "chr6:10-26")
        self.assertEqual(self.target_fasta.fetch("chr6", 10, 26), "GCAGCAGCAGCAGCAG")
        self.assertEqual(lifted_record["LocusStructure"], "(GCA)*")

    def test_rotates_each_motif_of_compound_locus(self):
        # the first repeat absorbs target chr6's extra leading G, so it reads GCAGCAG and its motif becomes GCA
        lifted_record, reason = self.lift(make_record(
            ["chr7:10-16", "chr7:16-25"], "(CAG)*AC(CAG)*", VariantType=["Repeat", "Repeat"]))
        self.assertIsNone(reason)
        self.assertEqual(lifted_record["ReferenceRegion"], ["chr6:10-17", "chr6:17-26"])
        self.assertEqual(lifted_record["LocusStructure"], "(GCA)*AC(CAG)*")

    def test_drops_locus_whose_anchor_is_not_covered_by_the_chain(self):
        # the left anchor skips source bases 9 and 8, which continue the CAG period, so it's source base 7, a gap
        self.assertEqual(self.lift(make_record("chr3:10-25")), (None, "anchor not covered by the chain file"))

    def test_moves_anchor_in_a_chain_gap_outward_to_the_nearest_covered_base(self):
        # the left anchor moves from source base 7 to base 6, which maps to target chr1 base 6. Source base 8 maps to
        # target base 7, so the lifted repeat starts at target base 7 plus the 2 skipped bases. This mapping of chr3
        # onto chr1 is arbitrary, so the target sequence there starts at a different point of the motif, and since
        # target chr1's repeat is longer, the lifted repeat also extends over the 2 whole GCA copies to its right.
        lifted_record, reason = self.lift(make_record("chr3:10-25"), max_anchor_search_distance=1)
        self.assertIsNone(reason)
        self.assertEqual(lifted_record["ReferenceRegion"], "chr1:9-30")
        self.assertEqual(self.target_fasta.fetch("chr1", 9, 30), "GCA" * 7)
        self.assertEqual(lifted_record["LocusStructure"], "(GCA)*")

    def test_skips_flanking_bases_that_continue_the_repeat(self):
        # Mapping source base 19 (the C of CAACAG) would put the right anchor inside target chr7's longer repeat and
        # truncate it to 3 CAGs. Skipping the CA that continues the period puts the anchor at source base 21, which
        # maps to target base 27, and trimming the 2 skipped bases gives the whole repeat of 5 CAGs.
        self.assertEqual(script.map_position(self.chain_index, "chr8", 19).position_0based, 19)
        lifted_record, reason = self.lift(make_record("chr8:10-19"))
        self.assertIsNone(reason)
        self.assertEqual(lifted_record["ReferenceRegion"], "chr7:10-25")
        self.assertEqual(self.target_fasta.fetch("chr7", 10, 25), "CAG" * 5)

    def test_moves_whole_copies_across_the_boundary_between_adjacent_repeats(self):
        # The chain maps the poly-A repeat's right anchor (source base 22, the first G) to target base 28, so the
        # poly-A repeat first absorbs 2 GAA copies, which are then moved to the GAA repeat.
        lifted_record, reason = self.lift(make_record(
            ["chr9:10-22", "chr9:22-31"], "(A)*(GAA)*", VariantType=["Repeat", "Repeat"]))
        self.assertIsNone(reason)
        self.assertEqual(lifted_record["ReferenceRegion"], ["chr8:10-22", "chr8:22-37"])
        self.assertEqual(lifted_record["LocusStructure"], "(A)*(GAA)*")

    def test_move_boundary_between_adjacent_repeats(self):
        # chr8: 12 As then 5 GAAs from base 22. Moves the 2 GAA copies at the end of the poly-A repeat to the right.
        self.assertEqual(script.move_boundary_between_adjacent_repeats(
            self.target_fasta, "chr8", 10, 28, 37, "A", "GAA"), 22)
        # chr9: 3 GGCCTGs from base 10, then 2 CGCCTGs. Moves the GGCCTG copy at the start of the right repeat left.
        self.assertEqual(script.move_boundary_between_adjacent_repeats(
            self.target_fasta, "chr9", 10, 22, 40, "GGCCTG", "CGCCTG"), 28)
        # already at the right place
        self.assertEqual(script.move_boundary_between_adjacent_repeats(
            self.target_fasta, "chr9", 10, 28, 40, "GGCCTG", "CGCCTG"), 28)
        # off by part of a motif: the right repeat's GAA motif is given in the phase that matches from base 24
        self.assertEqual(script.move_boundary_between_adjacent_repeats(
            self.target_fasta, "chr8", 10, 24, 37, "A", "AGA"), 22)

    def test_count_flanking_bases_that_continue_the_repeat(self):
        # before chr8's repeat: ...CAAG|CAGCAGCAG, so G and A continue it. After it: CAGCAGCAG|CAACAG, so C and A do.
        self.assertEqual(script.count_flanking_bases_that_continue_the_repeat(
            self.source_fasta, "chr8", 10, 19, "CAG", -1, 100), 2)
        self.assertEqual(script.count_flanking_bases_that_continue_the_repeat(
            self.source_fasta, "chr8", 10, 19, "CAG", 1, 100), 2)
        self.assertEqual(script.count_flanking_bases_that_continue_the_repeat(
            self.source_fasta, "chr8", 10, 19, "CAG", 1, 1), 1)

    def test_map_nearest_base_covered_by_chain(self):
        self.assertEqual(script.map_nearest_base_covered_by_chain(self.chain_index, "chr3", 7, -1, 1).position_0based, 6)
        self.assertEqual(script.map_nearest_base_covered_by_chain(self.chain_index, "chr3", 7, 1, 1).position_0based, 7)
        self.assertIsNone(script.map_nearest_base_covered_by_chain(self.chain_index, "chr3", 7, -1, 0))
        self.assertIsNone(script.map_nearest_base_covered_by_chain(self.chain_index, "chr3", 7, -1, 1, chain_id="1"))

    def test_find_loci_that_collide_with_loci_from_elsewhere(self):
        source_intervals = [
            ("chr1", 100, 120),    # 0: adjacent to locus 1 in the source
            ("chr1", 120, 140),    # 1
            ("chr1", 500, 520),    # 2: on another chromosome than locus 3 in the source
            ("chr9", 100, 120),    # 3
            ("chr1", 700, 740),    # 4: overlaps locus 5 in the source
            ("chr1", 710, 730),    # 5
            ("chr1", 2000, 2010),  # 6: 2 bases from locus 7 in the source
            ("chr1", 2012, 2020),  # 7
            ("chr1", 5000, 5010),  # 8: 4kb from locus 9 in the source
            ("chr1", 9000, 9010),  # 9
        ]
        target_intervals = [
            ("chr1", 100, 126),  # 0 and 1 overlap because both absorbed a target-side insertion between them
            ("chr1", 120, 146),
            ("chr2", 300, 320),  # 2 and 3 land on the same repeat
            ("chr2", 305, 320),
            ("chr1", 900, 940),  # 4 and 5 overlap in both references
            ("chr1", 910, 930),
            ("chr3", 50, 70),    # 6 and 7 land on one repeat, since the base between them matches it in the target
            ("chr3", 55, 70),
            ("chr4", 50, 60),    # 8 and 9 land on the same repeat
            ("chr4", 50, 60),
        ]
        self.assertEqual(script.find_loci_that_collide_with_loci_from_elsewhere(source_intervals, target_intervals),
                         {2, 3, 8, 9})

    def test_drops_locus_whose_purity_decreases(self):
        self.assertEqual(self.lift(make_record("chr4:10-25")), (None, "purity decreased by more than 0.2"))

    def test_drops_locus_with_Ns_in_the_new_reference(self):
        self.assertEqual(self.lift(make_record("chr5:10-25")), (None, "new reference sequence contains Ns"))

    def test_drops_locus_whose_repeat_is_absent_from_the_new_reference(self):
        self.assertEqual(self.lift(make_record("chr6:10-25")), (None, "repeat absent from the new reference"))

    def test_lifts_compound_locus_and_main_reference_region(self):
        # the 6bp target-side gap falls inside the first repeat, so the first repeat absorbs it
        lifted_record, reason = self.lift(make_record(
            ["chr1:10-16", "chr1:16-25"], "(CAG)*(CAG)*", MainReferenceRegion="chr1:16-25",
            VariantType=["Repeat", "Repeat"], VariantId=["TEST", "TEST_2"], Diseases=[{"Symbol": "X"}]))
        self.assertIsNone(reason)
        self.assertEqual(lifted_record["ReferenceRegion"], ["chr1:10-22", "chr1:22-31"])
        self.assertEqual(lifted_record["MainReferenceRegion"], "chr1:22-31")
        self.assertEqual(lifted_record["VariantId"], ["TEST", "TEST_2"])
        self.assertEqual(lifted_record["Diseases"], [{"Symbol": "X"}])

    def test_drops_compound_locus_on_reverse_strand(self):
        self.assertEqual(
            self.lift(make_record(["chr2:10-16", "chr2:16-25"], "(CAG)*(CAG)*", VariantType=["Repeat", "Repeat"])),
            (None, "locus with more than one repeat maps to the reverse strand"))

    def test_main_writes_sorted_catalog_and_dropped_loci(self):
        import json
        import sys
        from unittest import mock

        catalog_path = os.path.join(self.temp_dir, "catalog.json")
        with open(catalog_path, "wt") as f:
            json.dump([make_record("chr2:10-24"), make_record("chr4:10-25"), make_record("chr1:10-25")], f)
        chain_path = os.path.join(self.temp_dir, "test.over.chain")
        output_path = os.path.join(self.temp_dir, "output.json")
        dropped_loci_path = os.path.join(self.temp_dir, "dropped.tsv")
        argv = ["script", "--chain-file", chain_path,
                "--source-reference-fasta", os.path.join(self.temp_dir, "source.fa"),
                "--target-reference-fasta", os.path.join(self.temp_dir, "target.fa"),
                "--output-dropped-loci-tsv", dropped_loci_path, "-o", output_path, catalog_path]
        with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print"):
            script.main()

        with open(output_path) as f:
            self.assertEqual([r["ReferenceRegion"] for r in json.load(f)], ["chr1:10-31", "chr2:16-30"])
        with open(dropped_loci_path) as f:
            self.assertEqual(f.read().splitlines(), [
                "LocusId\tReferenceRegion\tReason", "TEST\tchr4:10-25\tpurity decreased by more than 0.2"])

    def test_main_keeps_only_the_first_of_loci_with_the_same_lifted_definition(self):
        import json
        import sys
        from unittest import mock

        # two loci with the same definition both lift to the whole CAG x 7 repeat at target chr1:10-31
        catalog_path = os.path.join(self.temp_dir, "catalog.json")
        with open(catalog_path, "wt") as f:
            json.dump([make_record("chr1:10-25"), dict(make_record("chr1:10-25"), LocusId="TEST2")], f)
        output_path = os.path.join(self.temp_dir, "output.json")
        dropped_loci_path = os.path.join(self.temp_dir, "dropped.tsv")
        argv = ["script", "--chain-file", os.path.join(self.temp_dir, "test.over.chain"),
                "--source-reference-fasta", os.path.join(self.temp_dir, "source.fa"),
                "--target-reference-fasta", os.path.join(self.temp_dir, "target.fa"),
                "--output-dropped-loci-tsv", dropped_loci_path, "-o", output_path, catalog_path]
        with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print"):
            script.main()

        with open(output_path) as f:
            self.assertEqual([(r["LocusId"], r["ReferenceRegion"]) for r in json.load(f)], [("TEST", "chr1:10-31")])
        with open(dropped_loci_path) as f:
            self.assertEqual(f.read().splitlines()[1:], [
                "TEST2\tchr1:10-25\thas the same definition in the new reference as another lifted locus"])

    def test_update_coordinates_in_locus_id(self):
        self.assertEqual(script.update_coordinates_in_locus_id(dict(
            make_record("chr12:132111361-132111451"), LocusId="EP400__12-132062524-132062611-CAG")),
            "EP400__12-132111361-132111451-CAG")
        # the motif comes from the lifted LocusStructure, which can be a different rotation
        self.assertEqual(script.update_coordinates_in_locus_id(dict(
            make_record("chr17:17755019-17755053", "(ATTTT)*"), LocusId="RAI1__17-17808416-17808460-TTTTA")),
            "RAI1__17-17755019-17755053-ATTTT")
        self.assertEqual(script.update_coordinates_in_locus_id(dict(
            make_record("chr6:45257567-45257618", "(GCN)*"), LocusId="RUNX2__6-45422750-45422801-GCN")),
            "RUNX2__6-45257567-45257618-GCN")
        # LocusIds without a coordinate-based id after "__" are left as they are
        for locus_id in "HTT", "1-10-25-CAG", "HOXA13_1":
            self.assertEqual(script.update_coordinates_in_locus_id(dict(make_record("chr1:10-31"), LocusId=locus_id)),
                             locus_id)

    def test_main_updates_coordinates_in_locus_ids(self):
        import json
        import sys
        from unittest import mock

        catalog_path = os.path.join(self.temp_dir, "catalog.json")
        with open(catalog_path, "wt") as f:
            json.dump([make_record("chr1:10-25"), dict(make_record("chr2:10-24"), LocusId="TEST__1-10-25-CAG")], f)
        output_path = os.path.join(self.temp_dir, "output.json")
        argv = ["script", "--chain-file", os.path.join(self.temp_dir, "test.over.chain"),
                "--source-reference-fasta", os.path.join(self.temp_dir, "source.fa"),
                "--target-reference-fasta", os.path.join(self.temp_dir, "target.fa"),
                "--update-coordinates-in-locus-ids", "-o", output_path, catalog_path]
        with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print"):
            script.main()

        with open(output_path) as f:
            self.assertEqual([(r["LocusId"], r["LocusIdBeforeLiftover"]) for r in json.load(f)], [
                ("TEST", "TEST"), ("TEST__2-16-30-TGC", "TEST__1-10-25-CAG")])


if __name__ == "__main__":
    unittest.main()
