import gzip
import os
import shutil
import simplejson as json
import sys
import tempfile
import unittest
from unittest import mock

import add_segdup_max_sequence_identity_annotations as script


def build_ucsc_genomic_super_dups_line(chrom, start_0based, end_1based, frac_match, other_copy=("chr2", 5000, 6000)):
    """Returns a 30-column line of UCSC's genomicSuperDups table with the given values and placeholders elsewhere."""
    other_chrom, other_start_0based, other_end_1based = other_copy
    fields = (["585", chrom, str(start_0based), str(end_1based), "name", "0", "+",
               other_chrom, str(other_start_0based), str(other_end_1based)]
              + ["0"] * 16 + [str(frac_match), "0.9", "0.01", "0.01"])
    return "\t".join(fields)


def build_vollger_2022_sd_bed_line(chrom, start_0based, end_1based, frac_match, other_copy=("chr2", 5000, 6000)):
    """Returns a 44-column line of a Vollger et al. 2022 SD BED file with the given values and placeholders elsewhere."""
    other_chrom, other_start_0based, other_end_1based = other_copy
    fields = ([chrom, str(start_0based), str(end_1based)] + ["0"] * 6
              + [other_chrom, str(other_start_0based), str(other_end_1based)] + ["0"] * 11
              + [str(frac_match)] + ["0"] * 20)
    return "\t".join(fields)


VOLLGER_2022_SD_BED_HEADER = "\t".join(
    ["#chr1", "start1", "end1", "name", "fakeScore", "strand1", "start1", "end1", "color", "chr2", "start2", "end2",
     "score", "strand2", "max_len", "aln_len", "indel_a", "indel_b", "alnB", "matchB", "mismatchB", "transitionsB",
     "transversions", "fracMatch", "fracMatchIndel", "jck", "k2K", "aln_gaps", "uppercaseA", "uppercaseB",
     "uppercaseMatches", "aln_matches", "aln_mismatches", "aln_gap_bases", "filter_score", "sat_bases", "unique_id",
     "original", "telo", "peri", "acro", "telo2", "peri2", "acro2"])


class TestAddSegdupMaxSequenceIdentityAnnotations(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def write_file(self, filename, lines):
        path = os.path.join(self.temp_dir, filename)
        with (gzip.open if path.endswith("gz") else open)(path, "wt") as f:
            f.write("\n".join(lines) + "\n")
        return path

    def test_get_masked_chry_pars(self):
        self.assertEqual(script.get_masked_chry_pars("hg38"), [("chrY", 10000, 2781479), ("chrY", 56887902, 57217415)])
        self.assertEqual(script.get_masked_chry_pars("chm13"), [("chrY", 0, 2458320), ("chrY", 62122809, 62460029)])

    def test_get_max_sequence_identity_of_overlapping_segdups(self):
        path = self.write_file("genomicSuperDups.txt.gz", [
            build_ucsc_genomic_super_dups_line("chr1", 100, 200, 0.95),
            build_ucsc_genomic_super_dups_line("chr1", 150, 300, 0.992918),
            build_ucsc_genomic_super_dups_line("chrUn_KI270742v1", 0, 100, 0.99),
        ])
        segdups_by_chrom, num_rows, num_rows_whose_other_copy_is_masked = script.read_segdups(
            path, "ucsc", script.get_masked_chry_pars("hg38"))
        self.assertEqual((num_rows, num_rows_whose_other_copy_is_masked), (3, 0))

        get_max_identity = script.get_max_sequence_identity_of_overlapping_segdups
        self.assertEqual(get_max_identity(segdups_by_chrom, "chr1", 110, 120), 0.95)
        self.assertEqual(get_max_identity(segdups_by_chrom, "chr1", 140, 160), 0.992918)  # overlaps both
        self.assertEqual(get_max_identity(segdups_by_chrom, "chr1", 90, 101), 0.95)  # 1 bp of overlap is enough
        self.assertEqual(get_max_identity(segdups_by_chrom, "chr1", 299, 400), 0.992918)
        self.assertIsNone(get_max_identity(segdups_by_chrom, "chr1", 300, 310))  # segdup intervals are half-open
        self.assertIsNone(get_max_identity(segdups_by_chrom, "chr1", 90, 100))
        self.assertIsNone(get_max_identity(segdups_by_chrom, "chr2", 100, 200))
        self.assertEqual(get_max_identity(segdups_by_chrom, "chrUn_KI270742v1", 10, 20), 0.99)

    def test_get_identity_range_label(self):
        get_label = script.get_identity_range_label
        self.assertEqual(get_label(None), "no segdup")
        self.assertEqual(get_label(0.9), "< 0.95")
        self.assertEqual(get_label(0.95), "0.95 to 0.98")
        self.assertEqual(get_label(0.98), "0.98 to 0.99")
        self.assertEqual(get_label(0.989999), "0.98 to 0.99")
        self.assertEqual(get_label(0.99), ">= 0.99")
        self.assertEqual(get_label(1.0), ">= 0.99")

    def test_read_segdups_from_vollger_2022_sd_bed(self):
        path = self.write_file("chm13v2.0_SD.full.bed", [
            VOLLGER_2022_SD_BED_HEADER,
            build_vollger_2022_sd_bed_line("chr1", 1, 13731, 0.957011),
            build_vollger_2022_sd_bed_line("chr9", 150595369, 150617116, 0.957011),
        ])
        segdups_by_chrom, _, _ = script.read_segdups(path, "vollger2022", script.get_masked_chry_pars("chm13"))
        self.assertEqual(sorted(segdups_by_chrom), ["chr1", "chr9"])
        self.assertEqual(script.get_max_sequence_identity_of_overlapping_segdups(
            segdups_by_chrom, "chr9", 150600000, 150600030), 0.957011)

    def test_read_segdups_skips_chm13_segdups_whose_other_copy_is_in_a_chry_par(self):
        # rows from chm13v2.0_SD.full.bed: chrX PAR1 paired with chrY PAR1, whose sequence chm13v2.0.maskedY masks
        path = self.write_file("chm13v2.0_SD.full.bed", [
            VOLLGER_2022_SD_BED_HEADER,
            build_vollger_2022_sd_bed_line("chrX", 1, 1307356, 0.990803, other_copy=("chrY", 622, 1328808)),
            build_vollger_2022_sd_bed_line("chrY", 622, 1328808, 0.990803, other_copy=("chrX", 1, 1307356)),
            build_vollger_2022_sd_bed_line("chrX", 153925814, 154259566, 0.998537,
                                           other_copy=("chrY", 62122789, 62457362)),
            build_vollger_2022_sd_bed_line("chrX", 1000, 3000, 0.95, other_copy=("chr2", 5000, 7000)),
            # the other copy ends 1 bp before chrY PAR2 starts, so it is kept
            build_vollger_2022_sd_bed_line("chr3", 100, 1300, 0.96, other_copy=("chrY", 62121609, 62122809)),
        ])
        segdups_by_chrom, num_rows, num_rows_whose_other_copy_is_masked = script.read_segdups(
            path, "vollger2022", script.get_masked_chry_pars("chm13"))
        self.assertEqual((num_rows, num_rows_whose_other_copy_is_masked), (5, 2))

        get_max_identity = script.get_max_sequence_identity_of_overlapping_segdups
        self.assertEqual(get_max_identity(segdups_by_chrom, "chrX", 1500, 1530), 0.95)
        self.assertIsNone(get_max_identity(segdups_by_chrom, "chrX", 154000000, 154000030))
        self.assertEqual(get_max_identity(segdups_by_chrom, "chrY", 1000, 1030), 0.990803)
        self.assertEqual(get_max_identity(segdups_by_chrom, "chr3", 200, 230), 0.96)

    def test_read_segdups_skips_hg38_segdups_whose_other_copy_is_in_a_chry_par(self):
        # a chrX PAR1 segdup paired with hg38's chrY PAR1 (chrY:10000-2781479), and one paired outside the chrY PARs
        path = self.write_file("hg38.chr_only.SDs.bed.gz", [
            VOLLGER_2022_SD_BED_HEADER,
            build_vollger_2022_sd_bed_line("chrX", 20000, 60000, 0.999, other_copy=("chrY", 20000, 60000)),
            build_vollger_2022_sd_bed_line("chrX", 20000, 30000, 0.96, other_copy=("chrY", 12500000, 12510000)),
        ])
        segdups_by_chrom, num_rows, num_rows_whose_other_copy_is_masked = script.read_segdups(
            path, "vollger2022", script.get_masked_chry_pars("hg38"))
        self.assertEqual((num_rows, num_rows_whose_other_copy_is_masked), (2, 1))
        self.assertEqual(script.get_max_sequence_identity_of_overlapping_segdups(
            segdups_by_chrom, "chrX", 25000, 25030), 0.96)

    def test_read_segdups_rejects_a_file_of_the_other_source(self):
        ucsc_path = self.write_file("genomicSuperDups.txt", [build_ucsc_genomic_super_dups_line("chr1", 100, 200, 0.95)])
        with self.assertRaises(ValueError):
            script.read_segdups(ucsc_path, "vollger2022", [])

        vollger_path = self.write_file("chm13v2.0_SD.full.bed", [build_vollger_2022_sd_bed_line("chr1", 100, 200, 0.95)])
        with self.assertRaises(ValueError):
            script.read_segdups(vollger_path, "ucsc", [])

    def test_read_segdups_rejects_a_header_without_frac_match_where_expected(self):
        header_fields = VOLLGER_2022_SD_BED_HEADER.split("\t")
        header_fields[23] = "fracMatchIndel"
        path = self.write_file("chm13v2.0_SD.full.bed", ["\t".join(header_fields)])
        with self.assertRaises(ValueError):
            script.read_segdups(path, "vollger2022", [])

    def test_read_segdups_rejects_an_identity_that_is_not_a_fraction(self):
        path = self.write_file("genomicSuperDups.txt", [build_ucsc_genomic_super_dups_line("chr1", 100, 200, 95)])
        with self.assertRaises(ValueError):
            script.read_segdups(path, "ucsc", [])

    def test_main(self):
        catalog_path = os.path.join(self.temp_dir, "catalog.json.gz")
        catalog_records = [
            {"LocusId": "1-140-160-AT", "ReferenceRegion": "chr1:140-160", "LocusStructure": "(AT)*",
             "VariantType": "Repeat", "Source": "SomeSource"},
            {"LocusId": "1-500-530-CAG", "ReferenceRegion": "chr1:500-530", "LocusStructure": "(CAG)*",
             "VariantType": "Repeat"},
        ]
        with gzip.open(catalog_path, "wt") as f:
            json.dump(catalog_records, f)

        # run the script once per source, as the hg38 pipeline does, so the output has both fields
        ucsc_path = self.write_file("genomicSuperDups.txt.gz", [
            build_ucsc_genomic_super_dups_line("chr1", 100, 200, 0.95),
            build_ucsc_genomic_super_dups_line("chr1", 150, 300, 0.992918),
        ])
        vollger_path = self.write_file("hg38.chr_only.SDs.bed.gz", [
            VOLLGER_2022_SD_BED_HEADER,
            build_vollger_2022_sd_bed_line("chr1", 120, 280, 0.9931),
        ])
        ucsc_output_path = os.path.join(self.temp_dir, "with_ucsc.json.gz")
        output_path = os.path.join(self.temp_dir, "output.json.gz")
        for segdups_path, segdups_source, input_path, output in [
                (ucsc_path, "ucsc", catalog_path, ucsc_output_path),
                (vollger_path, "vollger2022", ucsc_output_path, output_path)]:
            argv = ["script", "--segdups-file", segdups_path, "--segdups-source", segdups_source,
                    "--genome-version", "hg38", "-o", output, input_path]
            with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print"):
                script.main()

        with gzip.open(output_path, "rt") as f:
            output_records = json.load(f)

        self.assertEqual(output_records, [
            {**catalog_records[0], "UCSC_SegDupMaxSequenceIdentity": 0.992918,
             "Vollger2022_SegDupMaxSequenceIdentity": 0.9931},
            {**catalog_records[1], "UCSC_SegDupMaxSequenceIdentity": None,
             "Vollger2022_SegDupMaxSequenceIdentity": None},
        ])


if __name__ == "__main__":
    unittest.main()
