import gzip
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

import extract_MANE_transcripts_from_GENCODE_gtf as script

GENE_LINE = ('chr1\tHAVANA\tgene\t65419\t71585\t.\t+\t.\tgene_id "ENSG00000186092.7_8"; gene_type "protein_coding"; '
             'gene_name "OR4F5"; level 2;\n')
MANE_SELECT_TRANSCRIPT_LINE = (
    'chr1\tHAVANA\ttranscript\t65419\t71585\t.\t+\t.\tgene_id "ENSG00000186092.7_8"; '
    'transcript_id "ENST00000641515.2_3"; gene_name "OR4F5"; tag "basic"; tag "MANE_Select"; tag "Ensembl_canonical";\n')
MANE_SELECT_EXON_LINE_WITH_TAG_AT_END = (
    'chr1\tHAVANA\texon\t65419\t65433\t.\t+\t.\tgene_id "ENSG00000186092.7_8"; '
    'transcript_id "ENST00000641515.2_3"; gene_name "OR4F5"; tag "MANE_Select"\n')
MANE_PLUS_CLINICAL_CDS_LINE = (
    'chr1\tHAVANA\tCDS\t65565\t65573\t.\t+\t0\tgene_id "ENSG00000186092.7_8"; '
    'transcript_id "ENST00000641516.2_1"; gene_name "OR4F5"; tag "MANE_Plus_Clinical";\n')
OTHER_TRANSCRIPT_LINE = (
    'chr1\tHAVANA\ttranscript\t65419\t71585\t.\t+\t.\tgene_id "ENSG00000186092.7_8"; '
    'transcript_id "ENST00000999999.1_1"; gene_name "OR4F5"; tag "basic"; tag "MANE_Select_candidate";\n')


class TestExtractMANETranscriptsFromGENCODEGtf(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_is_MANE_transcript_record(self):
        self.assertTrue(script.is_MANE_transcript_record(MANE_SELECT_TRANSCRIPT_LINE))
        self.assertTrue(script.is_MANE_transcript_record(MANE_SELECT_EXON_LINE_WITH_TAG_AT_END))
        self.assertTrue(script.is_MANE_transcript_record(MANE_PLUS_CLINICAL_CDS_LINE))
        self.assertFalse(script.is_MANE_transcript_record(GENE_LINE))
        # a tag that starts with "MANE_Select" isn't the MANE_Select tag
        self.assertFalse(script.is_MANE_transcript_record(OTHER_TRANSCRIPT_LINE))

    def test_main(self):
        gtf_path = os.path.join(self.temp_dir, "gencode.gtf.gz")
        output_path = os.path.join(self.temp_dir, "mane.gtf.gz")
        with gzip.open(gtf_path, "wt") as f:
            f.write("##description: test\n")
            for line in GENE_LINE, MANE_SELECT_TRANSCRIPT_LINE, MANE_SELECT_EXON_LINE_WITH_TAG_AT_END, \
                    OTHER_TRANSCRIPT_LINE, MANE_PLUS_CLINICAL_CDS_LINE:
                f.write(line)

        with mock.patch.object(sys, "argv", ["script", "-o", output_path, gtf_path]), mock.patch("builtins.print"):
            script.main()

        with gzip.open(output_path, "rt") as f:
            self.assertEqual(f.read(), "##description: test\n" + MANE_SELECT_TRANSCRIPT_LINE +
                             MANE_SELECT_EXON_LINE_WITH_TAG_AT_END + MANE_PLUS_CLINICAL_CDS_LINE)


if __name__ == "__main__":
    unittest.main()
