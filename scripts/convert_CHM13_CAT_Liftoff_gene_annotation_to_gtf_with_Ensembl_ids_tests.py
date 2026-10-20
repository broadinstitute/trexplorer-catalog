import gzip
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

from str_analysis.utils.gtf_utils import compute_genomic_region_of_interval

import convert_CHM13_CAT_Liftoff_gene_annotation_to_gtf_with_Ensembl_ids as script


def make_gff3_line(source, feature_type, start_1based, end_1based, strand, **attributes):
    attributes_column = ";".join(f"{key}={value}" for key, value in attributes.items())
    return "\t".join(["chr1", source, feature_type, str(start_1based), str(end_1based), ".", strand, ".",
                      attributes_column])


# attributes of a CAT transcript of a GENCODE gene, a CAT novel isoform of the same gene, and a Liftoff rRNA gene copy
GENCODE_TRANSCRIPT_ATTRIBUTES = dict(
    source_transcript="ENST00000000001.2", source_gene="ENSG00000000001.5", source_gene_common_name="GENE1",
    biotype="protein_coding", transcript_biotype="protein_coding", transcript_id="CHM13_T0000001",
    gene_id="CHM13_G0000001", gene_name="GENE1")
NOVEL_ISOFORM_ATTRIBUTES = dict(
    source_transcript="N/A", source_gene="ENSG00000000001.5", source_gene_common_name="GENE1",
    biotype="protein_coding", transcript_biotype="protein_coding", transcript_id="CHM13_T0000002",
    gene_id="CHM13_G0000001", gene_name="MSTRG.1")
RRNA_GENE_COPY_ATTRIBUTES = dict(
    source_transcript="rRNA_T1", source_gene="LOFF_G0010000", source_gene_common_name="RNA45S1", biotype="rRNA",
    transcript_biotype="rRNA", transcript_id="LOFF_T0010000", gene_id="LOFF_G0010000", gene_name="RNA45S1")

GFF3_LINES = [
    "##gff-version 3",
    make_gff3_line("CAT", "gene", 100, 400, "+", ID="CHM13_G0000001", gene_id="CHM13_G0000001", gene_name="MSTRG.1"),
    make_gff3_line("CAT", "transcript", 100, 400, "+", ID="CHM13_T0000001", **GENCODE_TRANSCRIPT_ATTRIBUTES),
    make_gff3_line("CAT", "exon", 100, 200, "+", ID="exon:CHM13_T0000001:0", **GENCODE_TRANSCRIPT_ATTRIBUTES),
    make_gff3_line("CAT", "intron", 201, 299, "+", ID="intron:CHM13_T0000001:0", **GENCODE_TRANSCRIPT_ATTRIBUTES),
    make_gff3_line("CAT", "exon", 300, 400, "+", ID="exon:CHM13_T0000001:1", **GENCODE_TRANSCRIPT_ATTRIBUTES),
    make_gff3_line("CAT", "CDS", 150, 200, "+", ID="CDS:CHM13_T0000001:0", **GENCODE_TRANSCRIPT_ATTRIBUTES),
    # as in the real GFF3, the CDS records include the stop codon
    make_gff3_line("CAT", "CDS", 300, 353, "+", ID="CDS:CHM13_T0000001:1", **GENCODE_TRANSCRIPT_ATTRIBUTES),
    make_gff3_line("CAT", "start_codon", 150, 152, "+", ID="start:CHM13_T0000001", **GENCODE_TRANSCRIPT_ATTRIBUTES),
    make_gff3_line("CAT", "stop_codon", 351, 353, "+", ID="stop:CHM13_T0000001", **GENCODE_TRANSCRIPT_ATTRIBUTES),
    make_gff3_line("CAT", "transcript", 100, 900, "+", ID="CHM13_T0000002", **NOVEL_ISOFORM_ATTRIBUTES),
    make_gff3_line("CAT", "exon", 100, 900, "+", ID="exon:CHM13_T0000002:0", **NOVEL_ISOFORM_ATTRIBUTES),
    make_gff3_line("Liftoff", "transcript", 1000, 1100, "-", ID="LOFF_T0010000", **RRNA_GENE_COPY_ATTRIBUTES),
    make_gff3_line("Liftoff", "exon", 1000, 1100, "-", ID="exon:LOFF_T0010000:0", **RRNA_GENE_COPY_ATTRIBUTES),
]

GENCODE_TRANSCRIPT_GTF_ATTRIBUTES = ('gene_id "ENSG00000000001.5"; transcript_id "ENST00000000001.2"; '
                                     'gene_name "GENE1"; gene_type "protein_coding"; transcript_type "protein_coding";')
RRNA_GENE_COPY_GTF_ATTRIBUTES = ('gene_id "LOFF_G0010000"; transcript_id "LOFF_T0010000"; gene_name "RNA45S1"; '
                                 'gene_type "rRNA"; transcript_type "rRNA";')


class TestConvertCHM13CATLiftoffGeneAnnotationToGtf(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_parse_gff3_attributes(self):
        self.assertEqual(script.parse_gff3_attributes("ID=CHM13_T1;tag=basic%2CCCDS;gene_name=A-B;"), {
            "ID": "CHM13_T1", "tag": "basic,CCDS", "gene_name": "A-B"})

    def test_get_gtf_ids_of_transcript(self):
        self.assertEqual(script.get_gtf_ids_of_transcript({"ID": "CHM13_T0000001", **GENCODE_TRANSCRIPT_ATTRIBUTES}), {
            "gene_id": "ENSG00000000001.5", "transcript_id": "ENST00000000001.2", "gene_name": "GENE1",
            "gene_type": "protein_coding", "transcript_type": "protein_coding"})
        self.assertEqual(script.get_gtf_ids_of_transcript({"ID": "LOFF_T0010000", **RRNA_GENE_COPY_ATTRIBUTES}), {
            "gene_id": "LOFF_G0010000", "transcript_id": "LOFF_T0010000", "gene_name": "RNA45S1",
            "gene_type": "rRNA", "transcript_type": "rRNA"})
        self.assertIsNone(script.get_gtf_ids_of_transcript({"ID": "CHM13_T0000002", **NOVEL_ISOFORM_ATTRIBUTES}))
        self.assertIsNone(script.get_gtf_ids_of_transcript({
            **NOVEL_ISOFORM_ATTRIBUTES, "ID": "CHM13_T1000167", "source_gene": "None",
            "transcript_biotype": "unknown_likely_coding"}))
        with self.assertRaises(ValueError):
            script.get_gtf_ids_of_transcript({
                **GENCODE_TRANSCRIPT_ATTRIBUTES, "ID": "CHM13_T0000001", "source_gene_common_name": "GENE 1"})

    def test_remove_stop_codons_from_CDS_intervals(self):
        # a stop codon at the end of a + strand transcript's CDS, and at the start of a - strand transcript's CDS
        self.assertEqual(script.remove_stop_codons_from_CDS_intervals([(300, 353), (150, 200)], [(351, 353)]),
                         [(150, 200), (300, 350)])
        self.assertEqual(script.remove_stop_codons_from_CDS_intervals([(100, 200), (300, 400)], [(100, 102)]),
                         [(103, 200), (300, 400)])
        # a stop codon split by an intron, whose two records remove the end of one CDS record and all of the next one
        self.assertEqual(script.remove_stop_codons_from_CDS_intervals([(100, 200), (300, 301)], [(200, 200), (300, 301)]),
                         [(100, 199)])
        # a transcript without a stop codon record keeps its CDS
        self.assertEqual(script.remove_stop_codons_from_CDS_intervals([(100, 200)], []), [(100, 200)])

    def test_compute_UTR_intervals(self):
        # the exons, CDS and UTRs of GENCODE v50's ENST00000641515 (OR4F5), whose 3' UTR starts with the stop codon
        self.assertEqual(
            script.compute_UTR_intervals([(65419, 65433), (65520, 65573), (69037, 71585)],
                                         [(65565, 65573), (69037, 70005)]),
            [(65419, 65433), (65520, 65564), (70006, 71585)])
        # a CDS that spans whole exons leaves no UTR, and a non-coding transcript has no UTRs
        self.assertEqual(script.compute_UTR_intervals([(100, 200), (300, 400)], [(100, 200), (300, 400)]), [])
        self.assertEqual(script.compute_UTR_intervals([(100, 200)], []), [])

    def test_main(self):
        gff3_path = os.path.join(self.temp_dir, "annotation.gff3.gz")
        gtf_path = os.path.join(self.temp_dir, "annotation.gtf.gz")
        with gzip.open(gff3_path, "wt") as f:
            f.write("\n".join(GFF3_LINES) + "\n")

        with mock.patch.object(sys, "argv", ["script", "-o", gtf_path, gff3_path]), mock.patch("builtins.print"):
            script.main()

        with gzip.open(gtf_path, "rt") as f:
            gtf_lines = f.read().splitlines()

        # the novel isoform, and the gene, intron, start_codon and stop_codon records are left out
        self.assertEqual(gtf_lines, [
            f"chr1\tCAT\t{feature_type}\t{start_1based}\t{end_1based}\t.\t+\t.\t{GENCODE_TRANSCRIPT_GTF_ATTRIBUTES}"
            for feature_type, start_1based, end_1based in [
                ("transcript", 100, 400), ("exon", 100, 200), ("exon", 300, 400), ("CDS", 150, 200), ("CDS", 300, 350),
                ("UTR", 100, 149), ("UTR", 351, 400),
            ]
        ] + [
            f"chr1\tLiftoff\t{feature_type}\t1000\t1100\t.\t-\t.\t{RRNA_GENE_COPY_GTF_ATTRIBUTES}"
            for feature_type in ("transcript", "exon")
        ])

        # str_analysis reads the GTF, labels the UTRs as 5' and 3' UTRs, and drops the Ensembl version suffixes
        with mock.patch("builtins.print"):
            self.assertEqual(compute_genomic_region_of_interval("chr1", 120, 130, gtf_path),
                             ("5' UTR", "GENE1", "ENSG00000000001", "ENST00000000001"))
            self.assertEqual(compute_genomic_region_of_interval("chr1", 160, 170, gtf_path),
                             ("CDS", "GENE1", "ENSG00000000001", "ENST00000000001"))
            self.assertEqual(compute_genomic_region_of_interval("chr1", 360, 380, gtf_path),
                             ("3' UTR", "GENE1", "ENSG00000000001", "ENST00000000001"))
            self.assertEqual(compute_genomic_region_of_interval("chr1", 1020, 1030, gtf_path),
                             ("exon", "RNA45S1", "LOFF_G0010000", "LOFF_T0010000"))


if __name__ == "__main__":
    unittest.main()
