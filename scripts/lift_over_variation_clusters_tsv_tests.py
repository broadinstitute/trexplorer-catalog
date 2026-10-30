import collections
import gzip
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

import lift_over_variation_clusters_tsv as script
from lift_over_catalog import ChainBlock


def make_chain_index(blocks_by_source_chrom):
    """Builds the output of lift_over_catalog.parse_chain_file from lists of ChainBlocks."""
    return {chrom: ([block.source_start_0based for block in blocks], blocks)
            for chrom, blocks in blocks_by_source_chrom.items()}


# chr1:0-10,000 maps to chr1:1,000-11,000 of the target, and chr2:0-10,000 maps to chr5:0-10,000
CHAIN_INDEX = make_chain_index({
    "chr1": [ChainBlock(0, 10_000, "chr1", 1_000, "+", 1_000_000, "1")],
    "chr2": [ChainBlock(0, 10_000, "chr5", 0, "+", 1_000_000, "2")],
})


def target_locus(locus_id, chrom, start_0based, end_1based, motif):
    return script.TargetLocus(locus_id, chrom, start_0based, end_1based, [motif])


class TestLiftOverVariationClustersTsv(unittest.TestCase):

    def test_read_target_loci_by_source_locus_id_prefers_same_boundaries(self):
        records = [
            {"LocusId": "1-1100-1130-CAG", "ReferenceRegion": "chr1:1100-1130", "LocusStructure": "(CAG)*",
             "Hg38LocusId": "1-100-130-CAG", "Hg38LocusMatchType": "Hg38LocusIsNarrower"},
            {"LocusId": "1-1100-1124-CAG", "ReferenceRegion": "chr1:1100-1124", "LocusStructure": "(CAG)*",
             "Hg38LocusId": "1-100-130-CAG", "Hg38LocusMatchType": "SameBoundaries"},
            {"LocusId": "1-1200-1220-A", "ReferenceRegion": "chr1:1200-1220", "LocusStructure": "(A)*"},
            {"LocusId": "HTT", "ReferenceRegion": ["chr1:1300-1330", "chr1:1336-1350"],
             "LocusStructure": "(CAG)*CAACAG(CCG)*", "Hg38LocusId": "HTT", "Hg38LocusMatchType": "SameBoundaries"},
        ]
        self.assertEqual(script.read_target_loci_by_source_locus_id(records), {
            "1-100-130-CAG": target_locus("1-1100-1124-CAG", "chr1", 1100, 1124, "CAG"),
            "HTT": script.TargetLocus("HTT", "chr1", 1300, 1350, ["CAG", "CCG"]),
        })

    def test_parse_and_format_table_rows(self):
        row = script.parse_table_row(
            "ID=1-19175-19184-TCC;MOTIFS=TCC;STRUC=(TCC)n\tchr1:19175-19184\t4.0\t0\tchr1:19171-19184\n")
        self.assertEqual(row, script.TableRow("1-19175-19184-TCC", "chr1:19175-19184", "4.0", "0", "chr1:19171-19184"))
        self.assertEqual(script.parse_table_row("ID=1-10000-10108-TAACCC;MOTIFS=TAACCC;STRUC=(TAACCC)n\t"
                                                "chr1:10000-10108\t\tDEPTH\t\n").vc_end_offset, "DEPTH")
        self.assertEqual(
            script.format_row_without_variation_cluster(target_locus("1-1100-1124-CAG", "chr1", 1100, 1124, "CAG")),
            "ID=1-1100-1124-CAG;MOTIFS=CAG;STRUC=(CAG)n\tchr1:1100-1124\t0.0\t0\tchr1:1100-1124")

    def test_lift_variation_cluster(self):
        counters = collections.Counter()
        members = [target_locus("1-1100-1124-CAG", "chr1", 1100, 1124, "CAG"),
                   target_locus("1-1130-1150-A", "chr1", 1130, 1150, "A")]
        # the vc_region chr1:95-150 lifts to chr1:1095-1150
        self.assertEqual(script.lift_variation_cluster("chr1:95-150", members, CHAIN_INDEX, counters), [
            "ID=1-1100-1124-CAG;MOTIFS=CAG;STRUC=(CAG)n\tchr1:1100-1124\t5.0\t26\tchr1:1095-1150",
            "ID=1-1130-1150-A;MOTIFS=A;STRUC=(A)n\tchr1:1130-1150\t35.0\t0\tchr1:1095-1150",
        ])
        self.assertEqual(counters["variation clusters lifted over"], 1)

    def test_lift_variation_cluster_widens_it_to_contain_its_members(self):
        counters = collections.Counter()
        members = [target_locus("1-1090-1124-CAG", "chr1", 1090, 1124, "CAG")]
        self.assertEqual(script.lift_variation_cluster("chr1:95-150", members, CHAIN_INDEX, counters), [
            "ID=1-1090-1124-CAG;MOTIFS=CAG;STRUC=(CAG)n\tchr1:1090-1124\t0.0\t26\tchr1:1090-1150"])
        self.assertEqual(counters["variation clusters widened to contain the CHM13 loci of their members"], 1)

    def test_lift_variation_cluster_that_does_not_lift_over(self):
        counters = collections.Counter()
        members = [target_locus("1-1100-1124-CAG", "chr1", 1100, 1124, "CAG")]
        # chr3 isn't in the chain file, so the cluster's members are written as loci without a variation cluster
        self.assertEqual(script.lift_variation_cluster("chr3:95-150", members, CHAIN_INDEX, counters), [
            "ID=1-1100-1124-CAG;MOTIFS=CAG;STRUC=(CAG)n\tchr1:1100-1124\t0.0\t0\tchr1:1100-1124"])
        self.assertEqual(sum(v for k, v in counters.items() if k.startswith("variation clusters whose vc_region")), 1)

    def test_lift_variation_cluster_with_a_member_on_another_chromosome(self):
        counters = collections.Counter()
        members = [target_locus("5-100-124-CAG", "chr5", 100, 124, "CAG"),
                   target_locus("7-500-520-A", "chr7", 500, 520, "A")]
        self.assertEqual(script.lift_variation_cluster("chr2:95-150", members, CHAIN_INDEX, counters), [
            "ID=5-100-124-CAG;MOTIFS=CAG;STRUC=(CAG)n\tchr5:100-124\t5.0\t26\tchr5:95-150",
            "ID=7-500-520-A;MOTIFS=A;STRUC=(A)n\tchr7:500-520\t0.0\t0\tchr7:500-520",
        ])
        self.assertEqual(counters["members whose CHM13 locus doesn't overlap their lifted variation cluster"], 1)

    def test_lift_variation_cluster_with_a_far_away_member_on_the_same_chromosome(self):
        # Like hg38 chr9:63724055-63724597, whose only member lifted through another chain to a locus 30Mb away. The
        # member is written as a locus without a variation cluster, rather than stretching the cluster to reach it.
        counters = collections.Counter()
        members = [target_locus("1-1100-1124-CAG", "chr1", 1100, 1124, "CAG"),
                   target_locus("1-900000-900020-A", "chr1", 900_000, 900_020, "A")]
        self.assertEqual(script.lift_variation_cluster("chr1:95-150", members, CHAIN_INDEX, counters), [
            "ID=1-1100-1124-CAG;MOTIFS=CAG;STRUC=(CAG)n\tchr1:1100-1124\t5.0\t26\tchr1:1095-1150",
            "ID=1-900000-900020-A;MOTIFS=A;STRUC=(A)n\tchr1:900000-900020\t0.0\t0\tchr1:900000-900020",
        ])
        self.assertEqual(counters["members whose CHM13 locus doesn't overlap their lifted variation cluster"], 1)

    def test_lift_variation_cluster_whose_region_lifts_to_a_much_larger_one(self):
        # chr1:0-100 and chr1:100-20,000 are aligned to chr1:0-100 and chr1:20,100-40,000, so the target has a 20kb
        # insertion between the two flanks of the vc_region chr1:90-110
        chain_index = make_chain_index({"chr1": [
            ChainBlock(0, 100, "chr1", 0, "+", 1_000_000, "1"),
            ChainBlock(100, 20_000, "chr1", 20_100, "+", 1_000_000, "1"),
        ]})
        counters = collections.Counter()
        members = [target_locus("1-95-105-A", "chr1", 95, 105, "A")]
        self.assertEqual(script.lift_variation_cluster("chr1:90-110", members, chain_index, counters), [
            "ID=1-95-105-A;MOTIFS=A;STRUC=(A)n\tchr1:95-105\t0.0\t0\tchr1:95-105"])
        self.assertEqual(counters["variation clusters whose vc_region lifted to a region more than 10,000bp larger"], 1)

    def test_lift_variation_clusters_table(self):
        table_rows = [
            script.TableRow("1-90-94-AT", "chr1:90-94", "", "DEPTH", ""),
            script.TableRow("1-100-124-CAG", "chr1:100-124", "5.0", "26", "chr1:95-150"),
            script.TableRow("1-200-210-A", "chr1:200-210", "0.0", "0", "chr1:200-210"),
            script.TableRow("1-130-150-A", "chr1:130-150", "35.0", "0", "chr1:95-150"),
            script.TableRow("1-300-310-A", "chr1:300-310", "0.0", "0", "chr1:300-310"),
        ]
        target_loci = {
            "1-90-94-AT": target_locus("1-1090-1094-AT", "chr1", 1090, 1094, "AT"),
            "1-100-124-CAG": target_locus("1-1100-1124-CAG", "chr1", 1100, 1124, "CAG"),
            "1-130-150-A": target_locus("1-1130-1150-A", "chr1", 1130, 1150, "A"),
            "1-200-210-A": target_locus("1-1200-1210-A", "chr1", 1200, 1210, "A"),
        }
        counters = collections.Counter()
        self.assertEqual(script.lift_variation_clusters_table(table_rows, target_loci, CHAIN_INDEX, counters), [
            "ID=1-1090-1094-AT;MOTIFS=AT;STRUC=(AT)n\tchr1:1090-1094\t\tDEPTH\t",
            "ID=1-1100-1124-CAG;MOTIFS=CAG;STRUC=(CAG)n\tchr1:1100-1124\t5.0\t26\tchr1:1095-1150",
            "ID=1-1130-1150-A;MOTIFS=A;STRUC=(A)n\tchr1:1130-1150\t35.0\t0\tchr1:1095-1150",
            "ID=1-1200-1210-A;MOTIFS=A;STRUC=(A)n\tchr1:1200-1210\t0.0\t0\tchr1:1200-1210",
        ])
        self.assertEqual(counters["source table rows whose locus isn't in the lifted catalog"], 1)
        self.assertEqual(counters["variation clusters"], 1)
        self.assertEqual(counters["loci not in a variation cluster"], 1)

    def test_main(self):
        temp_dir = tempfile.mkdtemp()
        try:
            chain_path = os.path.join(temp_dir, "test.over.chain")
            with open(chain_path, "wt") as f:
                f.write("chain 1000 chr1 1000000 + 0 10000 chr1 1000000 + 1000 11000 1\n10000\n\n")
            catalog_path = os.path.join(temp_dir, "catalog.json.gz")
            with gzip.open(catalog_path, "wt") as f:
                f.write('[{"LocusId": "1-1100-1124-CAG", "ReferenceRegion": "chr1:1100-1124", '
                        '"LocusStructure": "(CAG)*", "Hg38LocusId": "1-100-124-CAG", '
                        '"Hg38LocusMatchType": "SameBoundaries"}]')
            table_path = os.path.join(temp_dir, "genome_clusters.tsv")
            with open(table_path, "wt") as f:
                f.write("region_info\toriginal_region\tvc_start_offset\tvc_end_offset\tvc_region\n")
                f.write("ID=1-100-124-CAG;MOTIFS=CAG;STRUC=(CAG)n\tchr1:100-124\t5.0\t26\tchr1:95-150\n")
            output_path = os.path.join(temp_dir, "lifted.tsv.gz")

            with mock.patch.object(sys, "argv", [
                    "lift_over_variation_clusters_tsv.py", "--chain-file", chain_path, "--lifted-catalog",
                    catalog_path, "-o", output_path, table_path]):
                script.main()

            with gzip.open(output_path, "rt") as f:
                self.assertEqual(f.read(), "region_info\toriginal_region\tvc_start_offset\tvc_end_offset\tvc_region\n"
                                           "ID=1-1100-1124-CAG;MOTIFS=CAG;STRUC=(CAG)n\tchr1:1100-1124\t5.0\t26\t"
                                           "chr1:1095-1150\n")
        finally:
            shutil.rmtree(temp_dir)


if __name__ == "__main__":
    unittest.main()
