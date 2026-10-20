import gzip
import os
import shutil
import simplejson as json
import sys
import tempfile
import unittest
from unittest import mock

import add_hg38_locus_ids_using_lifted_over_hg38_catalog as script


def make_record(reference_region, motif, hg38_locus_id=None, **other_fields):
    chrom, start_and_end = reference_region.split(":")
    start_0based, end_1based = start_and_end.split("-")
    record = {
        "LocusId": f"{chrom.replace('chr', '')}-{start_0based}-{end_1based}-{motif}",
        "ReferenceRegion": reference_region,
        "LocusStructure": f"({motif})*",
        "VariantType": "Repeat",
        **other_fields,
    }
    if hg38_locus_id is not None:
        record["Hg38LocusId"] = hg38_locus_id
    return record


# known disease-associated loci, whose TRExplorerV2 field holds the hg38 TRExplorer catalog's LocusId for each locus
KNOWN_DISEASE_ASSOCIATED_LOCI = [
    {"LocusId": "XYLT1", "Diseases": [{"Symbol": "BSS"}], "TRExplorerV2": "16-17470907-17470922-GCC"},
    {"LocusId": "EP400__1-500-530-AAG", "Diseases": [{"Symbol": "X"}], "TRExplorerV2": "1-610-630-AAG"},
    {"LocusId": "TMEM185A", "Diseases": [{"Symbol": "FRAXF"}]},
]


# hg38 loci lifted over to CHM13: CHM13 coordinates, with the hg38 LocusId in the Hg38LocusId field
LIFTED_HG38_RECORDS = [
    make_record("chr1:100-130", "CAG", "1-200-230-CAG"),     # same boundaries as a catalog locus
    make_record("chr1:97-133", "CAG", "1-197-233-CAG"),      # wider version of the same locus
    make_record("chr1:300-330", "AT", "1-400-430-AT"),       # different motif from the catalog locus
    make_record("chr1:500-520", "AAG", "1-600-620-AAG"),     # narrower, Jaccard similarity 20/30
    make_record("chr1:700-710", "AAG", "1-800-810-AAG"),     # narrower, Jaccard similarity 10/30
    make_record("chr1:903-933", "CAG", "1-1003-1033-CAG"),   # same size, shifted by 3bp
    make_record("chr2:96-140", "CAG", "2-196-240-CAG"),      # tied Jaccard similarity, smaller start
    make_record("chr2:100-144", "CAG", "2-200-244-CAG"),     # tied Jaccard similarity, larger start
    make_record("chr3:1000-1100", "AAAAAGGGGG", "3-2000-2100-AAAAAGGGGG"),  # VNTR with another motif of the same size
    make_record("chr3:3000-3100", "AAAAAGGGGG", "3-4000-4100-AAAAAGGGGG"),  # VNTR with another motif, same boundaries
    make_record("chr3:3005-3100", "AAAAACCCCC", "3-4005-4100-AAAAACCCCC"),  # VNTR with the same motif, narrower
    make_record("chr4:100-130", "CAG", "4-150-180-CAG"),     # the lifted locus that the merge kept
]


class TestAddHg38LocusIdsUsingLiftedOverHg38Catalog(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.lifted_over_catalog_path = os.path.join(self.temp_dir, "lifted_over.json.gz")
        with gzip.open(self.lifted_over_catalog_path, "wt") as f:
            json.dump(LIFTED_HG38_RECORDS, f)
        self.known_disease_associated_loci_path = os.path.join(self.temp_dir, "known_disease_associated_loci.json")
        with open(self.known_disease_associated_loci_path, "wt") as f:
            json.dump(KNOWN_DISEASE_ASSOCIATED_LOCI, f)

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def find_best_match(self, reference_region, motif, motif_length_match_sufficient_for_VNTRs=False):
        lifted_hg38_loci_by_chrom = script.read_lifted_over_hg38_loci(self.lifted_over_catalog_path)
        lifted_interval = script.find_best_matching_lifted_hg38_locus(
            lifted_hg38_loci_by_chrom,
            make_record(reference_region, motif),
            min_overlap_fraction=0.66,
            min_jaccard_similarity=0.66,
            motif_length_match_sufficient_for_VNTRs=motif_length_match_sufficient_for_VNTRs)
        return lifted_interval.data["Hg38LocusId"] if lifted_interval is not None else None

    def test_read_lifted_over_hg38_loci(self):
        lifted_hg38_loci_by_chrom = script.read_lifted_over_hg38_loci(self.lifted_over_catalog_path)
        self.assertEqual(sorted(lifted_hg38_loci_by_chrom), ["chr1", "chr2", "chr3", "chr4"])
        self.assertEqual(sum(len(t) for t in lifted_hg38_loci_by_chrom.values()), len(LIFTED_HG38_RECORDS))
        self.assertEqual(sorted(lifted_hg38_loci_by_chrom["chr4"])[0].data, {
            "ReferenceRegion": "chr4:100-130", "LocusStructure": "(CAG)*", "Hg38LocusId": "4-150-180-CAG"})

    def test_read_hg38_locus_ids_of_known_disease_associated_loci(self):
        self.assertEqual(
            script.read_hg38_locus_ids_of_known_disease_associated_loci(self.known_disease_associated_loci_path),
            {"XYLT1": "16-17470907-17470922-GCC", "EP400__1-500-530-AAG": "1-610-630-AAG"})

    def test_compute_jaccard_similarity(self):
        self.assertEqual(script.compute_jaccard_similarity(100, 130, 100, 130), 1)
        self.assertAlmostEqual(script.compute_jaccard_similarity(500, 530, 500, 520), 20/30)
        self.assertAlmostEqual(script.compute_jaccard_similarity(900, 930, 903, 933), 27/33)
        self.assertEqual(script.compute_jaccard_similarity(100, 130, 130, 160), 0)

    def test_get_hg38_locus_match_type(self):
        self.assertEqual(script.get_hg38_locus_match_type(100, 130, 100, 130), "SameBoundaries")
        self.assertEqual(script.get_hg38_locus_match_type(100, 130, 97, 133), "Hg38LocusIsWider")
        self.assertEqual(script.get_hg38_locus_match_type(500, 530, 500, 520), "Hg38LocusIsNarrower")
        self.assertEqual(script.get_hg38_locus_match_type(900, 930, 903, 933), "Hg38LocusIsShifted")

    def test_find_best_matching_lifted_hg38_locus(self):
        # the lifted locus with the same boundaries beats the wider one
        self.assertEqual(self.find_best_match("chr1:100-130", "CAG"), "1-200-230-CAG")
        # for a locus that neither lifted locus has the same boundaries as, the one with the higher Jaccard similarity wins
        self.assertEqual(self.find_best_match("chr1:98-132", "CAG"), "1-197-233-CAG")
        # a lifted locus with a different motif doesn't match
        self.assertIsNone(self.find_best_match("chr1:300-330", "CAG"))
        # Jaccard similarity 20/30 passes the 0.66 threshold, but 10/30 doesn't
        self.assertEqual(self.find_best_match("chr1:500-530", "AAG"), "1-600-620-AAG")
        self.assertIsNone(self.find_best_match("chr1:700-730", "AAG"))
        self.assertEqual(self.find_best_match("chr1:900-930", "CAG"), "1-1003-1033-CAG")
        # a tie in Jaccard similarity goes to the lifted locus with the smaller start
        self.assertEqual(self.find_best_match("chr2:100-140", "CAG"), "2-196-240-CAG")
        # VNTRs with different motifs of the same size match only with motif_length_match_sufficient_for_VNTRs
        self.assertIsNone(self.find_best_match("chr3:1000-1100", "AAAAACCCCC"))
        self.assertEqual(self.find_best_match("chr3:1000-1100", "AAAAACCCCC", motif_length_match_sufficient_for_VNTRs=True),
                         "3-2000-2100-AAAAAGGGGG")
        # a lifted locus with the same motif beats one with a different motif of the same size, even when the latter
        # has the same boundaries
        self.assertEqual(self.find_best_match("chr3:3000-3100", "AAAAACCCCC", motif_length_match_sufficient_for_VNTRs=True),
                         "3-4005-4100-AAAAACCCCC")
        # no lifted loci on the chromosome, or none overlapping
        self.assertIsNone(self.find_best_match("chr5:100-130", "CAG"))
        self.assertIsNone(self.find_best_match("chr1:5000-5030", "CAG"))

    def test_main(self):
        catalog_path = os.path.join(self.temp_dir, "catalog.json.gz")
        output_path = os.path.join(self.temp_dir, "output.json.gz")
        catalog_records = [
            make_record("chr1:100-130", "CAG"),
            make_record("chr1:300-330", "CAG"),
            # a lifted locus that matches takes precedence over the locus's KnownDiseaseAssociatedLocus label
            make_record("chr1:500-530", "AAG", KnownDiseaseAssociatedLocus="EP400__1-500-530-AAG"),
            # no lifted locus matches, so the Hg38LocusId comes from the KnownDiseaseAssociatedLocus label
            make_record("chr1:700-730", "AAG", KnownDiseaseAssociatedLocus="XYLT1"),
            # no lifted locus matches, and the label has no TRExplorerV2 value
            make_record("chr6:100-130", "CAG", KnownDiseaseAssociatedLocus="TMEM185A"),
            make_record("chr1:900-930", "CAG"),
            make_record("chr4:100-130", "CAG", "4-150-180-CAG"),
            make_record("chr5:10-40", "CAG", "5-20-50-CAG"),
        ]
        with gzip.open(catalog_path, "wt") as f:
            json.dump(catalog_records, f)

        argv = ["script",
                "--lifted-over-hg38-catalog", self.lifted_over_catalog_path,
                "--min-overlap-fraction", "0.66",
                "--min-jaccard-similarity", "0.66",
                "--motif-length-match-sufficient-for-VNTRs",
                "--known-disease-associated-loci", self.known_disease_associated_loci_path,
                "-o", output_path,
                catalog_path]
        with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print"):
            script.main()

        with gzip.open(output_path, "rt") as f:
            output_records = json.load(f)

        self.assertEqual(output_records, [
            {**catalog_records[0], "Hg38LocusId": "1-200-230-CAG", "Hg38LocusMatchType": "SameBoundaries"},
            {**catalog_records[1], "Hg38LocusId": None, "Hg38LocusMatchType": None},
            {**catalog_records[2], "Hg38LocusId": "1-600-620-AAG", "Hg38LocusMatchType": "Hg38LocusIsNarrower"},
            {**catalog_records[3], "Hg38LocusId": "16-17470907-17470922-GCC",
             "Hg38LocusMatchType": "SameKnownDiseaseAssociatedLocus"},
            {**catalog_records[4], "Hg38LocusId": None, "Hg38LocusMatchType": None},
            {**catalog_records[5], "Hg38LocusId": "1-1003-1033-CAG", "Hg38LocusMatchType": "Hg38LocusIsShifted"},
            # loci that already have an Hg38LocusId keep it, even when no overlapping lifted locus has that id
            {**catalog_records[6], "Hg38LocusMatchType": "SameBoundaries"},
            {**catalog_records[7], "Hg38LocusMatchType": None},
        ])


if __name__ == "__main__":
    unittest.main()
