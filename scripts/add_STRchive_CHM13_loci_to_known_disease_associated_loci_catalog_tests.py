import unittest

import add_STRchive_CHM13_loci_to_known_disease_associated_loci_catalog as script


def make_strchive_locus(strchive_id, gnomad, start_t2t=1000, stop_t2t=1030, motif="CAG"):
    return {
        "id": strchive_id, "gnomad": gnomad, "chrom": "chr3", "start_t2t": start_t2t, "stop_t2t": stop_t2t,
        "reference_motif_reference_orientation": [motif],
    }


class TestMakeRecordsForLociThatDidNotLiftOver(unittest.TestCase):

    def test_adds_STRchive_definitions_of_disease_loci_that_did_not_lift_over(self):
        hg38_records = [
            {"LocusId": "HTT", "Diseases": [{"Symbol": "HD"}]},                     # lifted over
            {"LocusId": "CNBP", "Diseases": [{"Symbol": "DM2"}], "Gene": "CNBP", "GeneRegion": "intron",
             "GeneId": "ENSG00000169714", "TRExplorerV2": "3-129172576-129172656-CAGG"},  # didn't lift over
            {"LocusId": "DMD", "Diseases": []},                                     # no diseases
            {"LocusId": "RUNX2__6-45422750-45422801-GCN", "Diseases": [{"Symbol": "CCD"}]},  # alternative definition
            {"LocusId": "TMEM185A", "Diseases": [{"Symbol": "FRAXF"}]},             # no STRchive CHM13 definition
        ]
        lifted_records = [{"LocusId": "HTT"}]
        strchive_loci = [
            make_strchive_locus("HD_HTT", ["HTT"]),
            make_strchive_locus("DM2_CNBP", ["CNBP"], motif="CAGG"),
            make_strchive_locus("DMD_DMD", ["DMD"]),
            make_strchive_locus("CCD_RUNX2", ["RUNX2__6-45422750-45422801-GCN"]),
            make_strchive_locus("FRAXF_TMEM185A", ["TMEM185A"], start_t2t=None),
        ]

        new_records, locus_ids_without_STRchive_definition = script.make_records_for_loci_that_did_not_lift_over(
            hg38_records, lifted_records, strchive_loci)

        self.assertEqual(new_records, [{
            "LocusId": "CNBP", "LocusStructure": "(CAGG)*", "ReferenceRegion": "chr3:1000-1030",
            "MainReferenceRegion": "chr3:1000-1030", "VariantType": "Repeat", "RepeatUnit": "CAGG", "Gene": "CNBP",
            "GeneRegion": "intron", "GeneId": "ENSG00000169714", "Diseases": [{"Symbol": "DM2"}],
            "TRExplorerV2": "3-129172576-129172656-CAGG", "STRchiveId": "DM2_CNBP",
        }])
        self.assertEqual(locus_ids_without_STRchive_definition, ["TMEM185A"])


if __name__ == "__main__":
    unittest.main()
