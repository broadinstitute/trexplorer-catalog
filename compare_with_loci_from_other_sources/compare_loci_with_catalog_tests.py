import collections
import intervaltree
import unittest
from compare_loci_with_catalog import compare_loci

class TestCompareLociWithTrexplorerCatalog(unittest.TestCase):
    
    def test_compare_loci(self):
        chrom = "1"

        start = 100
        end = 200
        canonical_motif = "AGC"
        trexplorer_loci = collections.defaultdict(intervaltree.IntervalTree)
        trexplorer_loci[chrom].add(intervaltree.Interval(start, end, data=canonical_motif))

        new_catalog = [
            [chrom, end, end + 10, canonical_motif],
            [chrom, start, end, canonical_motif],
            [chrom, start - 3, end + 3, canonical_motif],
            [chrom, start + 1, end + 1, "ACC"],
            [chrom, start + 1, end + 30, canonical_motif],
            [chrom, start, end + 50, canonical_motif],
            [chrom, start + 1, end + 60, canonical_motif],
            [chrom, start, end + 100, canonical_motif],
        ]
        df = compare_loci(trexplorer_loci, new_catalog)
        self.assertEqual(list(df.overlap), [
            "absent",
            "exact match", 
            "within 2 repeats",
            "within 2 repeats",
            "union/overlap < 1.5x",
            "1.5x <= union/overlap < 2x",
            "1.5x <= union/overlap < 2x",
            "2x <= union/overlap < 3x",
            
        ])
        self.assertEqual(list(df.motif_match), [
            "absent",
            "same motif", 
            "same motif",
            "same motif length",
            "same motif",
            "same motif",
            "same motif",
            "same motif",
        ])
        #self.assertEqual(len(df), 10)
