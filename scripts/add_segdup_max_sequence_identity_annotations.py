"""Add a UCSC_SegDupMaxSequenceIdentity or Vollger2022_SegDupMaxSequenceIdentity field to each locus: the highest
sequence identity among the segmental duplications (segdups) from that source that overlap it.

The segdups come from one of two sources, whose files both list each duplicated pair twice (once from each copy's
side), so a locus in either copy gets the pair's identity:
- ucsc: UCSC's hg38 genomicSuperDups table (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/
  genomicSuperDups.txt.gz), made with the WGAC method (Bailey et al. 2001, 2002): duplicated regions that are mostly
  1 kb or longer, with at least 90% identity. Besides the main chromosomes, it has pairs whose other copy is on an
  unplaced or unlocalized contig (chrUn_*, *_random), which are kept since those contigs are part of the reference
  that reads align to. It has no alt or fix patch contigs. UCSC has no table like it for CHM13: its CHM13 segdup
  track is the Vollger et al. 2022 set.
- vollger2022: the segdups of Vollger et al. 2022, made with SEDEF, also mostly 1 kb or longer with at least 90%
  identity. For CHM13 they're the T2T consortium's chm13v2.0_SD.full.bed (https://s3-us-west-2.amazonaws.com/
  human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.full.bed). For hg38 they're the paper's
  hg38.chr_only.SDs.bed.gz, from annotation_bed_files.zip in its Zenodo record (https://zenodo.org/records/5502036),
  which covers only the main chromosomes. The two files have the same columns.

The references that the hg38 and CHM13 catalogs are made from mask the chrY pseudoautosomal regions (PARs) with Ns, so
reads can't align to a copy there. Pairs whose other copy overlaps a chrY PAR by at least 1 bp are skipped, which leaves
the chrX PAR loci with the identities of their other segdups. The Vollger et al. 2022 files pair the chrX PARs with the
chrY PARs, while UCSC's hg38 table has no pairs whose other copy overlaps a chrY PAR. A pair whose other copy only
partly overlaps a chrY PAR is skipped whole, so the few loci next to a PAR boundary that only such a pair covers (2 in
the hg38 catalog, 1 in CHM13) get None.

The identity is each file's fracMatch column: the fraction of aligned bases that match, not counting indels. A segdup
counts if it overlaps the locus by at least 1 bp. Loci that overlap no segdup get None.
"""

import argparse
import collections
import gzip
import intervaltree

from add_pseudoautosomal_region_annotations import PSEUDOAUTOSOMAL_REGIONS_BY_GENOME_VERSION
from catalog_annotation_utils import write_catalog_with_a_field_added_to_each_record
from str_analysis.utils.misc_utils import parse_interval

# For each source of segdups: the field it fills, the 0-based column index of each value read from its file, and the
# number of columns the file's rows have
SEGDUP_SOURCES = {
    "ucsc": {"field_name": "UCSC_SegDupMaxSequenceIdentity", "chrom": 1, "start_0based": 2, "end_1based": 3,
             "other_chrom": 7, "other_start_0based": 8, "other_end_1based": 9, "fracMatch": 26, "num_columns": 30},
    "vollger2022": {"field_name": "Vollger2022_SegDupMaxSequenceIdentity", "chrom": 0, "start_0based": 1,
                    "end_1based": 2, "other_chrom": 9, "other_start_0based": 10, "other_end_1based": 11,
                    "fracMatch": 23, "num_columns": 44},
}


def get_masked_chry_pars(genome_version):
    """Returns the chrY PARs of a reference genome, which the catalog's reference masks with Ns.

    Args:
        genome_version (str): hg38 or chm13

    Returns:
        list: (chromosome, 0-based start, end) tuples, with chromosome names that start with chr as in the segdup files
    """
    return [(f"chr{chrom}", start_0based, end_1based) for chrom, start_0based, end_1based, _
            in PSEUDOAUTOSOMAL_REGIONS_BY_GENOME_VERSION[genome_version] if chrom == "Y"]


def read_segdups(segdups_path, segdups_source, masked_regions):
    """Reads the segdups of a UCSC genomicSuperDups table or a Vollger et al. 2022 SD BED file, skipping the ones whose
    other copy overlaps a masked region.

    Args:
        segdups_path (str): path of the segdup file, which can be gzipped
        segdups_source (str): a key of SEGDUP_SOURCES
        masked_regions (list): (chromosome, 0-based start, end) of the regions that the reference masks with Ns

    Returns:
        tuple: (dict that maps each chromosome to an IntervalTree of its segdups, with each segdup's fracMatch as the
            interval data, number of rows read, number of rows skipped because their other copy is masked)
    """
    columns = SEGDUP_SOURCES[segdups_source]
    segdups_by_chrom = collections.defaultdict(intervaltree.IntervalTree)
    num_rows = 0
    num_rows_whose_other_copy_is_masked = 0
    with (gzip.open if segdups_path.endswith("gz") else open)(segdups_path, "rt") as f:
        for line_number, line in enumerate(f, start=1):
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                # the Vollger et al. 2022 files' header line names their columns, so check that they're where expected
                if fields[columns["fracMatch"]] != "fracMatch":
                    raise ValueError(f"{segdups_path} header has {fields[columns['fracMatch']]} rather than "
                                     f"fracMatch in column {columns['fracMatch'] + 1}")
                continue

            if len(fields) != columns["num_columns"]:
                raise ValueError(f"{segdups_path} line {line_number} has {len(fields)} columns rather than the "
                                 f"{columns['num_columns']} of the {segdups_source} segdup files")

            frac_match = float(fields[columns["fracMatch"]])
            if not 0 < frac_match <= 1:
                raise ValueError(f"{segdups_path} line {line_number} has fracMatch {frac_match}, which is not a "
                                 f"fraction")

            num_rows += 1
            other_chrom = fields[columns["other_chrom"]]
            other_start_0based = int(fields[columns["other_start_0based"]])
            other_end_1based = int(fields[columns["other_end_1based"]])
            if any(other_chrom == masked_chrom and other_start_0based < masked_end and other_end_1based > masked_start
                   for masked_chrom, masked_start, masked_end in masked_regions):
                num_rows_whose_other_copy_is_masked += 1
                continue

            segdups_by_chrom[fields[columns["chrom"]]].addi(
                int(fields[columns["start_0based"]]), int(fields[columns["end_1based"]]), frac_match)

    return segdups_by_chrom, num_rows, num_rows_whose_other_copy_is_masked


def get_max_sequence_identity_of_overlapping_segdups(segdups_by_chrom, chrom, start_0based, end_1based):
    """Returns the highest fracMatch among the segdups that overlap the given interval.

    Args:
        segdups_by_chrom (dict): segdups returned by read_segdups
        chrom (str): chromosome of the interval
        start_0based (int): start of the interval
        end_1based (int): end of the interval

    Returns:
        float: the highest fracMatch, or None if no segdup overlaps the interval
    """
    if chrom not in segdups_by_chrom:
        return None

    return max((segdup.data for segdup in segdups_by_chrom[chrom].overlap(start_0based, end_1based)), default=None)


def get_identity_range_label(identity):
    """Returns the label of the range of segdup identities that the main() summary counts an identity in.

    Args:
        identity (float): a value of the field this script adds, or None for a locus that overlaps no segdup

    Returns:
        str: "no segdup", "< 0.95", "0.95 to 0.98", "0.98 to 0.99" or ">= 0.99"
    """
    if identity is None:
        return "no segdup"
    if identity >= 0.99:
        return ">= 0.99"
    if identity >= 0.98:
        return "0.98 to 0.99"
    if identity >= 0.95:
        return "0.95 to 0.98"
    return "< 0.95"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--segdups-file", required=True,
                        help="Path of the segdup file, eg. genomicSuperDups.txt.gz, hg38.chr_only.SDs.bed.gz or "
                             "chm13v2.0_SD.full.bed")
    parser.add_argument("--segdups-source", required=True, choices=sorted(SEGDUP_SOURCES),
                        help="ucsc for UCSC's genomicSuperDups table, or vollger2022 for a Vollger et al. 2022 SD BED "
                             "file. This determines the field that the script adds.")
    parser.add_argument("--genome-version", required=True, choices=sorted(PSEUDOAUTOSOMAL_REGIONS_BY_GENOME_VERSION),
                        help="Reference genome of the catalog, which determines the coordinates of the masked chrY "
                             "PARs")
    parser.add_argument("-o", "--output-catalog-json-path", required=True, help="Path of the output catalog JSON file")
    parser.add_argument("catalog_json_path", help="Path of the input catalog JSON file")
    args = parser.parse_args()

    segdups_by_chrom, num_rows, num_rows_whose_other_copy_is_masked = read_segdups(
        args.segdups_file, args.segdups_source, get_masked_chry_pars(args.genome_version))
    print(f"Read {num_rows:,d} segdups from {args.segdups_file}, and skipped {num_rows_whose_other_copy_is_masked:,d} "
          f"of them whose other copy overlaps a chrY PAR, which the reference masks with Ns")

    field_name = SEGDUP_SOURCES[args.segdups_source]["field_name"]
    counts_by_identity = write_catalog_with_a_field_added_to_each_record(
        args.catalog_json_path, args.output_catalog_json_path, field_name,
        lambda record: get_max_sequence_identity_of_overlapping_segdups(
            segdups_by_chrom, *parse_interval(record["ReferenceRegion"])))

    counts_by_identity_range = collections.Counter()
    for identity, count in counts_by_identity.items():
        counts_by_identity_range[get_identity_range_label(identity)] += count

    print(f"Wrote {sum(counts_by_identity_range.values()):,d} loci to {args.output_catalog_json_path}. "
          f"{field_name} counts:")
    for identity_range in "no segdup", "< 0.95", "0.95 to 0.98", "0.98 to 0.99", ">= 0.99":
        print(f"{counts_by_identity_range[identity_range]:12,d}  {identity_range}")


if __name__ == "__main__":
    main()
