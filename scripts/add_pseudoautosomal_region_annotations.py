"""Add a PseudoautosomalRegion field to each locus: PAR1 or PAR2 for a locus inside a pseudoautosomal region (PAR).

The two PARs at the ends of chrX and chrY have the same sequence on both chromosomes, so males have two copies of the
TRs in them, as for autosomal TRs. A locus gets the name of the PAR that contains its whole reference region, or None
if it lies outside both PARs, including a locus that crosses a PAR boundary. The PARs of chrY are masked with Ns in the
references that the hg38 and CHM13 catalogs are made from, so in practice only chrX loci get a value.
"""

import argparse

from catalog_annotation_utils import write_catalog_with_a_field_added_to_each_record
from str_analysis.utils.misc_utils import parse_interval

FIELD_NAME = "PseudoautosomalRegion"

# The (chromosome, 0-based start, end, PAR name) of each PAR, by reference genome
PSEUDOAUTOSOMAL_REGIONS_BY_GENOME_VERSION = {
    # From the Genome Reference Consortium's GRCh38 assembly regions file, with its 1-based starts made 0-based:
    # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_regions.txt
    "hg38": [
        ("X", 10000, 2781479, "PAR1"),
        ("X", 155701382, 156030895, "PAR2"),
        ("Y", 10000, 2781479, "PAR1"),
        ("Y", 56887902, 57217415, "PAR2"),
    ],
    # From the T2T consortium's BED file of the CHM13v2.0 PARs:
    # https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_PAR.bed
    "chm13": [
        ("X", 0, 2394410, "PAR1"),
        ("X", 153925834, 154259566, "PAR2"),
        ("Y", 0, 2458320, "PAR1"),
        ("Y", 62122809, 62460029, "PAR2"),
    ],
}


def get_pseudoautosomal_region(pseudoautosomal_regions, chrom, start_0based, end_1based):
    """Returns the name of the PAR that contains the given interval.

    Args:
        pseudoautosomal_regions (list): the (chromosome, 0-based start, end, PAR name) of each PAR, with chromosome
            names that don't start with "chr", as in PSEUDOAUTOSOMAL_REGIONS_BY_GENOME_VERSION
        chrom (str): chromosome of the interval, with or without the "chr" prefix
        start_0based (int): start of the interval
        end_1based (int): end of the interval

    Returns:
        str: the name of the PAR (eg. PAR1) that contains the whole interval, or None if no PAR does
    """
    chrom = chrom.removeprefix("chr")
    for par_chrom, par_start_0based, par_end_1based, par_name in pseudoautosomal_regions:
        if chrom == par_chrom and par_start_0based <= start_0based and end_1based <= par_end_1based:
            return par_name

    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--genome-version", required=True, choices=sorted(PSEUDOAUTOSOMAL_REGIONS_BY_GENOME_VERSION),
                        help="Reference genome of the catalog, which determines the PAR coordinates")
    parser.add_argument("-o", "--output-catalog-json-path", required=True, help="Path of the output catalog JSON file")
    parser.add_argument("catalog_json_path", help="Path of the input catalog JSON file")
    args = parser.parse_args()

    pseudoautosomal_regions = PSEUDOAUTOSOMAL_REGIONS_BY_GENOME_VERSION[args.genome_version]
    counts_by_region = write_catalog_with_a_field_added_to_each_record(
        args.catalog_json_path, args.output_catalog_json_path, FIELD_NAME,
        lambda record: get_pseudoautosomal_region(
            pseudoautosomal_regions, *parse_interval(record["ReferenceRegion"])))

    print(f"Wrote {sum(counts_by_region.values()):,d} loci to {args.output_catalog_json_path}. {FIELD_NAME} counts:")
    for region_name, count in sorted(counts_by_region.items(), key=lambda x: -x[1]):
        print(f"{count:12,d}  {region_name}")


if __name__ == "__main__":
    main()
