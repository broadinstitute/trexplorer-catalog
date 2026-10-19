"""Remove the catalog loci that are within a given distance of either end of their chromosome.

Genotyping tools need flanking sequence on both sides of a locus: ExpansionHunter fails on a locus less than 1000bp from
a chromosome end (with an error like 'Unable to extract chrM:-793-207'), and trgt validate rejects a locus whose 50bp
flanks extend past a chromosome end. hg38's chromosome ends are runs of Ns, so its catalog has no such loci, but CHM13's
chromosome ends are telomeric repeats, which the CHM13 perfect repeats source includes.
"""

import argparse
import gzip
import ijson
import simplejson as json

from str_analysis.utils.misc_utils import parse_interval


def is_near_chromosome_end(record, chromosome_sizes, min_distance_to_chromosome_end):
    """Returns True if any of the record's ReferenceRegions is within min_distance_to_chromosome_end of either end of
    its chromosome.

    Args:
        record (dict): catalog record in ExpansionHunter format
        chromosome_sizes (dict): chromosome name => size in bases
        min_distance_to_chromosome_end (int): minimum number of bases between a locus and either chromosome end

    Returns:
        bool
    """
    reference_regions = record["ReferenceRegion"]
    for reference_region in reference_regions if isinstance(reference_regions, list) else [reference_regions]:
        chrom, start_0based, end_1based = parse_interval(reference_region)
        if start_0based < min_distance_to_chromosome_end or (
                end_1based > chromosome_sizes[chrom] - min_distance_to_chromosome_end):
            return True
    return False


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__.split("\n")[0])
    parser.add_argument("-R", "--reference-fasta", required=True, help="Reference genome FASTA, which must have a "
                        ".fai index. Only its .fai index is read")
    parser.add_argument("--min-distance-to-chromosome-end", type=int, default=1000, help="Remove loci that are closer "
                        "than this many bases to either end of their chromosome")
    parser.add_argument("-o", "--output-catalog-json-path", required=True, help="Path of the output catalog JSON file")
    parser.add_argument("catalog_json_path", help="Path of the input catalog JSON file")
    args = parser.parse_args()

    chromosome_sizes = {}
    with open(f"{args.reference_fasta}.fai") as f:
        for line in f:
            fields = line.split("\t")
            chromosome_sizes[fields[0]] = int(fields[1])

    input_count = 0
    output_count = 0
    fopen = gzip.open if args.catalog_json_path.endswith("gz") else open
    f2open = gzip.open if args.output_catalog_json_path.endswith("gz") else open
    with fopen(args.catalog_json_path, "rb") as f, f2open(args.output_catalog_json_path, "wt") as f2:
        f2.write("[")
        for record in ijson.items(f, "item"):
            input_count += 1
            if is_near_chromosome_end(record, chromosome_sizes, args.min_distance_to_chromosome_end):
                print(f"Removing {record['LocusId']} @ {record['ReferenceRegion']}")
                continue
            if output_count > 0:
                f2.write(", ")
            f2.write(json.dumps(record, use_decimal=True, indent=4))
            output_count += 1
        f2.write("]")

    print(f"Removed {input_count - output_count:,d} out of {input_count:,d} loci that are less than "
          f"{args.min_distance_to_chromosome_end:,d}bp from a chromosome end. Wrote {output_count:,d} loci to "
          f"{args.output_catalog_json_path}")


if __name__ == "__main__":
    main()
