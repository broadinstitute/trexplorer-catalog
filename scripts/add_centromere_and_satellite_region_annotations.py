"""Add a CentromereAndSatelliteRegion field to each locus, from a centromeric satellite (CenSat) annotation BED file.

The CenSat annotations from the Miga lab (https://github.com/hloucks/CenSatData), made with the same workflow for hg38
and CHM13, divide each genome's centromeres and satellite arrays into regions of these classes:
- alpha satellite: higher-order repeat arrays (active_hor for the active array, hor, and dhor for diverged ones),
  monomeric alpha satellite (mon), and alpha satellite that can't be put into one of those bins (mixedAlpha)
- human satellites 1A, 1B, 2 and 3 (HSat1A, HSat1B, HSat2, HSat3, and HSat2_3 for regions of both HSat2 and HSat3)
- beta and gamma satellites (bSat, gSat), and other centromeric satellites (cenSat)
- rDNA arrays (rDNA)
- centromeric transition regions (ct): the sequence between the satellite arrays of a centromere
- assembly gaps (GAP). A region annotated with more than one class has the classes joined by commas, eg. GAP,HSat3.

Each locus gets the class of the region it overlaps the most, without the subfamily details in parentheses (eg.
active_hor rather than active_hor(S1C1/5/19H1L)), or None if it overlaps no region.
"""

import argparse
import collections
import gzip
import ijson
import intervaltree
import simplejson as json

from str_analysis.utils.misc_utils import parse_interval

FIELD_NAME = "CentromereAndSatelliteRegion"


def parse_region_class(region_name):
    """Returns the class of a CenSat region, without the subfamily details in parentheses.

    Args:
        region_name (str): name column of the CenSat BED file, eg. active_hor(S1C1/5/19H1L)

    Returns:
        str: the region's class, eg. active_hor
    """
    return region_name.split("(", 1)[0]


def read_centromere_and_satellite_regions(bed_path):
    """Reads the regions of a CenSat BED file.

    Args:
        bed_path (str): path of the CenSat BED file, which can be gzipped and start with a track line

    Returns:
        dict: maps each chromosome to an IntervalTree of its regions, with each region's class as the interval data
    """
    regions_by_chrom = collections.defaultdict(intervaltree.IntervalTree)
    with (gzip.open if bed_path.endswith("gz") else open)(bed_path, "rt") as f:
        for line in f:
            if line.startswith(("track", "#")):
                continue
            fields = line.rstrip("\n").split("\t")
            regions_by_chrom[fields[0]].addi(int(fields[1]), int(fields[2]), parse_region_class(fields[3]))

    return regions_by_chrom


def get_class_of_region_with_largest_overlap(regions_by_chrom, chrom, start_0based, end_1based):
    """Returns the class of the region that overlaps the given interval the most.

    Args:
        regions_by_chrom (dict): regions returned by read_centromere_and_satellite_regions
        chrom (str): chromosome of the interval
        start_0based (int): start of the interval
        end_1based (int): end of the interval

    Returns:
        str: the class of the region with the largest overlap, or of the leftmost one among regions with the same
            overlap, or None if no region overlaps the interval
    """
    if chrom not in regions_by_chrom:
        return None

    overlapping_regions = regions_by_chrom[chrom].overlap(start_0based, end_1based)
    if not overlapping_regions:
        return None

    return max(overlapping_regions, key=lambda region: (
        min(region.end, end_1based) - max(region.begin, start_0based), -region.begin)).data


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--centromere-and-satellite-regions-bed", required=True,
                        help="Path of a CenSat annotation BED file, eg. hg38.cenSatv1.0.bed")
    parser.add_argument("-o", "--output-catalog-json-path", required=True, help="Path of the output catalog JSON file")
    parser.add_argument("catalog_json_path", help="Path of the input catalog JSON file")
    args = parser.parse_args()

    regions_by_chrom = read_centromere_and_satellite_regions(args.centromere_and_satellite_regions_bed)
    print(f"Read {sum(len(t) for t in regions_by_chrom.values()):,d} regions from "
          f"{args.centromere_and_satellite_regions_bed}")

    # stream the records one at a time, since loading a whole annotated catalog takes about 20GB of memory
    counts_by_class = collections.Counter()
    with (gzip.open if args.catalog_json_path.endswith("gz") else open)(args.catalog_json_path, "rb") as f, \
            (gzip.open if args.output_catalog_json_path.endswith("gz") else open)(
                args.output_catalog_json_path, "wt") as f2:
        f2.write("[")
        for i, record in enumerate(ijson.items(f, "item")):
            chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
            record[FIELD_NAME] = get_class_of_region_with_largest_overlap(
                regions_by_chrom, chrom, start_0based, end_1based)
            counts_by_class[record[FIELD_NAME]] += 1

            if i > 0:
                f2.write(", ")
            f2.write(json.dumps(record, use_decimal=True, indent=4))
        f2.write("]")

    print(f"Wrote {sum(counts_by_class.values()):,d} loci to {args.output_catalog_json_path}. {FIELD_NAME} counts:")
    for region_class, count in sorted(counts_by_class.items(), key=lambda x: -x[1]):
        print(f"{count:12,d}  {region_class}")


if __name__ == "__main__":
    main()
