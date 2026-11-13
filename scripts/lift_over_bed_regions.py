"""Lift the regions of a BED file over to another reference genome with a UCSC chain file, keeping each region's name.

This is for annotation regions that can be megabases long, like the centromere and satellite regions of the Miga lab's
hg38 CenSat BED file, rather than for repeat loci (scripts/lift_over_catalog.py lifts those by their flanking bases, and
drops a locus whose two sides map through different chains). Each region is lifted to the target bases that the chain
file aligns to its bases on the same chromosome, so the source and target references must name their chromosomes the
same way (as hg38 and hg19 do). Pieces of a region that a chain aligns to another chromosome are dropped: they come from
short alignments of satellite sequence to paralogous arrays elsewhere (eg. 8,112 pieces of the hg38 CenSat regions,
7.1 Mb in total, would otherwise land on other chromosomes in hg19, and label loci there with unrelated classes).
Within a chain, the target bases between two aligned pieces of the region (a chain gap, eg. an insertion in the target
reference) are included too when there are at most --max-chain-gap-to-fill of them, so a region that one chain aligns
with only small gaps lifts to one region. Larger chain gaps split the lifted region, as do different chains. Without
that limit, an 11 kb hg38 cenSat region on chr9 would lift to a 2.1 Mb hg19 region, since the chain that aligns it
skips over megabases of hg19 sequence in the segmental duplications next to chr9's centromere. Parts of a region that
no chain aligns are dropped, as is a region that no chain aligns at all. Most of an hg38 centromere model, which hg19
lacks since its centromeres are gaps, is dropped this way, but the parts of it that a chain aligns to hg19's
pericentromeric satellite sequence lift there.

Chains also align satellite sequence to paralogous copies far away on the same chromosome (eg. GRCh38's falsely
duplicated copies of 21q22 on chr21's short arm are 'ct' regions, and lift onto the real 21q22 of hg19, where they'd
label every locus in genes like CBS and U2AF1 'ct'). So each lifted piece is also lifted back with --reverse-chain-file,
and kept only if at least half of its bases map back into the input regions. A piece of the real 21q22 maps back to the
real 21q22 of hg38, which has no CenSat region, and is dropped. A piece of hg19's pericentromeric satellite sequence
maps back into an hg38 satellite region, though not always the one it was lifted from (satellite arrays align to many
copies of themselves), and is kept.

The output BED has 4 columns: chromosome, start, end and the input region's name. It's sorted by position.
"""

import argparse
import bisect
import collections
import gzip
import intervaltree

from lift_over_catalog import map_position, parse_chain_file

# The smallest fraction of a lifted piece's bases that must map back into the input regions for the piece to be kept
MIN_FRACTION_OF_BASES_THAT_MAP_BACK_INTO_THE_REGIONS = 0.5


def lift_region(chain_index, chrom, start_0based, end_1based, max_chain_gap_to_fill):
    """Lifts a region to the target bases that the chain file aligns to it on the same chromosome, filling small chain
    gaps.

    Args:
        chain_index (dict): output of lift_over_catalog.parse_chain_file
        chrom (str): source chromosome of the region
        start_0based (int): 0-based start of the region in the source reference
        end_1based (int): 1-based end of the region in the source reference
        max_chain_gap_to_fill (int): the most target bases between two aligned pieces of the region in the same chain
            that the lifted region includes rather than splitting at them

    Returns:
        list: (target chrom, target start_0based, target end_1based) tuples, sorted by target position
    """
    if chrom not in chain_index:
        return []

    block_starts, blocks = chain_index[chrom]
    target_pieces_by_chain_id = collections.defaultdict(list)
    # the blocks of a source chromosome don't overlap, so the first block that can overlap the region is the last one
    # that starts at or before the region's start
    for block in blocks[max(0, bisect.bisect_right(block_starts, start_0based) - 1):]:
        if block.source_start_0based >= end_1based:
            break
        if block.target_chrom != chrom:
            continue
        overlap_start_0based = max(start_0based, block.source_start_0based)
        overlap_end_1based = min(end_1based, block.source_end)
        if overlap_end_1based <= overlap_start_0based:
            continue

        first_base_target_position = map_position(chain_index, chrom, overlap_start_0based).position_0based
        last_base_target_position = map_position(chain_index, chrom, overlap_end_1based - 1).position_0based
        target_pieces_by_chain_id[block.chain_id].append((
            block.target_chrom,
            min(first_base_target_position, last_base_target_position),
            max(first_base_target_position, last_base_target_position) + 1))

    lifted_regions = []
    for target_pieces in target_pieces_by_chain_id.values():
        # a chain's pieces are all on one target chromosome
        target_pieces.sort()
        current_chrom, current_start_0based, current_end_1based = target_pieces[0]
        for _, piece_start_0based, piece_end_1based in target_pieces[1:]:
            if piece_start_0based - current_end_1based <= max_chain_gap_to_fill:
                current_end_1based = max(current_end_1based, piece_end_1based)
            else:
                lifted_regions.append((current_chrom, current_start_0based, current_end_1based))
                current_start_0based, current_end_1based = piece_start_0based, piece_end_1based
        lifted_regions.append((current_chrom, current_start_0based, current_end_1based))

    return sorted(lifted_regions)


def count_bases_that_map_back_into_the_regions(reverse_chain_index, source_regions_by_chrom, chrom, start_0based,
                                               end_1based):
    """Lifts a lifted piece back to the source reference, and counts how many of its bases land in the input regions.

    The piece is lifted back without filling chain gaps, so that only the source bases that the piece's own bases align
    to are counted.

    Args:
        reverse_chain_index (dict): output of lift_over_catalog.parse_chain_file for the chain file from the target back
            to the source reference
        source_regions_by_chrom (dict): maps each source chromosome to an IntervalTree of the input regions, with
            overlapping regions merged
        chrom (str): chromosome of the lifted piece, which lift_region keeps the same as its region's
        start_0based (int): 0-based start of the lifted piece in the target reference
        end_1based (int): 1-based end of the lifted piece in the target reference

    Returns:
        int: the number of the piece's bases that map back into the input regions
    """
    if chrom not in source_regions_by_chrom:
        return 0
    return sum(min(back_end_1based, region.end) - max(back_start_0based, region.begin)
               for _, back_start_0based, back_end_1based in lift_region(
                   reverse_chain_index, chrom, start_0based, end_1based, max_chain_gap_to_fill=0)
               for region in source_regions_by_chrom[chrom].overlap(back_start_0based, back_end_1based))


def lift_bed_regions(chain_index, reverse_chain_index, bed_lines, max_chain_gap_to_fill):
    """Lifts the regions of a BED file, keeping the lifted pieces whose bases mostly map back into the input regions.

    Args:
        chain_index (dict): output of lift_over_catalog.parse_chain_file
        reverse_chain_index (dict): output of lift_over_catalog.parse_chain_file for the chain file from the target back
            to the source reference
        bed_lines (iterable): lines of the input BED file, whose track and comment lines are skipped
        max_chain_gap_to_fill (int): see lift_region

    Returns:
        2-tuple: (sorted list of (chrom, start_0based, end_1based, name) tuples, collections.Counter of how many regions
            lifted to 0, 1 or more regions, and how many lifted pieces were dropped since they don't map back)
    """
    regions = [line.rstrip("\n").split("\t") for line in bed_lines
               if line.strip() and not line.startswith(("track", "#"))]
    source_regions_by_chrom = collections.defaultdict(intervaltree.IntervalTree)
    for fields in regions:
        source_regions_by_chrom[fields[0]].addi(int(fields[1]), int(fields[2]))
    for source_regions in source_regions_by_chrom.values():
        source_regions.merge_overlaps()

    lifted_regions = []
    counters = collections.Counter()
    for fields in regions:
        target_spans = []
        for target_chrom, start_0based, end_1based in lift_region(
                chain_index, fields[0], int(fields[1]), int(fields[2]), max_chain_gap_to_fill):
            if count_bases_that_map_back_into_the_regions(
                    reverse_chain_index, source_regions_by_chrom, target_chrom, start_0based, end_1based
            ) >= MIN_FRACTION_OF_BASES_THAT_MAP_BACK_INTO_THE_REGIONS * (end_1based - start_0based):
                target_spans.append((target_chrom, start_0based, end_1based))
            else:
                counters["lifted pieces dropped since less than half of their bases map back into the regions"] += 1
        counters["regions"] += 1
        if not target_spans:
            counters["regions with no lifted piece left, so were dropped"] += 1
        elif len(target_spans) > 1:
            counters["regions split by large chain gaps or several chains into more than one region"] += 1
        lifted_regions.extend((target_chrom, start_0based, end_1based, fields[3])
                              for target_chrom, start_0based, end_1based in target_spans)

    return sorted(lifted_regions), counters


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--chain-file", required=True, help="UCSC chain file from the source to the target reference, "
                        "eg. hg38ToHg19.over.chain.gz")
    parser.add_argument("--reverse-chain-file", required=True, help="UCSC chain file from the target back to the "
                        "source reference, eg. hg19ToHg38.over.chain.gz, for checking that each lifted piece maps "
                        "back into the input regions")
    parser.add_argument("--max-chain-gap-to-fill", type=int, default=10_000, help="The most target bases between two "
                        "aligned pieces of a region in the same chain that the lifted region includes rather than "
                        "being split at them")
    parser.add_argument("-o", "--output-bed-path", required=True, help="Path of the output BED file, which is gzipped "
                        "if it ends in .gz")
    parser.add_argument("bed_path", help="Path of the input BED file, which can be gzipped. Its 4th column is the name "
                        "that each lifted region keeps.")
    args = parser.parse_args()

    chain_index = parse_chain_file(args.chain_file)
    reverse_chain_index = parse_chain_file(args.reverse_chain_file)
    with (gzip.open if args.bed_path.endswith("gz") else open)(args.bed_path, "rt") as f:
        lifted_regions, counters = lift_bed_regions(chain_index, reverse_chain_index, f, args.max_chain_gap_to_fill)

    with (gzip.open if args.output_bed_path.endswith("gz") else open)(args.output_bed_path, "wt") as f:
        for chrom, start_0based, end_1based, name in lifted_regions:
            f.write(f"{chrom}\t{start_0based}\t{end_1based}\t{name}\n")

    print(f"Wrote {len(lifted_regions):,d} regions to {args.output_bed_path}")
    for key, count in sorted(counters.items(), key=lambda x: -x[1]):
        print(f"{count:12,d}  {key}")


if __name__ == "__main__":
    main()
