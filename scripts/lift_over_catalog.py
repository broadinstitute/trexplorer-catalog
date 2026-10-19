"""Lift a tandem repeat catalog over to another reference genome using a UCSC chain file.

This follows the approach of BCFtools/liftover (Genovese et al. 2024, https://doi.org/10.1093/bioinformatics/btae038).
Rather than mapping every base of a repeat, it maps the two bases that flank the repeat (the anchors), and the lifted
repeat is everything between the anchors' new positions. The two reference genomes come from different people, so a
repeat often has a different length in each, and a chain file usually puts its gap inside or at the edge of such a
repeat. Mapping only the flanking bases keeps those loci.

As in BCFtools/liftover, the anchors are the bases just outside the maximally extended repeat: when the bases next to
the repeat continue its period (like the CA of the CAACAG after HTT's CAG repeat), the aligner could have placed the
length difference after them, so mapping the base right next to the repeat could land inside the longer repeat of the
new reference and truncate it. The anchors therefore skip past those bases (up to MAX_FLANKING_BASES_TO_SKIP, and never
into another repeat of the same locus), and the same number of bases is then trimmed from the lifted repeat. Since the
chain's alignments allow mismatches, it can also place a length difference beyond those bases (for example, FGF14's GAA
repeat is 84bp longer in CHM13, and starts with AAAG interruptions in hg38). So, on each side where no whole copy of the
motif continues the old repeat, the lifted repeat is also extended over the adjacent whole copies of the motif in the new
reference. The chain can likewise put the boundary between two adjacent repeats of one locus in the wrong place (for
example, 3 of FXN's 9 GAA copies land at the end of its poly-A repeat), so that boundary is moved to where the two
repeats' sequences best match their motifs.

When a chain gap also covers an anchor, the anchor is moved outward, away from the repeat, to the nearest base within
--max-anchor-search-distance that is covered by the same chain as the other anchor. The lifted repeat then also
includes the new reference's version of the flanking bases that were skipped, so it can be slightly wider than the
repeat, but never cuts into it. (BCFtools/liftover instead places such an anchor by realigning the sequence locally,
which can move it into the repeat when the skipped flanking bases resemble the repeat.)

A locus is dropped when:
- an anchor isn't covered by the chain file even after moving it outward, or the two anchors map through different
  chains
- the repeat is absent from the new reference (the anchors map next to each other), or is shorter than one motif
- the new reference sequence contains Ns (for example, a hard-masked chrY PAR)
- the repeat's purity in the new reference is more than --max-purity-decrease below its purity in the old reference,
  which is what a mis-mapped locus or a locus that is not a repeat in the new reference looks like. Purity is computed
  for the best rotation of the motif, since the repeat can start at a different point of the motif in the new reference
- it has more than one repeat (eg. "(CAG)*CAACAG(CCG)*"), or sequence outside its repeat, and maps to the reverse
  strand, since rebuilding such a LocusStructure on the other strand isn't supported
- its new interval overlaps the new interval of a locus that was on another chromosome, or more than
  MIN_DISTANCE_BETWEEN_PARALOGOUS_LOCI away, in the old reference. This happens when the chain maps two paralogous
  copies onto the same repeat, and then there's no way to tell which of the two old loci the repeat corresponds to, so
  both are dropped. (Loci that were close together in the old reference can also land on overlapping intervals, for
  example two G runs separated by a base that is a G in the new reference, and those are kept, unless their lifted
  definitions are identical, in which case only the first is kept.)

Each repeat's motif is replaced by its rotation (or, on the reverse strand, the rotation of its reverse complement)
that best matches the repeat's sequence in the new reference, keeping the original motif when no other rotation matches
better.

LocusIds are left as they are, and also copied to the --original-locus-id-field. A coordinate-based LocusId (like
"1-10000-10108-TAACCC") therefore still names the old coordinates, so regenerate those with
str_analysis.annotate_and_filter_str_catalog --set-locus-id. The exception is a LocusId that ends in a coordinate-based
id after "__", like the known disease-associated loci catalog's alternative definitions (eg.
"EP400__12-132062524-132062611-CAG"): with --update-coordinates-in-locus-ids, that id is replaced by the lifted
coordinates and motif (eg. "EP400__12-132111361-132111451-CAG").

The chain file's target is the old reference and its query is the new reference, as in UCSC .over.chain files. For
hg38 to CHM13, use hg38ToHs1.over.chain.gz rather than hg38-chm13v2.over.chain.gz, which has erroneous mappings
caused by a chaintools bug (see the paper above).
"""

import argparse
import bisect
import collections
import gzip
import pysam
import re

from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator, parse_motifs_from_locus_structure
from str_analysis.utils.export_json import export_json
from str_analysis.utils.find_motif_utils import compute_best_phase_repeat_purity
from str_analysis.utils.misc_utils import parse_interval, reverse_complement

CORE_FIELDS = ("LocusId", "LocusStructure", "ReferenceRegion", "VariantType", "VariantId")

# The most bases next to a repeat that the anchors skip because they continue the repeat's period
MAX_FLANKING_BASES_TO_SKIP = 1000

# Two loci that land on overlapping intervals in the new reference are treated as paralogous copies only if they were
# on different chromosomes, or more than this many bases apart, in the old reference
MIN_DISTANCE_BETWEEN_PARALOGOUS_LOCI = 1000

ChainBlock = collections.namedtuple("ChainBlock", [
    "source_start_0based", "source_end", "target_chrom", "target_start_0based_on_strand", "target_strand",
    "target_chrom_size", "chain_id"])

MappedPosition = collections.namedtuple("MappedPosition", ["chrom", "position_0based", "strand", "chain_id"])

# A LocusId that ends in a coordinate-based id after "__", like "EP400__12-132062524-132062611-CAG"
LOCUS_ID_WITH_COORDINATES_REGEX = re.compile(r"^(.+__)[^_-]+-\d+-\d+-[A-Z]+$")


def parse_chain_file(chain_file_path):
    """Parses a UCSC chain file into the aligned blocks of each source chromosome, sorted by start.

    The chain format is described at https://genome.ucsc.edu/goldenPath/help/chain.html. In an .over.chain file,
    every source base maps to at most one target position, so the blocks of a source chromosome don't overlap.

    Args:
        chain_file_path (str): path of the chain file, optionally gzipped

    Returns:
        dict: source chromosome => (list of block start positions, list of ChainBlock tuples), both sorted by start
    """
    blocks_by_source_chrom = collections.defaultdict(list)
    fopen = gzip.open if chain_file_path.endswith("gz") else open
    with fopen(chain_file_path, "rt") as f:
        for line in f:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "chain":
                (_, _, source_chrom, _, source_strand, source_start, _, target_chrom, target_chrom_size,
                 target_strand, target_start, _, chain_id) = fields
                if source_strand != "+":
                    raise ValueError(f"Expected the source strand to be '+' in chain header: {line.strip()}")
                source_position = int(source_start)
                target_position = int(target_start)
                target_chrom_size = int(target_chrom_size)
                continue

            block_size = int(fields[0])
            blocks_by_source_chrom[source_chrom].append(ChainBlock(
                source_position, source_position + block_size, target_chrom, target_position, target_strand,
                target_chrom_size, chain_id))
            if len(fields) == 3:
                source_position += block_size + int(fields[1])
                target_position += block_size + int(fields[2])

    chain_index = {}
    for source_chrom, blocks in blocks_by_source_chrom.items():
        blocks.sort()
        chain_index[source_chrom] = ([block.source_start_0based for block in blocks], blocks)

    return chain_index


def map_position(chain_index, chrom, position_0based):
    """Maps one base of the source reference to the target reference.

    Args:
        chain_index (dict): output of parse_chain_file
        chrom (str): source chromosome
        position_0based (int): 0-based position of the base in the source reference

    Returns:
        MappedPosition: the base's position in the target reference, or None if no chain block covers it
    """
    if chrom not in chain_index:
        return None

    block_starts, blocks = chain_index[chrom]
    i = bisect.bisect_right(block_starts, position_0based) - 1
    if i < 0 or position_0based >= blocks[i].source_end:
        return None

    block = blocks[i]
    position_on_strand = block.target_start_0based_on_strand + (position_0based - block.source_start_0based)
    if block.target_strand == "-":
        target_position_0based = block.target_chrom_size - 1 - position_on_strand
    else:
        target_position_0based = position_on_strand

    return MappedPosition(block.target_chrom, target_position_0based, block.target_strand, block.chain_id)


def map_nearest_base_covered_by_chain(chain_index, chrom, position_0based, step, max_distance, chain_id=None):
    """Maps the base at position_0based, or if no chain block covers it, the nearest base in the direction of step that
    a chain block covers.

    Args:
        chain_index (dict): output of parse_chain_file
        chrom (str): source chromosome
        position_0based (int): 0-based position of the base in the source reference
        step (int): -1 to search leftward, 1 to search rightward
        max_distance (int): how many bases to search beyond position_0based
        chain_id (str): if specified, only accept a base covered by this chain

    Returns:
        MappedPosition: the mapped base, or None if no base within max_distance is covered (by chain_id, if specified)
    """
    for distance in range(max_distance + 1):
        mapped_position = map_position(chain_index, chrom, position_0based + step * distance)
        if mapped_position is not None and (chain_id is None or mapped_position.chain_id == chain_id):
            return mapped_position
    return None


def lift_interval(chain_index, chrom, start_0based, end_1based, max_anchor_search_distance=0):
    """Lifts a repeat interval by mapping the base on each side of it, and taking everything between them.

    Args:
        chain_index (dict): output of parse_chain_file
        chrom (str): source chromosome
        start_0based (int): 0-based start of the repeat in the source reference
        end_1based (int): 1-based end of the repeat in the source reference
        max_anchor_search_distance (int): when a chain gap covers an anchor, how far outward from the repeat to look
            for a base covered by the same chain as the other anchor

    Returns:
        2-tuple: ((target chrom, target start_0based, target end_1based, target strand), None) if the interval
            lifts over, or (None, the reason it doesn't)
    """
    left_anchor = map_position(chain_index, chrom, start_0based - 1)
    right_anchor = map_position(chain_index, chrom, end_1based)
    if left_anchor is None:
        left_anchor = map_nearest_base_covered_by_chain(
            chain_index, chrom, start_0based - 1, -1, max_anchor_search_distance,
            chain_id=right_anchor.chain_id if right_anchor else None)
    if right_anchor is None:
        right_anchor = map_nearest_base_covered_by_chain(
            chain_index, chrom, end_1based, 1, max_anchor_search_distance,
            chain_id=left_anchor.chain_id if left_anchor else None)
    if left_anchor is None or right_anchor is None:
        return None, "anchor not covered by the chain file"
    if left_anchor.chain_id != right_anchor.chain_id:
        return None, "anchors map through different chains"

    if left_anchor.strand == "+":
        target_start_0based, target_end_1based = left_anchor.position_0based + 1, right_anchor.position_0based
    else:
        target_start_0based, target_end_1based = right_anchor.position_0based + 1, left_anchor.position_0based

    if target_end_1based <= target_start_0based:
        return None, "repeat absent from the new reference"

    return (left_anchor.chrom, target_start_0based, target_end_1based, left_anchor.strand), None


def count_flanking_bases_that_continue_the_repeat(fasta, chrom, start_0based, end_1based, motif, step, max_bases):
    """Counts the bases right before (step=-1) or right after (step=1) a repeat that continue its period.

    Args:
        fasta (pysam.FastaFile): reference genome
        chrom (str): chromosome
        start_0based (int): 0-based start of the repeat
        end_1based (int): 1-based end of the repeat
        motif (str): the repeat's motif, rotated so that the repeat sequence matches it from the repeat's first base
        step (int): -1 to count bases before the repeat, 1 to count bases after it
        max_bases (int): the most bases to count

    Returns:
        int: the number of consecutive flanking bases that match the motif's period
    """
    if step == 1:
        flanking_sequence = fasta.fetch(chrom, end_1based, end_1based + max_bases).upper()
        expected_bases = (motif[(end_1based - start_0based + i) % len(motif)] for i in range(len(flanking_sequence)))
    else:
        flanking_sequence = fasta.fetch(chrom, max(0, start_0based - max_bases), start_0based).upper()[::-1]
        expected_bases = (motif[(-1 - i) % len(motif)] for i in range(len(flanking_sequence)))

    count = 0
    for base, expected_base in zip(flanking_sequence, expected_bases):
        if base != expected_base:
            break
        count += 1
    return count


def move_boundary_between_adjacent_repeats(fasta, chrom, left_start_0based, boundary, right_end_1based, left_motif,
                                           right_motif):
    """Moves the boundary between two adjacent repeats to where their sequences best match their motifs.

    The new boundary maximizes the number of bases left of it that match the left repeat's motif plus the number of
    bases right of it that match the right repeat's motif, keeping at least one motif length in each repeat. Only
    boundaries that are a whole number of motifs from the left repeat's start, or from the current boundary in the right
    repeat's motif, are considered, so that motifs that differ in only a few bases (like NOP56's GGCCTG and CGCCTG)
    aren't split mid-copy. Among equally good boundaries, the one closest to the current boundary is kept.

    Args:
        fasta (pysam.FastaFile): reference genome
        chrom (str): chromosome
        left_start_0based (int): 0-based start of the left repeat
        boundary (int): 1-based end of the left repeat, which is also the 0-based start of the right repeat
        right_end_1based (int): 1-based end of the right repeat
        left_motif (str): the left repeat's motif, rotated to match its sequence from its first base
        right_motif (str): the right repeat's motif, rotated to match its sequence from its first base

    Returns:
        int: the new boundary
    """
    sequence = fasta.fetch(chrom, left_start_0based, right_end_1based).upper()
    right_offset = boundary - left_start_0based
    # matches_before[i] and right_matches_before[i] count the bases in sequence[:i] that match each repeat's motif
    left_matches_before = [0]
    right_matches_before = [0]
    for i, base in enumerate(sequence):
        left_matches_before.append(left_matches_before[-1] + (base == left_motif[i % len(left_motif)]))
        right_matches_before.append(
            right_matches_before[-1] + (base == right_motif[(i - right_offset) % len(right_motif)]))

    def matching_bases(split):
        return left_matches_before[split] + right_matches_before[-1] - right_matches_before[split]

    possible_splits = [
        split for split in range(len(left_motif), len(sequence) - len(right_motif) + 1)
        if split % len(left_motif) == 0 or (split - right_offset) % len(right_motif) == 0
    ]
    if not possible_splits:
        return boundary
    best_split = max(possible_splits, key=lambda split: (matching_bases(split), -abs(split - right_offset)))
    return left_start_0based + best_split


def lift_record(record, chain_index, source_fasta, target_fasta, max_purity_decrease, keep_only_core_fields,
                original_locus_id_field, max_anchor_search_distance=0):
    """Lifts one catalog record over to the target reference.

    Args:
        record (dict): catalog record in ExpansionHunter format
        chain_index (dict): output of parse_chain_file
        source_fasta (pysam.FastaFile): source reference genome
        target_fasta (pysam.FastaFile): target reference genome
        max_purity_decrease (float): drop the record if any of its repeats has a purity in the target reference that
            is more than this much lower than its purity in the source reference
        keep_only_core_fields (bool): keep only the fields in CORE_FIELDS, plus MainReferenceRegion if present
        original_locus_id_field (str): name of the field that receives the record's LocusId before liftover
        max_anchor_search_distance (int): when a chain gap covers an anchor, how far outward from the repeat to look
            for a base covered by the same chain as the other anchor

    Returns:
        2-tuple: (lifted record, None), or (None, the reason the record was dropped)
    """
    reference_regions = record["ReferenceRegion"]
    is_compound_locus = isinstance(reference_regions, list)
    if not is_compound_locus:
        reference_regions = [reference_regions]
    motifs = parse_motifs_from_locus_structure(record["LocusStructure"])
    if len(motifs) != len(reference_regions):
        raise ValueError(f"Number of motifs ({len(motifs)}) != number of reference regions "
                         f"({len(reference_regions)}) in record: {record}")

    source_intervals = [parse_interval(reference_region) for reference_region in reference_regions]
    lifted_regions = []  # [target chrom, start_0based, end_1based, strand] of each repeat
    source_purities = []
    whole_copies_continue_in_source = []  # whether a whole copy of the motif continues the source repeat (left, right)
    for i, ((chrom, start_0based, end_1based), motif) in enumerate(zip(source_intervals, motifs)):
        # search every rotation, including for motifs longer than compute_best_phase_repeat_purity's default limit
        source_sequence = source_fasta.fetch(chrom, start_0based, end_1based).upper()
        source_purity, _, source_phase = compute_best_phase_repeat_purity(
            source_sequence, motif, max_motif_length_for_phase_search=len(motif))
        source_motif = motif[source_phase:] + motif[:source_phase]

        # skip the flanking bases that continue the repeat's period, but not past the previous or next repeat of a
        # locus with more than one repeat
        bases_skipped_on_the_left = count_flanking_bases_that_continue_the_repeat(
            source_fasta, chrom, start_0based, end_1based, source_motif, -1,
            max(0, start_0based - source_intervals[i - 1][2]) if i > 0 else MAX_FLANKING_BASES_TO_SKIP)
        bases_skipped_on_the_right = count_flanking_bases_that_continue_the_repeat(
            source_fasta, chrom, start_0based, end_1based, source_motif, 1,
            max(0, source_intervals[i + 1][1] - end_1based) if i < len(source_intervals) - 1
            else MAX_FLANKING_BASES_TO_SKIP)

        lifted_interval, reason = lift_interval(
            chain_index, chrom, start_0based - bases_skipped_on_the_left, end_1based + bases_skipped_on_the_right,
            max_anchor_search_distance)
        if lifted_interval is None:
            return None, reason
        target_chrom, target_start_0based, target_end_1based, target_strand = lifted_interval
        if target_strand == "-":
            bases_skipped_on_the_left, bases_skipped_on_the_right = bases_skipped_on_the_right, bases_skipped_on_the_left
        target_start_0based += bases_skipped_on_the_left
        target_end_1based -= bases_skipped_on_the_right
        if target_end_1based <= target_start_0based:
            return None, "repeat absent from the new reference"

        lifted_regions.append([target_chrom, target_start_0based, target_end_1based, target_strand])
        source_purities.append(source_purity)
        # after the swap above, left and right are in the target's orientation
        whole_copies_continue_in_source.append(
            (bases_skipped_on_the_left >= len(motif), bases_skipped_on_the_right >= len(motif)))

    if len({(chrom, strand) for chrom, _, _, strand in lifted_regions}) > 1:
        return None, "repeats map to different chromosomes or strands"

    target_strand = lifted_regions[0][3]
    has_single_repeat_structure = not is_compound_locus and record["LocusStructure"] in (
        f"({motifs[0]})*", f"({motifs[0]})+")
    if target_strand == "-":
        if is_compound_locus:
            return None, "locus with more than one repeat maps to the reverse strand"
        if not has_single_repeat_structure:
            return None, "LocusStructure has sequence outside the repeat and maps to the reverse strand"
    else:
        for (_, _, previous_end_1based, _), (_, start_0based, _, _) in zip(lifted_regions, lifted_regions[1:]):
            if start_0based < previous_end_1based:
                return None, "repeats overlap or are out of order in the new reference"

    # Since the chain's alignments allow mismatches, it can also place a length difference beyond the skipped flanking
    # bases, leaving whole copies of the motif just outside the lifted repeat. Where no whole copy continues the source
    # repeat, extend the lifted repeat over such copies, without reaching the previous or next repeat of the locus.
    motifs_on_target_strand = [reverse_complement(motif) if target_strand == "-" else motif for motif in motifs]
    target_motifs = []  # each repeat's motif rotated to match its sequence in the target reference from its first base
    for i, (target_chrom, target_start_0based, target_end_1based, _) in enumerate(lifted_regions):
        motif_length = len(motifs[i])
        _, _, phase = compute_best_phase_repeat_purity(
            target_fasta.fetch(target_chrom, target_start_0based, target_end_1based).upper(), motifs_on_target_strand[i],
            max_motif_length_for_phase_search=motif_length)
        target_motif = motifs_on_target_strand[i][phase:] + motifs_on_target_strand[i][:phase]
        target_motifs.append(target_motif)

        continues_on_the_left, continues_on_the_right = whole_copies_continue_in_source[i]
        if not continues_on_the_left:
            target_start_0based -= count_flanking_bases_that_continue_the_repeat(
                target_fasta, target_chrom, target_start_0based, target_end_1based, target_motif, -1,
                max(0, target_start_0based - lifted_regions[i - 1][2]) if i > 0 else MAX_FLANKING_BASES_TO_SKIP,
            ) // motif_length * motif_length
        if not continues_on_the_right:
            target_end_1based += count_flanking_bases_that_continue_the_repeat(
                target_fasta, target_chrom, target_start_0based, target_end_1based, target_motif, 1,
                max(0, lifted_regions[i + 1][1] - target_end_1based) if i < len(lifted_regions) - 1
                else MAX_FLANKING_BASES_TO_SKIP,
            ) // motif_length * motif_length
        lifted_regions[i][1:3] = [target_start_0based, target_end_1based]

    # The boundary between two adjacent repeats of a locus is where the chain puts it, which can be off when the chain
    # puts a length difference on the wrong side of it (like FXN, where copies of the GAA repeat land at the end of the
    # poly-A repeat before it), so move it to where the two repeats' sequences best match their motifs.
    for i in range(len(lifted_regions) - 1):
        if lifted_regions[i][2] == lifted_regions[i + 1][1]:
            lifted_regions[i][2] = lifted_regions[i + 1][1] = move_boundary_between_adjacent_repeats(
                target_fasta, lifted_regions[i][0], lifted_regions[i][1], lifted_regions[i][2],
                lifted_regions[i + 1][2], target_motifs[i], target_motifs[i + 1])

    new_motifs = []
    for i, (target_chrom, target_start_0based, target_end_1based, _) in enumerate(lifted_regions):
        motif_length = len(motifs[i])
        target_sequence = target_fasta.fetch(target_chrom, target_start_0based, target_end_1based).upper()
        if "N" in target_sequence:
            return None, "new reference sequence contains Ns"
        if len(target_sequence) < motif_length:
            return None, "repeat is shorter than one motif in the new reference"

        target_purity, _, phase = compute_best_phase_repeat_purity(
            target_sequence, motifs_on_target_strand[i], max_motif_length_for_phase_search=motif_length)
        if target_purity < source_purities[i] - max_purity_decrease:
            return None, f"purity decreased by more than {max_purity_decrease}"
        new_motifs.append(motifs_on_target_strand[i][phase:] + motifs_on_target_strand[i][:phase])

    new_motifs_iterator = iter(new_motifs)
    new_locus_structure = re.sub(r"\(([^)]*)\)", lambda match: f"({next(new_motifs_iterator)})",
                                 record["LocusStructure"])

    new_reference_regions = [f"{chrom}:{start_0based}-{end_1based}" for chrom, start_0based, end_1based, _ in lifted_regions]

    lifted_record = {
        key: value for key, value in record.items() if not keep_only_core_fields or key in CORE_FIELDS
    }
    lifted_record["LocusStructure"] = new_locus_structure
    lifted_record["ReferenceRegion"] = new_reference_regions if is_compound_locus else new_reference_regions[0]

    if "MainReferenceRegion" in record:
        if record["MainReferenceRegion"] not in reference_regions:
            raise ValueError(f"MainReferenceRegion isn't one of the record's ReferenceRegions: {record}")
        lifted_record["MainReferenceRegion"] = new_reference_regions[
            reference_regions.index(record["MainReferenceRegion"])]

    lifted_record[original_locus_id_field] = record["LocusId"]

    return lifted_record, None


def get_outer_interval(reference_regions):
    """Returns the (chrom, start_0based, end_1based) interval that spans all of a record's ReferenceRegions."""
    intervals = [parse_interval(reference_region) for reference_region in (
        reference_regions if isinstance(reference_regions, list) else [reference_regions])]
    return intervals[0][0], min(start for _, start, _ in intervals), max(end for _, _, end in intervals)


def update_coordinates_in_locus_id(lifted_record):
    """Returns the lifted record's LocusId with the coordinate-based id after its "__" replaced by the lifted
    coordinates and motif, or the LocusId unchanged if it doesn't end in a coordinate-based id after "__".

    For example, "EP400__12-132062524-132062611-CAG" becomes "EP400__12-132111361-132111451-CAG" for a record lifted to
    chr12:132111361-132111451.

    Args:
        lifted_record (dict): catalog record returned by lift_record

    Returns:
        str: the updated LocusId
    """
    match = LOCUS_ID_WITH_COORDINATES_REGEX.match(lifted_record["LocusId"])
    if not match:
        return lifted_record["LocusId"]
    if isinstance(lifted_record["ReferenceRegion"], list):
        raise ValueError(f"Can't update the coordinates in the LocusId of a locus with more than one repeat: "
                         f"{lifted_record}")

    chrom, start_0based, end_1based = parse_interval(lifted_record["ReferenceRegion"])
    motif = parse_motifs_from_locus_structure(lifted_record["LocusStructure"])[0]
    return f"{match.group(1)}{re.sub('^chr', '', chrom)}-{start_0based}-{end_1based}-{motif}"


def find_loci_that_collide_with_loci_from_elsewhere(source_intervals, target_intervals):
    """Finds lifted loci whose target interval overlaps the target interval of a locus that was on another chromosome,
    or more than MIN_DISTANCE_BETWEEN_PARALOGOUS_LOCI away, in the source reference.

    Args:
        source_intervals (list): (chrom, start_0based, end_1based) of each locus in the source reference
        target_intervals (list): (chrom, start_0based, end_1based) of each locus in the target reference, in the same
            order as source_intervals

    Returns:
        set: indices of the colliding loci
    """
    colliding_indices = set()
    overlap_candidates = []  # (target end_1based, index) of earlier loci on the current chromosome
    previous_chrom = None
    for (chrom, start_0based, end_1based), i in sorted((interval, i) for i, interval in enumerate(target_intervals)):
        if chrom != previous_chrom:
            overlap_candidates = []
            previous_chrom = chrom
        overlap_candidates = [(other_end_1based, j) for other_end_1based, j in overlap_candidates
                              if other_end_1based > start_0based]
        source_chrom, source_start_0based, source_end_1based = source_intervals[i]
        for _, j in overlap_candidates:
            other_source_chrom, other_source_start_0based, other_source_end_1based = source_intervals[j]
            if source_chrom != other_source_chrom or max(source_start_0based, other_source_start_0based) - min(
                    source_end_1based, other_source_end_1based) > MIN_DISTANCE_BETWEEN_PARALOGOUS_LOCI:
                colliding_indices.update((i, j))
        overlap_candidates.append((end_1based, i))

    return colliding_indices


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__.split("\n")[0])
    parser.add_argument("--chain-file", required=True, help="UCSC chain file from the source to the target reference, "
                        "eg. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz")
    parser.add_argument("--source-reference-fasta", required=True, help="FASTA of the reference the catalog is on")
    parser.add_argument("--target-reference-fasta", required=True, help="FASTA of the reference to lift over to")
    parser.add_argument("--max-purity-decrease", type=float, default=0.2, help="Drop a locus when the purity of any "
                        "of its repeats in the target reference is more than this much below its purity in the source "
                        "reference")
    parser.add_argument("--max-anchor-search-distance", type=int, default=20, help="When a chain gap covers the base "
                        "next to a repeat, look up to this many bases further outward for a base covered by the same "
                        "chain as the base on the other side of the repeat")
    parser.add_argument("--keep-only-core-fields", action="store_true", help="Keep only " + ", ".join(CORE_FIELDS) +
                        " and MainReferenceRegion in each record. Use this for an annotated catalog whose other fields "
                        "describe the source reference")
    parser.add_argument("--original-locus-id-field", default="LocusIdBeforeLiftover",
                        help="Name of the field that receives each record's LocusId before liftover")
    parser.add_argument("--update-coordinates-in-locus-ids", action="store_true", help="In LocusIds that end in a "
                        "coordinate-based id after '__' (eg. EP400__12-132062524-132062611-CAG), replace that id with "
                        "the lifted coordinates and motif")
    parser.add_argument("--output-dropped-loci-tsv", help="Optional path of a TSV file that lists each dropped locus "
                        "and the reason it was dropped")
    parser.add_argument("-o", "--output-catalog-json-path", required=True, help="Path of the output catalog JSON file")
    parser.add_argument("catalog_path", help="Path of the input catalog (JSON or BED)")
    args = parser.parse_args()

    print(f"Parsing {args.chain_file}")
    chain_index = parse_chain_file(args.chain_file)
    source_fasta = pysam.FastaFile(args.source_reference_fasta)
    target_fasta = pysam.FastaFile(args.target_reference_fasta)
    target_contig_order = {contig: i for i, contig in enumerate(target_fasta.references)}

    counters = collections.Counter()
    dropped_loci = []  # (LocusId, ReferenceRegion, reason) of each dropped locus
    lifted_records = []
    source_intervals = []
    for record in get_variant_catalog_iterator(args.catalog_path):
        counters["input loci"] += 1
        lifted_record, reason = lift_record(
            record, chain_index, source_fasta, target_fasta, args.max_purity_decrease, args.keep_only_core_fields,
            args.original_locus_id_field, args.max_anchor_search_distance)
        if lifted_record is None:
            dropped_loci.append((record["LocusId"], record["ReferenceRegion"], reason))
            continue

        lifted_records.append(lifted_record)
        source_intervals.append(get_outer_interval(record["ReferenceRegion"]))

    # Also keep only the first of the loci that lift to exactly the same definition, which can happen for loci that
    # were close together in the old reference
    colliding_indices = find_loci_that_collide_with_loci_from_elsewhere(
        source_intervals, [get_outer_interval(r["ReferenceRegion"]) for r in lifted_records])
    lifted_definitions = set()
    kept_records = []
    for i, lifted_record in enumerate(lifted_records):
        lifted_definition = (str(lifted_record["ReferenceRegion"]), lifted_record["LocusStructure"])
        if i in colliding_indices:
            reason = "maps onto the same repeat as a locus from elsewhere in the old reference"
        elif lifted_definition in lifted_definitions:
            reason = "has the same definition in the new reference as another lifted locus"
        else:
            lifted_definitions.add(lifted_definition)
            kept_records.append(lifted_record)
            continue
        chrom, start_0based, end_1based = source_intervals[i]
        dropped_loci.append((lifted_record[args.original_locus_id_field], f"{chrom}:{start_0based}-{end_1based}",
                             reason))
    lifted_records = kept_records

    if args.update_coordinates_in_locus_ids:
        for lifted_record in lifted_records:
            lifted_record["LocusId"] = update_coordinates_in_locus_id(lifted_record)

    counters["lifted loci"] = len(lifted_records)
    for _, _, reason in dropped_loci:
        counters[f"dropped: {reason}"] += 1

    if args.output_dropped_loci_tsv:
        with open(args.output_dropped_loci_tsv, "wt") as dropped_loci_tsv:
            dropped_loci_tsv.write("LocusId\tReferenceRegion\tReason\n")
            for locus_id, reference_region, reason in dropped_loci:
                dropped_loci_tsv.write(f"{locus_id}\t{reference_region}\t{reason}\n")

    def sort_key(lifted_record):
        reference_regions = lifted_record["ReferenceRegion"]
        chrom, start_0based, end_1based = parse_interval(
            reference_regions[0] if isinstance(reference_regions, list) else reference_regions)
        return target_contig_order[chrom], start_0based, end_1based, lifted_record["LocusStructure"], lifted_record["LocusId"]

    export_json(sorted(lifted_records, key=sort_key), args.output_catalog_json_path)

    for key, count in sorted(counters.items(), key=lambda item: -item[1]):
        print(f"{count:12,d} ({count / counters['input loci']:6.1%}) {key}")
    print(f"Wrote {len(lifted_records):,d} loci to {args.output_catalog_json_path}")


if __name__ == "__main__":
    main()
