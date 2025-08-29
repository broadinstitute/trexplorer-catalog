import argparse
import collections
import gzip
import os.path


def _compute_introns_for_chrom(chrom, transcript_to_exons):
    """Converts transcript information to a list of introns.

    Args:
        chrom (str): chromosome name
        transcript_to_exons (dict): maps (transcript_id, gene_id, gene_name) keys to
            (exon_number, start_1based, end_1based, strand) values.

    Return:
        list: 0-based intron coordinates for the given transcripts
    """
    introns_0based = []
    for exons_in_transcript in transcript_to_exons.values():
        prev_exon_3prime = None
        for exon_number, start_1based, end_1based, strand in sorted(exons_in_transcript):
            if strand == "+":
                current_exon_5prime = int(start_1based)
                if exon_number > 1:
                    assert prev_exon_3prime < current_exon_5prime
                    intron_start_1based = prev_exon_3prime + 1
                    intron_end_1based = current_exon_5prime - 1
                prev_exon_3prime = int(end_1based)
            elif strand == "-":
                current_exon_5prime = int(end_1based)
                if exon_number > 1:
                    assert prev_exon_3prime > current_exon_5prime
                    intron_start_1based = current_exon_5prime + 1
                    intron_end_1based = prev_exon_3prime - 1
                prev_exon_3prime = int(start_1based)
            else:
                print(f"ERROR: strand == {strand}")
                continue

            if exon_number > 1 and intron_start_1based > intron_end_1based:
                print(f"WARNING: intron_start_1based > intron_end_1based: {chrom}:{intron_start_1based}-{intron_end_1based}")

            if exon_number > 1:
                introns_0based.append((chrom, intron_start_1based - 1, intron_end_1based))

    return introns_0based


def parse_exons_and_introns_from_gff(gff_path):
    """Returns a set of (chrom, start, end) 3-tuples where each tuple represents 0-based coordinates for all
    introns in the given gff file.

    Args:
        gff_path (str): Path of a gencode gff file

    Return:
          2-tuple:
            set: 3-tuples representing 0-based exon coordinates
            set: 3-tuples representing 0-based intron coordinates
    """
    if ".gff" not in gff_path:
        print(f"WARNING: filename doesn't contain '.gff': {gff_path}")

    print(f"Parsing: {gff_path}")
    # example line: chr1	HAVANA	exon	11869	12227	.	+	.	ID=exon:ENST00000456328.2:1;Parent=ENST00000456328.2;gene_id=ENSG00000223972.5;transcript_id=ENST00000456328.2;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;transcript_type=processed_transcript;transcript_name=DDX11L1-002;exon_number=1;exon_id=ENSE00002234944.1;level=2;transcript_support_level=1;tag=basic;havana_gene=OTTHUMG00000000961.2;havana_transcript=OTTHUMT00000362751.1
    exons_0based = []
    introns_0based = []

    fopen = gzip.open if gff_path.endswith("gz") else open
    with fopen(gff_path, "rt") as f:
        prev_chrom = None
        transcript_to_exons = collections.defaultdict(list)
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "exon":
                continue
            chrom = fields[0]
            if "_" in chrom:
                continue

            if chrom != prev_chrom:
                if prev_chrom is not None:
                    introns_for_prev_chrom = _compute_introns_for_chrom(prev_chrom, transcript_to_exons)
                    print(f"Parsed {len(introns_for_prev_chrom)} introns from {prev_chrom}")
                    introns_0based.extend(introns_for_prev_chrom)
                transcript_to_exons = collections.defaultdict(list)
                prev_chrom = chrom

            start_1based = int(fields[3])
            end_1based = int(fields[4])
            exons_0based.append((chrom, start_1based - 1, end_1based))

            strand = fields[6]
            annot = {}
            for key_value in fields[8].strip("; ").split(";"):
                if "=" in key_value:
                    key, value = key_value.split("=", 1)
                else:
                    key = key_value
                    value = None
                annot[key] = value

            missing_keys = [f'"{k}"' for k in ("exon_number", "transcript_id", "gene_id", "gene_name") if k not in annot]
            if missing_keys:
                raise ValueError(f"Missing keys {', '.join(missing_keys)} in gff line '{line}'. \nAre you sure {gff_path} is a Gencode .gff file?")

            exon_number = annot['exon_number']
            transcript_id = annot['transcript_id'].split(".")[0]
            gene_id = annot['gene_id'].split(".")[0]
            gene_name = annot['gene_name']
            transcript_to_exons[(transcript_id, gene_id, gene_name)].append((int(exon_number), start_1based, end_1based, strand))

        introns_for_prev_chrom = _compute_introns_for_chrom(prev_chrom, transcript_to_exons)
        print(f"Found {len(introns_for_prev_chrom)} introns on {prev_chrom}")
        introns_0based.extend(introns_for_prev_chrom)

    return set(exons_0based), set(introns_0based)


def get_sum_of_all_interval_lengths(bed_path):
    """Returns the sum of the sizes of all intervals in the given BED file."""

    fopen = gzip.open if bed_path.endswith("gz") else open
    with fopen(bed_path, "rt") as f:
        total = 0
        for line in f:
            fields = line.rstrip("\n").split("\t")
            start_0based = int(fields[1])
            end_1based = int(fields[2])
            total += end_1based - start_0based

    return total


def merge_and_write_to_bed(output_path, chrom_start_0based_end_list):
    """Writes the given list of (chrom, start_0based, end_1based) tuples to a bed file at output_path
    and then bgzips and tabixes the file.

    Args:
        output_path (str): Path to write the BED file to. Must end with .bed.
        chrom_start_0based_end_list (list): List of (chrom, start_0based, end_1based) tuples to write to the BED file.

    Returns:
        str: Path to the bgzipped BED file (f"{output_path}.gz")

    """
    if not output_path.endswith(".bed"):
        raise ValueError(f"Output path must end with .bed: {output_path}")

    print(f"Writing {len(chrom_start_0based_end_list):,d} records to: {output_path}")

    with open(output_path, "wt") as f:

        for chrom, start_0based, end_1based in sorted(chrom_start_0based_end_list):
            f.write(f"{chrom}\t{start_0based}\t{end_1based}\n")

    os.system(f"bedtools merge -i {output_path} > {output_path}.merged")
    os.system(f"mv {output_path}.merged {output_path}")

    os.system(f"bgzip -f {output_path}")
    os.system(f"tabix -f {output_path}.gz")

    return f"{output_path}.gz"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-k", "--keep-bed-files", action="store_true", help="Don't delete intermediate .bed files after computing statistics")
    p.add_argument("gff_path", help="Path of GFF or GFF3 file from Gencode, MANE, or other gene models")
    args = p.parse_args()

    # parse input GFF file into lists of exons and introns
    filename_prefix = os.path.basename(args.gff_path.split(".gff")[0])
    exons_0based, introns_0based = parse_exons_and_introns_from_gff(args.gff_path)

    # output exons and introns to separate BED files
    introns_bed_path = f"{filename_prefix}.introns.bed"
    exons_bed_path = f"{filename_prefix}.exons.bed"

    merge_and_write_to_bed(introns_bed_path, introns_0based)
    merge_and_write_to_bed(exons_bed_path, exons_0based)

    # compute 3rd BED file with intronic regions that do not overlap exons
    os.system(f"bedtools subtract -a {introns_bed_path}.gz -b {exons_bed_path}.gz | sort -k1,1 -k2,2n | bgzip > {filename_prefix}.introns_minus_exons.bed.gz")
    os.system(f"tabix -f {filename_prefix}.introns_minus_exons.bed.gz")
    print(f"Wrote {filename_prefix}.introns_minus_exons.bed.gz")

    # compute and print statistics
    total_exon_bases = get_sum_of_all_interval_lengths(f"{exons_bed_path}.gz")
    total_intron_bases = get_sum_of_all_interval_lengths(f"{filename_prefix}.introns_minus_exons.bed.gz")

    genome_size = 3*10**9
    print(f"Total exon bases: {total_exon_bases:,d} out of {genome_size:,d} ({total_exon_bases/genome_size:.2%} of genome)")
    print(f"Total intron bases (excluding exons): {total_intron_bases:,d} out of {genome_size:,d} ({total_intron_bases/genome_size:.2%} of genome)")

    # delete intermediate BED files unless user requested to keep them
    if not args.keep_bed_files:
        for bed_path in [f"{filename_prefix}.introns_minus_exons.bed.gz", f"{introns_bed_path}.gz", f"{exons_bed_path}.gz"]:
            os.remove(bed_path)
            os.remove(f"{bed_path}.tbi")


if __name__ == "__main__":
    main()
