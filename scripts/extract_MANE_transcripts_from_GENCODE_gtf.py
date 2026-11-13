"""Write the records of the MANE transcripts in a GENCODE GTF file to a new GTF file.

MANE (https://www.ncbi.nlm.nih.gov/refseq/MANE/) only publishes GTF files for GRCh38, so the hg19 catalog pipeline makes
its MANE gene annotations from GENCODE's GRCh37 mapping of the GENCODE release that the hg38 catalog uses (eg.
gencode.v50lift37.annotation.gtf.gz), which tags each MANE transcript's records with "MANE_Select" or
"MANE_Plus_Clinical", as the GRCh38 MANE GTF has them. This keeps the comment lines and every transcript, exon, CDS, UTR,
start codon and stop codon record with one of those tags. The gene records, which don't have tags, are left out, as
str_analysis.annotate_and_filter_str_catalog skips gene records anyway.
"""

import argparse
import collections
import gzip

MANE_TAGS = ('tag "MANE_Select"', 'tag "MANE_Plus_Clinical"')


def is_MANE_transcript_record(gtf_line):
    """Returns True if a GTF line is a record of a MANE transcript.

    Args:
        gtf_line (str): a line of a GENCODE GTF file that isn't a comment

    Returns:
        bool: True if the line has the MANE_Select or MANE_Plus_Clinical tag. Each tag in MANE_TAGS ends with its
            closing quote, so a longer tag like MANE_Select_candidate doesn't match.
    """
    return any(mane_tag in gtf_line for mane_tag in MANE_TAGS)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("-o", "--output-gtf-path", required=True, help="Path of the output GTF file, which is gzipped "
                        "if it ends in .gz")
    parser.add_argument("gencode_gtf_path", help="Path of the GENCODE GTF file, which can be gzipped")
    args = parser.parse_args()

    counters = collections.Counter()
    with (gzip.open if args.gencode_gtf_path.endswith("gz") else open)(args.gencode_gtf_path, "rt") as input_file, \
            (gzip.open if args.output_gtf_path.endswith("gz") else open)(args.output_gtf_path, "wt") as output_file:
        for line in input_file:
            if line.startswith("#"):
                output_file.write(line)
                continue
            if is_MANE_transcript_record(line):
                output_file.write(line)
                counters[line.split("\t")[2]] += 1

    print(f"Wrote {sum(counters.values()):,d} MANE transcript records to {args.output_gtf_path}:")
    for feature_type, count in sorted(counters.items(), key=lambda x: -x[1]):
        print(f"{count:12,d}  {feature_type}")


if __name__ == "__main__":
    main()
