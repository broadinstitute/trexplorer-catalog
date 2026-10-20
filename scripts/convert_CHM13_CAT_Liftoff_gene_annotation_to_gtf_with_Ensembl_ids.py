"""Convert the T2T consortium's CHM13 gene annotation (GENCODE v35 mapped onto CHM13 by CAT and Liftoff) from GFF3 to a
GTF whose gene and transcript ids are the Ensembl ids of the GENCODE genes and transcripts.

The GFF3 (chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz) gives each gene and transcript its own id: CHM13_G... and
CHM13_T... for CAT, LOFF_G... and LOFF_T... for Liftoff. It keeps the GENCODE gene and transcript that each one came
from in its source_gene and source_transcript attributes, and the GENCODE gene name in source_gene_common_name. UCSC's
catLiftOffGenesV1.gtf.gz is made from the same file, but keeps only the CHM13 ids. This script writes a GTF in the
format that str_analysis.annotate_and_filter_str_catalog reads, so that the catalog's GencodeGeneId and
GencodeTranscriptId fields hold ENSG and ENST ids, as in the hg38 catalog:
- Transcripts with an ENST source_transcript and an ENSG source_gene get those ids and the source_gene_common_name gene
  name. A GENCODE transcript that maps to several places in CHM13 (eg. the paralogous gene copies that Liftoff added)
  keeps the same ENST id at each place.
- The copies of the rRNA genes in the rDNA arrays (RNA45S, RNA18S, RNA5-8S and RNA28S, whose transcript biotype is
  rRNA) have no Ensembl ids, so they keep their LOFF_G... and LOFF_T... ids. This keeps the loci in the rDNA arrays
  annotated.
- All other transcripts are left out: the novel isoforms that CAT predicted for GENCODE genes, which have an ENSG
  source_gene but no ENST source_transcript, and the CAT predictions that have neither.

The GFF3's CDS records include the stop codon, as is usual in GFF3, while GENCODE's GTF ends the CDS before the stop
codon and counts the stop codon as part of the 3' UTR. So the script removes each transcript's stop codon (from its
stop_codon records) from its CDS records. CAT transcripts have exon and CDS records but no UTR records, so the GTF then
gets a UTR record for each part of a coding transcript's exons that lies outside its CDS, which, as in GENCODE's GTF,
includes the stop codon. Each UTR record is written as a 5UTR or 3UTR record, labeled from its position relative to its
own copy's CDS. Plain UTR records would leave the labeling to str_analysis, which compares a UTR with the first CDS
record that has the same transcript id, and so with another copy's CDS when several copies share an ENST id. The GFF3's
own UTR records, which only some Liftoff transcripts have, are replaced by these. Its gene, intron, start_codon and
stop_codon records aren't written to the GTF, since str_analysis doesn't use them.
"""

import argparse
import collections
import gzip
import re
import urllib.parse

# characters that str_analysis's GTF parser can't read inside an attribute value
INVALID_GTF_ATTRIBUTE_VALUE_REGEX = re.compile(r'[\s;"]')


def parse_gff3_attributes(attributes_column):
    """Parses the attributes column of a GFF3 record.

    Args:
        attributes_column (str): the 9th column, eg. "ID=CHM13_T0000001;gene_name=AC114498.1"

    Returns:
        dict: attribute name => value, with GFF3's percent-encoded characters (eg. %2C for a comma) decoded
    """
    attributes = {}
    for key_and_value in attributes_column.strip().strip(";").split(";"):
        key, value = key_and_value.split("=", 1)
        attributes[key] = urllib.parse.unquote(value)
    return attributes


def get_gtf_ids_of_transcript(transcript_attributes):
    """Picks the ids and names to write to the GTF for a transcript, or None for a transcript to leave out.

    Args:
        transcript_attributes (dict): attributes of the GFF3 transcript record, as returned by parse_gff3_attributes

    Returns:
        dict: gene_id, transcript_id, gene_name, gene_type and transcript_type for the GTF, or None if the transcript
            has no Ensembl ids and isn't an rRNA gene copy
    """
    source_gene = transcript_attributes.get("source_gene", "")
    source_transcript = transcript_attributes.get("source_transcript", "")
    if source_gene.startswith("ENSG") and source_transcript.startswith("ENST"):
        gtf_ids = {
            "gene_id": source_gene,
            "transcript_id": source_transcript,
            "gene_name": transcript_attributes["source_gene_common_name"],
        }
    elif transcript_attributes.get("transcript_biotype") == "rRNA":
        gtf_ids = {
            "gene_id": transcript_attributes["gene_id"],
            "transcript_id": transcript_attributes["transcript_id"],
            "gene_name": transcript_attributes["gene_name"],
        }
    else:
        return None

    gtf_ids["gene_type"] = transcript_attributes["biotype"]
    gtf_ids["transcript_type"] = transcript_attributes["transcript_biotype"]
    for key, value in gtf_ids.items():
        if not value or INVALID_GTF_ATTRIBUTE_VALUE_REGEX.search(value):
            raise ValueError(f"Can't write {key} value '{value}' of transcript {transcript_attributes['ID']} to a GTF")

    return gtf_ids


def remove_stop_codons_from_CDS_intervals(CDS_intervals, stop_codon_intervals):
    """Removes the stop codon bases from a transcript's GFF3 CDS records, which include them, so that the CDS ends
    before the stop codon as in GENCODE's GTF.

    Args:
        CDS_intervals (list): (start_1based, end_1based) of each GFF3 CDS record of the transcript
        stop_codon_intervals (list): (start_1based, end_1based) of each of the transcript's stop_codon records. A stop
            codon split by an intron has two records.

    Returns:
        list: (start_1based, end_1based) of each remaining part of the CDS, sorted by start
    """
    remaining_CDS_intervals = sorted(CDS_intervals)
    for stop_codon_start_1based, stop_codon_end_1based in stop_codon_intervals:
        CDS_intervals_without_this_stop_codon = []
        for CDS_start_1based, CDS_end_1based in remaining_CDS_intervals:
            if CDS_start_1based < stop_codon_start_1based:
                CDS_intervals_without_this_stop_codon.append(
                    (CDS_start_1based, min(CDS_end_1based, stop_codon_start_1based - 1)))
            if CDS_end_1based > stop_codon_end_1based:
                CDS_intervals_without_this_stop_codon.append(
                    (max(CDS_start_1based, stop_codon_end_1based + 1), CDS_end_1based))
        remaining_CDS_intervals = CDS_intervals_without_this_stop_codon

    return remaining_CDS_intervals


def compute_5UTR_and_3UTR_intervals(exon_intervals, CDS_intervals, strand):
    """Returns the parts of a transcript's exons that lie outside its CDS, each labeled as a 5' or 3' UTR.

    Given a CDS that ends before the stop codon, as returned by remove_stop_codons_from_CDS_intervals, these include
    the stop codon, like GENCODE's UTR records.

    Args:
        exon_intervals (list): (start_1based, end_1based) of each exon
        CDS_intervals (list): (start_1based, end_1based) of each part of the CDS
        strand (str): the transcript's strand, "+" or "-"

    Returns:
        list: (feature_type, start_1based, end_1based) of each UTR part, where feature_type is the GTF feature type
            "5UTR" or "3UTR", sorted by start, or an empty list for a non-coding transcript
    """
    if not CDS_intervals:
        return []

    if strand not in ("+", "-"):
        raise ValueError(f"Can't tell the 5' UTR from the 3' UTR of a coding transcript on strand '{strand}'")

    UTR_feature_type_before_CDS, UTR_feature_type_after_CDS = ("5UTR", "3UTR") if strand == "+" else ("3UTR", "5UTR")
    CDS_start_1based = min(start_1based for start_1based, _ in CDS_intervals)
    CDS_end_1based = max(end_1based for _, end_1based in CDS_intervals)
    UTR_intervals = []
    for exon_start_1based, exon_end_1based in sorted(exon_intervals):
        if exon_start_1based < CDS_start_1based:
            UTR_intervals.append(
                (UTR_feature_type_before_CDS, exon_start_1based, min(exon_end_1based, CDS_start_1based - 1)))
        if exon_end_1based > CDS_end_1based:
            UTR_intervals.append(
                (UTR_feature_type_after_CDS, max(exon_start_1based, CDS_end_1based + 1), exon_end_1based))

    return UTR_intervals


def format_gtf_line(transcript_record, feature_type, start_1based, end_1based, gtf_ids):
    """Returns one line of the output GTF.

    Args:
        transcript_record (dict): chrom, source and strand of the transcript
        feature_type (str): transcript, exon, CDS, 5UTR or 3UTR
        start_1based (int): start of the record
        end_1based (int): end of the record
        gtf_ids (dict): attributes returned by get_gtf_ids_of_transcript

    Returns:
        str: the GTF line, ending with a newline
    """
    attributes_column = " ".join(f'{key} "{value}";' for key, value in gtf_ids.items())
    return "\t".join([
        transcript_record["chrom"], transcript_record["source"], feature_type, str(start_1based), str(end_1based), ".",
        transcript_record["strand"], ".", attributes_column]) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].replace("\n", " "))
    parser.add_argument("-o", "--output-gtf-path", required=True, help="Path of the output GTF, eg. *.gtf.gz")
    parser.add_argument("gff3_path", help="Path of chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz")
    args = parser.parse_args()

    # gather each transcript's record, exons, CDS and stop codon, keyed by the GFF3's own transcript id, which is unique
    # (unlike the ENST id of a GENCODE transcript that maps to several places)
    transcripts = {}
    exon_intervals_by_transcript_id = collections.defaultdict(list)
    CDS_intervals_by_transcript_id = collections.defaultdict(list)
    stop_codon_intervals_by_transcript_id = collections.defaultdict(list)
    with (gzip.open if args.gff3_path.endswith("gz") else open)(args.gff3_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            feature_type = fields[2]
            if feature_type not in ("transcript", "exon", "CDS", "stop_codon"):
                continue

            attributes = parse_gff3_attributes(fields[8])
            interval = (int(fields[3]), int(fields[4]))
            if feature_type == "transcript":
                transcripts[attributes["transcript_id"]] = {
                    "chrom": fields[0], "source": fields[1], "interval": interval, "strand": fields[6],
                    "attributes": attributes,
                }
            elif feature_type == "exon":
                exon_intervals_by_transcript_id[attributes["transcript_id"]].append(interval)
            elif feature_type == "CDS":
                CDS_intervals_by_transcript_id[attributes["transcript_id"]].append(interval)
            else:
                stop_codon_intervals_by_transcript_id[attributes["transcript_id"]].append(interval)

    transcript_ids_without_a_transcript_record = (
        set(exon_intervals_by_transcript_id) | set(CDS_intervals_by_transcript_id)) - set(transcripts)
    if transcript_ids_without_a_transcript_record:
        raise ValueError(f"{len(transcript_ids_without_a_transcript_record):,d} transcript ids have exon or CDS records "
                         f"but no transcript record, eg. {sorted(transcript_ids_without_a_transcript_record)[:3]}")

    counters = collections.Counter()
    with (gzip.open if args.output_gtf_path.endswith("gz") else open)(args.output_gtf_path, "wt") as f:
        for transcript_id, transcript_record in transcripts.items():
            attributes = transcript_record["attributes"]
            gtf_ids = get_gtf_ids_of_transcript(attributes)
            if gtf_ids is None:
                if attributes.get("source_gene", "").startswith("ENSG"):
                    counters["transcripts left out: novel isoforms of GENCODE genes, with no ENST id"] += 1
                else:
                    counters["transcripts left out: no ENSG or ENST id, and not an rRNA gene copy"] += 1
                continue

            if gtf_ids["transcript_id"].startswith("ENST"):
                counters["transcripts written with their Ensembl ids"] += 1
            else:
                counters["transcripts written with their LOFF ids: rRNA gene copies"] += 1

            CDS_intervals_without_stop_codon = remove_stop_codons_from_CDS_intervals(
                CDS_intervals_by_transcript_id[transcript_id], stop_codon_intervals_by_transcript_id[transcript_id])
            f.write(format_gtf_line(transcript_record, "transcript", *transcript_record["interval"], gtf_ids))
            for feature_type, start_1based, end_1based in (
                [("exon", *interval) for interval in sorted(exon_intervals_by_transcript_id[transcript_id])]
                + [("CDS", *interval) for interval in CDS_intervals_without_stop_codon]
                + compute_5UTR_and_3UTR_intervals(
                    exon_intervals_by_transcript_id[transcript_id], CDS_intervals_without_stop_codon,
                    transcript_record["strand"])
            ):
                f.write(format_gtf_line(transcript_record, feature_type, start_1based, end_1based, gtf_ids))
                counters[f"{feature_type} records written"] += 1

    print(f"Wrote {args.output_gtf_path}:")
    for description, count in sorted(counters.items()):
        print(f"{count:12,d}  {description}")


if __name__ == "__main__":
    main()
