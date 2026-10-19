"""Write a copy of an annotated catalog that keeps only the fields ExpansionHunter reads.

Each record keeps LocusId, ReferenceRegion, VariantType and LocusStructure. chrM loci are skipped because
ExpansionHunter prints an error like 'Unable to extract chrM:-793-207 from hg38.fa' for them, and so are loci with
more than 5 Ns in their flanks.
"""

import argparse
import gzip
import ijson
import json

EXPANSION_HUNTER_FIELDS = {"LocusId", "ReferenceRegion", "VariantType", "LocusStructure"}
MAX_NS_IN_FLANKS = 5


def is_chrM_locus(record):
    """Returns True if the record's LocusId says it's on the mitochondrial chromosome."""
    return record["LocusId"].startswith("M-") or record["LocusId"].startswith("chrM-")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("-o", "--output-catalog-json-path", required=True,
                        help="Path of the output ExpansionHunter catalog JSON file")
    parser.add_argument("annotated_catalog_json_path", help="Path of the annotated catalog JSON file. Its records "
                        "must have the NsInFlanks field")
    args = parser.parse_args()

    f = gzip.open(args.annotated_catalog_json_path, "rt")
    out = gzip.open(args.output_catalog_json_path, "wt")
    i = 0
    total = 0
    out.write("[")
    for record in ijson.items(f, "item", use_float=True):
        total += 1
        if is_chrM_locus(record):
            print(f"Skipping chrM locus: {record['LocusId']}")
            continue
        if record["NsInFlanks"] > MAX_NS_IN_FLANKS:
            continue
        if i > 0:
            out.write(", ")
        i += 1
        out.write(json.dumps({
            k: v for k, v in record.items() if k in EXPANSION_HUNTER_FIELDS
        }, indent=4))
    out.write("]")
    out.close()
    f.close()
    print(f"Wrote {i:,d} out of {total:,d} ({i/total:0.1%}) loci to {args.output_catalog_json_path}")


if __name__ == "__main__":
    main()
