"""Convert an annotated catalog JSON file to a TSV file with the core columns first.

The VariantType field is dropped. By default, core columns that the catalog lacks are added as empty columns, so the
TSV always has the same leading columns. --omit-core-columns-missing-from-catalog leaves them out instead, which suits
a catalog that never gets some of these annotations (for example, a CHM13 catalog has no hg38 allele frequencies).
"""

import argparse
import gzip
import json
import pandas as pd
from pprint import pformat

CORE_COLUMNS = [
    'LocusId', 'ReferenceRegion', 'LocusStructure', 'CanonicalMotif', 'TRsInRegion',
    'Source', 'GencodeGeneRegion', 'GencodeGeneId', 'GencodeGeneName', 'GencodeTranscriptId',
    'RefseqGeneRegion', 'RefseqGeneId', 'RefseqGeneName', 'RefseqTranscriptId',
    'ManeGeneRegion', 'ManeGeneId', 'ManeGeneName', 'ManeTranscriptId',
    'KnownDiseaseAssociatedMotif',  'KnownDiseaseAssociatedLocus', 'NsInFlanks',
    'LeftFlankMappability', 'FlanksAndLocusMappability', 'RightFlankMappability',
    'FoundInTRExplorerV1:KnownDiseaseAssociatedLoci', 'FoundInTRExplorerV1:Illumina174kPolymorphicTRs',
    'FoundInTRExplorerV1:PerfectRepeatsInReference', 'FoundInTRExplorerV1:PolymorphicTRsInT2TAssemblies',
    'NumRepeatsInReference', 'ReferenceRepeatPurity',
    'AlleleFrequenciesFromIllumina174k', 'StdevFromIllumina174k',
    'AlleleFrequenciesFromT2TAssemblies', 'StdevFromT2TAssemblies',
    'VariationCluster', 'VariationClusterMotifs', 'VariationClusterSizeDiff',
    #'AlleleFrequenciesFromHPRC256', 'StdevFromHPRC256',
    'AoU1027_Stdev', 'AoU1027_Median', 'AoU1027_MaxAllele',
]

DROP_COLUMNS = ['VariantType', ]


def order_catalog_columns(df, omit_core_columns_missing_from_catalog=False):
    """Returns df with the core columns first, followed by the other columns in their original order.

    Args:
        df (pd.DataFrame): catalog table, one row per locus
        omit_core_columns_missing_from_catalog (bool): leave out core columns that df lacks rather than adding them
            as empty columns

    Returns:
        pd.DataFrame: the reordered table, without the VariantType column
    """
    if omit_core_columns_missing_from_catalog:
        core_columns = [c for c in CORE_COLUMNS if c in df.columns]
    else:
        core_columns = CORE_COLUMNS
        for c in set(core_columns) - set(df.columns):
            df[c] = None

    return df[core_columns + [c for c in df.columns if c not in (core_columns + DROP_COLUMNS)]]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--omit-core-columns-missing-from-catalog", action="store_true",
                        help="Leave out core columns the catalog lacks, instead of adding them as empty columns")
    parser.add_argument("-o", "--output-tsv-path", required=True, help="Path of the output TSV file")
    parser.add_argument("catalog_json_path", help="Path of the annotated catalog JSON file")
    args = parser.parse_args()

    with gzip.open(args.catalog_json_path, "rt") as f:
        data = json.load(f)
    print(f"Loaded {len(data):,d} records from {args.catalog_json_path}")
    print("Writing it to TSV file...")
    df = order_catalog_columns(pd.DataFrame(data), args.omit_core_columns_missing_from_catalog)
    df.to_csv(args.output_tsv_path, sep="\t", index=False)
    print(f"Wrote {len(df):,d} rows to {args.output_tsv_path} with columns: {pformat(list(df.columns))}")


if __name__ == "__main__":
    main()
