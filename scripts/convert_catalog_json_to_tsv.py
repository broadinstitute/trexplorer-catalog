"""Convert an annotated catalog JSON file to a TSV file with the core columns first.

The VariantType field is dropped. By default, core columns that the catalog lacks are added as empty columns, so the
TSV always has the same leading columns. --omit-core-columns-missing-from-catalog leaves them out instead, which suits
a catalog that never gets some of these annotations (for example, a CHM13 catalog has no hg38 allele frequencies).

The catalog is read one record at a time, twice, since loading the whole hg38 annotated catalog at once takes about
60 GB of memory. The first pass finds the columns, in the order they first appear, and the types of the values in each.
The second pass writes the records in chunks, with each column converted to the type that pandas would give it if the
whole catalog were loaded into one table. That keeps the output the same as writing the whole table at once: for
example, a column of integers that some loci lack is written as 12.0 rather than 12.
"""

import argparse
import gzip
from pprint import pformat

import ijson
import pandas as pd

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

# The number of records converted to a table and written at a time
RECORDS_PER_CHUNK = 100_000


def iterate_catalog_records(catalog_json_path):
    """Yields the records of a catalog JSON file one at a time, with numbers read as int or float as json.load does.

    Args:
        catalog_json_path (str): path of the catalog JSON file, which can be gzipped

    Yields:
        dict: a catalog record
    """
    with (gzip.open if catalog_json_path.endswith("gz") else open)(catalog_json_path, "rb") as f:
        yield from ijson.items(f, "item", use_float=True)


def get_columns_and_their_value_types(records):
    """Returns the catalog's columns in the order they first appear, and the types of the values in each column.

    Args:
        records (iterable): catalog records (dicts)

    Returns:
        tuple: (list of column names, dict that maps each column name to a dict of {value type: an example value},
            which includes None as a key if some records lack the column or have None in it)
    """
    example_value_by_type_by_column = {}
    num_records_with_column = {}
    num_records = 0
    for record in records:
        num_records += 1
        for column, value in record.items():
            example_value_by_type = example_value_by_type_by_column.setdefault(column, {})
            example_value_by_type.setdefault(None if value is None else type(value), value)
            num_records_with_column[column] = num_records_with_column.get(column, 0) + 1

    for column, example_value_by_type in example_value_by_type_by_column.items():
        if num_records_with_column[column] < num_records:
            example_value_by_type[None] = None

    return list(example_value_by_type_by_column), example_value_by_type_by_column


def get_dtype_that_pandas_infers(example_value_by_type):
    """Returns the dtype that pd.DataFrame gives a column whose values have these types.

    Args:
        example_value_by_type (dict): {value type: an example value} for the column, with None as a key if some
            records lack the column or have None in it

    Returns:
        numpy.dtype: eg. int64 for a column of integers, float64 for a column of integers that some records lack, or
            object for a column of strings
    """
    return pd.DataFrame([{"column": value} for value in example_value_by_type.values()])["column"].dtype


def get_ordered_tsv_columns(catalog_columns, omit_core_columns_missing_from_catalog=False):
    """Returns the TSV's columns: the core columns first, followed by the catalog's other columns in their original
    order, without the VariantType column.

    Args:
        catalog_columns (list): the catalog's columns, in the order they first appear
        omit_core_columns_missing_from_catalog (bool): leave out core columns that the catalog lacks rather than
            adding them as empty columns

    Returns:
        list: column names
    """
    if omit_core_columns_missing_from_catalog:
        core_columns = [c for c in CORE_COLUMNS if c in catalog_columns]
    else:
        core_columns = CORE_COLUMNS

    return core_columns + [c for c in catalog_columns if c not in (core_columns + DROP_COLUMNS)]


def convert_records_to_table(records, tsv_columns, dtype_by_column):
    """Returns a table of these records, with the TSV's columns and each column converted to the given dtype.

    Args:
        records (list): catalog records (dicts)
        tsv_columns (list): the TSV's columns. Columns that a record lacks are empty in its row.
        dtype_by_column (dict): the dtype of each column that the catalog has. Other columns are left empty.

    Returns:
        pd.DataFrame: one row per record
    """
    # dtype=object keeps pandas from inferring each column's type from this chunk alone
    df = pd.DataFrame(records, dtype=object).reindex(columns=tsv_columns)
    for column, dtype in dtype_by_column.items():
        if column in tsv_columns:
            df[column] = df[column].astype(dtype)

    return df


def write_catalog_tsv(catalog_json_path, output_tsv_path, omit_core_columns_missing_from_catalog=False):
    """Writes the records of a catalog JSON file to a TSV file, reading the catalog one record at a time.

    Args:
        catalog_json_path (str): path of the annotated catalog JSON file, which can be gzipped
        output_tsv_path (str): path of the output TSV file, which is gzipped if it ends with gz
        omit_core_columns_missing_from_catalog (bool): leave out core columns that the catalog lacks rather than
            adding them as empty columns

    Returns:
        tuple: (number of rows written, list of the TSV's columns)
    """
    catalog_columns, example_value_by_type_by_column = get_columns_and_their_value_types(
        iterate_catalog_records(catalog_json_path))
    dtype_by_column = {column: get_dtype_that_pandas_infers(example_value_by_type)
                       for column, example_value_by_type in example_value_by_type_by_column.items()}
    tsv_columns = get_ordered_tsv_columns(catalog_columns, omit_core_columns_missing_from_catalog)

    num_rows = 0
    with (gzip.open if output_tsv_path.endswith("gz") else open)(output_tsv_path, "wt") as f:
        chunk = []
        for record in iterate_catalog_records(catalog_json_path):
            chunk.append(record)
            if len(chunk) == RECORDS_PER_CHUNK:
                convert_records_to_table(chunk, tsv_columns, dtype_by_column).to_csv(
                    f, sep="\t", index=False, header=(num_rows == 0))
                num_rows += len(chunk)
                chunk = []

        if chunk or num_rows == 0:
            convert_records_to_table(chunk, tsv_columns, dtype_by_column).to_csv(
                f, sep="\t", index=False, header=(num_rows == 0))
            num_rows += len(chunk)

    return num_rows, tsv_columns


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--omit-core-columns-missing-from-catalog", action="store_true",
                        help="Leave out core columns the catalog lacks, instead of adding them as empty columns")
    parser.add_argument("-o", "--output-tsv-path", required=True, help="Path of the output TSV file")
    parser.add_argument("catalog_json_path", help="Path of the annotated catalog JSON file")
    args = parser.parse_args()

    num_rows, tsv_columns = write_catalog_tsv(
        args.catalog_json_path, args.output_tsv_path, args.omit_core_columns_missing_from_catalog)
    print(f"Wrote {num_rows:,d} rows to {args.output_tsv_path} with columns: {pformat(tsv_columns)}")


if __name__ == "__main__":
    main()
