"""Count and report the number of missing values in each field or column of a JSON or TSV file in a memory-efficient way."""

import argparse
import collections
import ijson
import gzip
import pandas as pd
import tqdm


parser = argparse.ArgumentParser()
parser.add_argument("-n", type=int, help="Number of records to process")
parser.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
parser.add_argument("json_or_tsv_path", help="Path of the JSON file to print missing values percentages for")
args = parser.parse_args()


missing_value_counters = collections.Counter()

is_json = is_tsv = False
if ".json" in args.json_or_tsv_path:
    is_json = True
elif ".tsv" in args.json_or_tsv_path:
    is_tsv = True
else:
    parser.error(f"Unsupported file format: {args.json_or_tsv_path}")

fopen = gzip.open if args.json_or_tsv_path.endswith("gz") else open
with fopen(args.json_or_tsv_path, "rt") as f:
    record_counter = 0
    if is_json:
        iterator = ijson.items(f, "item")
    else:
        header = next(f).rstrip('\n').split("\t")
        iterator = f

    if args.show_progress_bar:
        iterator = tqdm.tqdm(iterator, total=args.n, unit=" records", unit_scale=True)
    for record in iterator:
        if is_tsv:
            record = record.rstrip('\n').split("\t")
            record = dict(zip(header, record))

        for key in record.keys():
            if key not in missing_value_counters:
                print(f"Adding field: {key}")
                missing_value_counters[key] = record_counter
                    
            if record[key] is None or record[key] == "" or pd.isna(record[key]):
                missing_value_counters[key] += 1

        if is_json:
            for key in set(missing_value_counters.keys()) - set(record.keys()):
                missing_value_counters[key] += 1

        record_counter += 1
        if args.n and record_counter >= args.n:
            break
                
print(f"Missing value counters and percentages in {args.json_or_tsv_path}:")
for key, missing_count in sorted(missing_value_counters.items(), key=lambda x: x[1], reverse=True):
    print(f"{missing_count:10,d} out of {record_counter:10,d} ({missing_count / record_counter:6.1%}) records don't have, and "
          f"{record_counter - missing_count:10,d} out of {record_counter:10,d} ({(record_counter - missing_count) / record_counter:6.1%}) "
          f"records do have {key}")
