#!/usr/bin/env python3
"""
Generate a table showing counts of TRs by source, broken down into STRs and VNTRs.

STRs are defined as having motif length 1-6bp.
VNTRs are defined as having motif length 7+bp.
"""

import argparse
import gzip
import ijson
import os
from collections import defaultdict
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "json_path",
        help="Path to the JSON.gz file with annotations"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Count by source: {source: {"str": count, "vntr": count}}
    counts = defaultdict(lambda: {"str": 0, "vntr": 0})

    # Get file size for progress bar
    file_size = os.path.getsize(args.json_path)

    # Parse the JSON file using ijson for streaming
    with open(args.json_path, "rb") as raw_f:
        with tqdm(total=file_size, unit="B", unit_scale=True, desc="Parsing") as pbar:
            last_pos = 0
            with gzip.open(raw_f, "rb") as f:
                for record in ijson.items(f, "item"):
                    source = record.get("Source", "Unknown")
                    motif = record.get("Motif", "")

                    if len(motif) <= 6:
                        counts[source]["str"] += 1
                    else:
                        counts[source]["vntr"] += 1

                    # Update progress bar based on compressed bytes read
                    current_pos = raw_f.tell()
                    pbar.update(current_pos - last_pos)
                    last_pos = current_pos

    # Define explicit source order
    source_order = [
        "TRExplorerV1:KnownDiseaseAssociatedLoci",
        "TRExplorerV1:Illumina174kPolymorphicTRs",
        "TRExplorerV1:PerfectRepeatsInReference",
        "TRExplorerV1:PolymorphicTRsInT2TAssemblies",
        "TRExplorerV2:KnownDiseaseAssociatedLociV2",
        "TRExplorerV2:PolymorphicTRsInT2TAssembliesV2",
        "TRExplorerV2:KnownFunctionalVNTRs",
        "TRExplorerV2:HipSTRCatalog",
        "TRExplorerV2:AdottoTRsFromDanzi2025",
        "TRExplorerV2:ClinvarIndelsThatAreTRs2025",
        "TRExplorerV2:Tanudisastro2025",
        "TRExplorerV2:Manigbas2024",
        "TRExplorerV2:Sulovari2021",
        "TRExplorerV2:Garg2021",
        "TRExplorerV2:Mukamel2021",
        "TRExplorerV2:Annear2021",
        "TRExplorerV2:Hause2016",
        "TRExplorerV2:VamosV3",
    ]

    # Use explicit order, then append any sources not in the list
    sorted_sources = [s for s in source_order if s in counts]
    for s in sorted(counts.keys()):
        if s not in sorted_sources:
            sorted_sources.append(s)

    # Calculate totals across all sources
    total_str = sum(counts[s]["str"] for s in counts)
    total_vntr = sum(counts[s]["vntr"] for s in counts)
    grand_total = total_str + total_vntr

    # Print the table
    print(f"{'Source':<65} {'Total':>10} {'STRs':>12} {'VNTRs':>12}")
    print("-" * 100)

    prev_version = None
    for source in sorted_sources:
        # Determine version (V1, V2, or other)
        if "V1:" in source:
            curr_version = "V1"
        elif "V2:" in source:
            curr_version = "V2"
        else:
            curr_version = "other"

        # Print blank line when transitioning between versions
        if prev_version is not None and prev_version != curr_version:
            print()

        str_count = counts[source]["str"]
        vntr_count = counts[source]["vntr"]
        total = str_count + vntr_count
        print(f"{source:<65} {total:>10,} {str_count:>12,} {vntr_count:>12,}")

        prev_version = curr_version

    # Print the total row at the bottom
    print("-" * 100)
    print(f"{'TOTAL':<65} {grand_total:>10,} {total_str:>12,} {total_vntr:>12,}")


if __name__ == "__main__":
    main()
