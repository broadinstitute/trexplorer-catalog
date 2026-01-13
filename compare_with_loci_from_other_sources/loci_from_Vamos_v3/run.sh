#!/usr/bin/env bash

# This script processes Vamos v3 TR catalog data:
#   1. Converts the raw Vamos catalog to a standardized BED format
#   2. Filters loci based on quality criteria (motif size, purity, etc.)
#   3. Uploads the filtered catalog to Google Cloud Storage

set -ex

python3 convert_to_bed.py --version 3
python3 -u select_loci_to_include_in_catalog.py

gsutil -m cp vamosGenomicTR_v3.0.loci_to_include_in_catalog.bed.gz* gs://tandem-repeat-catalog/v2.0/
