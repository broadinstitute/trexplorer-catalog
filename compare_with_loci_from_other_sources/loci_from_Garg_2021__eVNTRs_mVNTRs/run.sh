set -ex

python3 generate_bed.py

gcloud storage cp Garg_2021_eVNTRs_mVNTRs.bed.gz*  gs://tandem-repeat-catalog/v2.0/
