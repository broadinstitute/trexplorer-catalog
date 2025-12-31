set -ex

python3 generate_bed.py

gsutil -m cp Sulovari_2019_human_specific_STRs.bed.gz*  gs://tandem-repeat-catalog/v2.0/
