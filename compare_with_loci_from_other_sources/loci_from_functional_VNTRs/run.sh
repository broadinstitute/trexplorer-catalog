set -ex

./compare.sh

python3 select_loci_to_include_in_catalog.py

gsutil -m cp functional_vntrs.loci_to_include_in_catalog.bed.gz  gs://tandem-repeat-catalog/v2.0/functional_VNTRs.loci_to_include_in_catalog.bed.gz
