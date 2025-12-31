set -ex

#./compare.sh

python3 -u select_loci_to_include_in_catalog.py

gsutil -m cp hg38.hipstr_reference.loci_to_include_in_catalog.bed.gz*  gs://tandem-repeat-catalog/v2.0/
