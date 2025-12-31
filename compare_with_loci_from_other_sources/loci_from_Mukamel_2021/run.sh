set -ex

#./compare.sh

python3 select_loci_to_include_in_catalog.py

gsutil -m cp Mukamel_2021.loci_to_include_in_catalog.bed.gz*  gs://tandem-repeat-catalog/v2.0/
