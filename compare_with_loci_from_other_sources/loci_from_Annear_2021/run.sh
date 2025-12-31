set -ex

python3 generate_bed.py

#./compare.sh

python3 select_loci_to_include_in_catalog.py

gsutil -m cp Annear_2021.loci_to_include_in_catalog.bed.gz*  gs://tandem-repeat-catalog/v2.0/



