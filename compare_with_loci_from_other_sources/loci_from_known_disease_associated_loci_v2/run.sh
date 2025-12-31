set -ex

python3 generate_bed_file.py

#./compare.sh

python3 select_loci_to_include_in_catalog.py

gsutil -m cp known_disease_associated_loci_v2.loci_to_include_in_catalog.bed.gz*  gs://tandem-repeat-catalog/v2.0/
