set -ex

./generate_bed_files.sh

#./compare.sh

python3 select_loci_to_include_in_catalog.py 

gsutil -m cp Tanudisastro_2025.loci_to_include_in_catalog.bed.gz*   gs://tandem-repeat-catalog/v2.0/
