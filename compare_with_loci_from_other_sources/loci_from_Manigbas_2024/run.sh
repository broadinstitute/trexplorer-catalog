set -ex

python3 generate_bed.py

#./compare.sh

#python3 select_loci_to_include_in_catalog.py

gsutil -m cp Manigbas_2024_phenotype_associated_TRs.bed.gz*  gs://tandem-repeat-catalog/v2.0/



