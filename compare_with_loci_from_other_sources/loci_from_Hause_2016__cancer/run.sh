set -ex

python3 generate_bed.py


gcloud storage cp Hause_2016_cancer_MSI_loci.bed.gz*  gs://tandem-repeat-catalog/v2.0/
