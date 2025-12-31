set -ex

python3 generate_vcf.py

python3 -u -m str_analysis.filter_vcf_to_tandem_repeats catalog -R ~/hg38.fa --trf-executable-path trf Gymrek_2016_eSTRs.vcf.gz --verbose

gsutil -m cp Gymrek_2016_eSTRs.tandem_repeats.bed.gz   gs://tandem-repeat-catalog/v2.0/
