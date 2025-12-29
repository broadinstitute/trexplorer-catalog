set -ex

python3 generate_clinvar_tandem_repeat_catalog.py
./compare.sh
python3 select_loci_to_include_in_catalog.py

gsutil -m cp clinvar_2025_11_03.merged.tandem_repeats.detailed.loci_to_include_in_catalog.bed.gz gs://tandem-repeat-catalog/v2.0/clinvar_2025_11_03.loci_to_include_in_catalog.bed.gz
