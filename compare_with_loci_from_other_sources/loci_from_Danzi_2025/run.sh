

python3 generate_bed_file_from_table.py
./compare.sh
python3 select_loci_to_include_in_catalog.py

gsutil -m cp Danzi_2025_OE_or_LPSStdev_outliers.narrow_adotto_boundaries.only_matching_motif_lengths.loci_to_include_in_catalog.bed.gz  gs://tandem-repeat-catalog/v2.0/Danzi_2025.adotto.loci_to_include_in_catalog.bed.gz
