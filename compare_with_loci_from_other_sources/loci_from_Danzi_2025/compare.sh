
catalog=TRExplorer_v2:../../results__2025-11-22/1_to_1000bp_motifs/repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz

echo ================================================================================================
echo Comparing  to ${catalog}
echo ================================================================================================
set -ex
python3 ../compare_loci_with_catalog.py --catalog ${catalog}  \
  Danzi_2025_outliers__wide_boundaries:Danzi_2025_OE_or_LPSStdev_outliers.wide_boundaries.bed.gz \
  Danzi_2025_outliers__narrow_boundaries_all:Danzi_2025_OE_or_LPSStdev_outliers.narrow_adotto_boundaries.all_overlapping_loci.bed.gz \
  Danzi_2025_outliers__narrow_boundaries_only_matching_motif_lengths:Danzi_2025_OE_or_LPSStdev_outliers.narrow_adotto_boundaries.only_matching_motif_lengths.bed.gz \
  Danzi_2025_outliers__narrow_boundaries_only_matching_motif_lengths_replaced_with_LPS_motif:Danzi_2025_OE_or_LPSStdev_outliers.narrow_adotto_boundaries.only_matching_motif_lengths_replaced_with_LPS_motif.bed.gz \
  adotto_v1.2:adotto_tr_catalog_v1.2.bed.gz \
  --print-stats 2 \
  --write-bed-files-with-subsets

