catalog=TRExplorer_v2:../../results__2025-12-30/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz
echo ================================================================================================
echo Comparing  to ${catalog}
echo ================================================================================================
set -ex

#    adotto_v1.2:adotto_tr_catalog_v1.2.bed.gz \
#    HipSTR:hg38.hipstr_reference.catalog.bed.gz \

#python3 ../compare_loci_with_catalog.py --catalog ${catalog}  --print-stats 2  --write-bed-files-with-subsets  HipSTR:hg38.hipstr_reference.catalog.bed.gz
#    Platinum_v1:human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.catalog.bed.gz \
#    PopSTR_v2:popstr_catalog_v2.bed.gz \
#    UCSC_simple_repeat_track:simpleRepeat_track_from_UCSC.bed.gz

for catalog2 in  \
  adotto_v1.2:adotto_tr_catalog_v1.2.bed.gz \
  HipSTR:hg38.hipstr_reference.catalog.bed.gz \
  Platinum_v1:human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.catalog.bed.gz \
  PopSTR_v2:popstr_catalog_v2.bed.gz \
  UCSC_simple_repeat_track:simpleRepeat_track_from_UCSC.bed.gz
do
  python3 ../compute_motif_stats.py $(echo ${catalog} | cut -f 2 -d : ) --known-loci ../loci_from_known_disease_associated_loci_v2/known_disease_associated_loci_v2.bed --known-loci functional_vntrs.bed
done

set +ex

echo "Done"
