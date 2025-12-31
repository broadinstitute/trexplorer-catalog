#catalog=Adotto_v2_1:../adotto_tr_catalog_v1.2.bed.gz
#catalog=Platinum_v1:../human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.catalog.bed.gz
#catalog=Vamos_v2_1:../vamos_catalog.ori.v2.1.bed.gz
#catalog=TRExplorer_v1:../../results__2025-09-04/release_draft_2025-09-04/repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz

catalog=TRExplorer_v2:../../results__2025-12-30/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz

echo ================================================================================================
echo Comparing Mukamel_2021 to ${catalog}
echo ================================================================================================
set -ex
python3 ../compare_loci_with_catalog.py Mukamel_2021:./vntrs_in_ST1.bed --print-stats 2 --catalog ${catalog}
set +ex
echo

