#catalog=Adotto_v2_1:../adotto_tr_catalog_v1.2.bed.gz
#catalog=Platinum_v1:../human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.catalog.bed.gz
#catalog=Vamos_v2_1:../vamos_catalog.ori.v2.1.bed.gz


for catalog in \
    TRExplorer_v2:../../results__2025-11-28/1_to_1000bp_motifs/repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz \
; do
    echo ================================================================================================
    echo Comparing Tanudisastro_2024_sc-eTRs to ${catalog}
    echo ================================================================================================
    set -ex
    python3 ../compare_loci_with_catalog.py  --print-stats 2 --catalog ${catalog}  sc-eTRs_from_TableS1:TableS1v0.1.bed.gz   sc-eTRs_from_TableS3:TableS3v0.1.bed.gz  sc-eTRs-latest:fm_etrs_coordinates.bed.gz
    set +ex
    echo
done

