#catalog=Adotto_v2_1:../adotto_tr_catalog_v1.2.bed.gz
#catalog=Platinum_v1:../human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.catalog.bed.gz
#catalog=Vamos_v2_1:../vamos_catalog.ori.v2.1.bed.gz


for catalog in \
    TRExplorer_v2:../../results__2025-12-24/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz \
    TRExplorer_v1:../../results__2025-09-04/release_draft_2025-09-04/repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz \
; do
    echo ================================================================================================
    echo Comparing clinvar to ${catalog}
    echo ================================================================================================
    set -ex
    python3 ../compare_loci_with_catalog.py  --print-stats 2 --catalog ${catalog}   clinvar_TRs:./clinvar_2025_11_03.tandem_repeats.bed.gz  clinvar_TRs_excluding_TRF:./clinvar_2025_11_03.tandem_repeats.excluding_TRF.bed.gz 
    set +ex
    echo
done

