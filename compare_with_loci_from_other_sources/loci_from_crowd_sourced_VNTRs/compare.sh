

label="$(cat functional_vntrs.bed | wc -l)_known_VNTRs"

for catalog in TRExplorer_v2:../../results__2025-11-22/1_to_1000bp_motifs/repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz; do
#for catalog in \
#    Vamos_v2_1:../vamos_catalog.ori.v2.1.bed.gz \
#    Adotto_v2_1:../adotto_tr_catalog_v1.2.bed.gz \
#    Platinum_v1:../human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.catalog.bed.gz \
#    TRExplorer_v1:../../results__2025-09-04/release_draft_2025-09-04/repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz \
#    TRExplorer_v2:../../results__2025-11-22/1_to_1000bp_motifs/repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz \
#; do
    echo ================================================================================================
    echo Comparing ${label} to ${catalog}
    echo ================================================================================================
    set -ex
    python3 ../compare_loci_with_catalog.py ${label}:functional_vntrs.bed --print-stats 3 --catalog ${catalog}
    set +ex
    echo
done

echo "Done"
