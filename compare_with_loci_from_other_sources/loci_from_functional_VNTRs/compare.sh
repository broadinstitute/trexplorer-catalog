

label="$(cat functional_vntrs.bed | wc -l)_known_VNTRs"

#for catalog in TRExplorer_v2:../../results__2025-12-07/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz; do

    #Platinum_v1:../human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.catalog.bed.gz \
    #Adotto_v2_1:../adotto_tr_catalog_v1.2.bed.gz \
    #Vamos_v2_1:../vamos_catalog.ori.v2.1.bed.gz \

#TRExplorer_v1:../../results__2025-09-04/release_draft_2025-09-04/repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz \
#TRExplorer_v2:../../results__2026-01-01/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz
#for catalog in TRExplorer_v2:../../results__2026-01-01/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.final_merged.bed.gz;

    #TRExplorer_v1:../../results__2026-01-01/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.merged_v1.json.gz \
    #TRExplorer_v2a:../../results__2026-01-01/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.merged_v2a.json.gz \
    #TRExplorer_v2b:../../results__2026-01-01/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.merged_v2b.json.gz \

for catalog in TRExplorer_v2:../..//results__2026-01-01/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.final_merged.json.gz;
do
    echo ================================================================================================
    echo Comparing ${label} to ${catalog}
    echo ================================================================================================
    set -ex
    python3 ../compare_loci_with_catalog.py ${label}:functional_vntrs.bed --print-stats 3  --overlap-by-motif-max-x 30 --catalog ${catalog}
    
    #python3 ../compute_motif_stats.py $(echo ${catalog} | cut -f 2 -d : ) --known-loci ../loci_from_known_disease_associated_loci_v2/known_disease_associated_loci_v2.bed --known-loci functional_vntrs_v2.bed
    set +ex
    echo
done

echo "Done"
