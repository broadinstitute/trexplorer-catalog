#TRExplorer_v2:../../results__2025-12-30/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz

for catalog in TRExplorer_v1:../../results__2024-10-01/1_to_1000bp_motifs/repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz; do
    echo ================================================================================================
    echo Comparing known disease-associated loci to ${catalog}
    echo ================================================================================================
    set -ex
    python3 ../compare_loci_with_catalog.py  --print-stats 2 --catalog ${catalog}  known_disease_associated_loci_v2:known_disease_associated_loci_v2.bed
    set +ex
    echo
done

