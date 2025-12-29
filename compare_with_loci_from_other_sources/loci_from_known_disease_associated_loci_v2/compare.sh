#catalog=Adotto_v2_1:../adotto_tr_catalog_v1.2.bed.gz
#catalog=Platinum_v1:../human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.catalog.bed.gz

for catalog in Vamos_v2_1:../vamos_catalog.ori.v2.1.bed.gz \
		   TRExplorer_v2:../../results__2025-12-28/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz; do
    echo ================================================================================================
    echo Comparing known disease-associated loci to ${catalog}
    echo ================================================================================================
    set -ex
    python3 ../compare_loci_with_catalog.py  --print-stats 2 --catalog ${catalog}  known_disease_associated_loci_v2:known_disease_associated_loci_v2.bed
    set +ex
    echo
done

