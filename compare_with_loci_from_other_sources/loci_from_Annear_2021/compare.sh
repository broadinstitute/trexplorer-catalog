#catalog=Adotto_v2_1:../adotto_tr_catalog_v1.2.bed.gz
#catalog=Platinum_v1:../human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.catalog.bed.gz
#catalog=Vamos_v2_1:../vamos_catalog.ori.v2.1.bed.gz


catalog=TRExplorer_v2:../../results__2025-12-30/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz
echo ================================================================================================
echo Comparing Annear 2021 CCG loci to ${catalog}
echo ================================================================================================
set -ex
python3 ../compare_loci_with_catalog.py  --print-stats 2 --catalog ${catalog}  Annear_2021_CCG_loci:Annear_2021_loci.bed.gz --write-bed-files-with-subsets

python3 -m plot_dot_plot  -R ~/hg38.fa  loci_from_Annar_2021_CCG_loci.absent_from_TRExplorer_v2.bed.gz --output-dir plots
