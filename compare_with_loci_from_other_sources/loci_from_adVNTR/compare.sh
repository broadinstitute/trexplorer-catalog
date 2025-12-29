catalog=TRExplorer_v2:../../results__2025-12-28/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz

echo ================================================================================================
echo Comparing  to ${catalog}
echo ================================================================================================
set -ex
python3 ../compare_loci_with_catalog.py --catalog ${catalog}  adVNTR_gene_proximal_and_phenotype_associated_VNTRs:adVNTR_loci.bed.gz adotto_known_pathogenic_VNTRs:adotto_known_pathogenic_loci.bed.gz --print-stats 2  --write-bed-files-with-subsets
set +ex
echo

echo "Done"
