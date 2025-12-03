

for catalog in TRExplorer_v2:../../results__2025-11-22/1_to_1000bp_motifs/repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz;
do
    echo ================================================================================================
    echo Comparing  to ${catalog}
    echo ================================================================================================
    set -ex
    python3 ../compare_loci_with_catalog.py --catalog ${catalog}  adVNTR_gene_proximal_and_phenotype_associated_VNTRs:adVNTR_loci.bed.gz adotto_known_pathogenic_VNTRs:adotto_known_pathogenic_loci.bed.gz --print-stats 2  --write-bed-files-with-subsets
    set +ex
    echo
done

echo "Done"
