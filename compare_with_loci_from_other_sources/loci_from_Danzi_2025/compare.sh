
tr_loci_bed=Danzi_2025_OE_or_LPSStdev_outliers.bed.gz

for catalog in TRExplorer_v2:../../results__2025-11-22/1_to_1000bp_motifs/repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz;
do
    echo ================================================================================================
    echo Comparing  to ${catalog}
    echo ================================================================================================
    set -ex
    python3 ../compare_loci_with_catalog.py --catalog ${catalog}  Danzi_2025_OE_or_LPSStdev_outliers:${tr_loci_bed} --print-stats 2  --write-bed-files-of-loci-with-same-boundardies-but-different-motifs
    set +ex
    echo
done

echo "Done"
