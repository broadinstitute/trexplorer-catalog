catalog=TRExplorer_v2:../../results__2025-12-24/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz
hipstr_catalog=HipSTR:hg38.hipstr_reference.catalog.bed.gz
echo ================================================================================================
echo Comparing  to ${catalog}
echo ================================================================================================
set -ex

python3 ../compare_loci_with_catalog.py --catalog ${catalog}  --print-stats 2  --write-bed-files-with-subsets ${hipstr_catalog}
python3 ../compute_motif_stats.py $(echo ${hipstr_catalog} | cut -f 2 -d : )

set +ex

echo "Done"
