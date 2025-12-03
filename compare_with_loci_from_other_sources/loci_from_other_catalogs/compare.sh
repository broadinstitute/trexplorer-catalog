
tr_loci_bed=Danzi_2025_OE_or_LPSStdev_outliers.bed.gz

for catalog in TRExplorer_v2:../../results__2025-11-22/1_to_1000bp_motifs/repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz;
do
    echo ================================================================================================
    echo Comparing  to ${catalog}
    echo ================================================================================================
    set -ex
    python3 ../compare_loci_with_catalog.py --catalog ${catalog}  --print-stats 2  --write-bed-files-with-subsets \
	    adotto_v1.2:adotto_tr_catalog_v1.2.bed.gz \
	    HipSTR:hg38.hipstr_reference.catalog.bed.gz \
	    Platinum_v1:human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.catalog.bed.gz \
	    PopSTR_v2:popstr_catalog_v2.bed.gz \
	    UCSC_simple_repeat_track:simpleRepeat_track_from_UCSC.bed.gz
    set +ex
    echo
done

echo "Done"
