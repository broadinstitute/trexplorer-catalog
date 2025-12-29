set -ex


for d in \
    loci_from_adVNTR \
    loci_from_Annear_2021 \
    loci_from_clinvar_2025_11_03 \
    loci_from_functional_VNTRs \
    loci_from_Danzi_2025 \
    loci_from_known_disease_associated_loci_v2 \
    loci_from_Mukamel_2021 \
    loci_from_other_catalogs \
    loci_from_Tanudisastro_2024_sc-eTRs \
    ; do

    echo $d
    cd $d 
    ./run.sh > run.sh.log
    cd -

done


python3 generate_combined_catalog_of_selected_loci.py
