set -e

# loci_from_adVNTR \
# loci_from_Annear_2021 \
#    loci_from_known_disease_associated_loci_v2 \
#    loci_from_other_catalogs \
#
	
for d in \
    ./loci_from_Hause_2016__cancer \
    ./loci_from_Danzi_2025 \
    ./loci_from_Manigbas_2024 \
    ./loci_from_clinvar_2025_11_03 \
    ./loci_from_Gymrek_2016__eSTRs \
    ./loci_from_known_disease_associated_loci_v2 \
    ./loci_from_Sulovari_2019__human-specific_STRs \
    ./loci_from_Tanudisastro_2024__sc-eTRs \
    ./loci_from_functional_VNTRs \
    ./loci_from_Garg_2021__eVNTRs_mVNTRs \
    ./loci_from_Annear_2021 \
    ./loci_from_HipSTR_catalog \
    ./loci_from_Mukamel_2021 \
; do

    echo =========================================
    echo $d
    cd $d

    set -x
    ls run.sh
    grep gsutil run.sh

    ./run.sh 2>&1 | tee run.sh.log 
    set +x
    cd ..

done

echo Done

#python3 generate_combined_catalog_of_selected_loci.py
