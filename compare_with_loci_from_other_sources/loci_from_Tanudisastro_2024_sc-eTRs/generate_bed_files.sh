set -ex

cut -f 1,2,3,4 TableS1v0.1.tsv | grep -v chromosome | bedtools sort -i - | bgzip > TableS1v0.1.bed.gz;  tabix -f TableS1v0.1.bed.gz
cut -f 1,2,3,4 TableS3v0.1.tsv | grep -v chromosome | bedtools sort -i - | bgzip > TableS3v0.1.bed.gz;  tabix -f TableS3v0.1.bed.gz
cat fm_etrs_coordinates.csv | tr ',' '\t' | cut -f 1,2,3,4 | grep -v simplest_motif | bedtools sort -i - | bgzip > fm_etrs_coordinates.bed.gz;  tabix -f fm_etrs_coordinates.bed.gz
