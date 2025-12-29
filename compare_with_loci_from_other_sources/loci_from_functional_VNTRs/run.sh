set -ex

./compare.sh

python3 select_loci_to_include_in_catalog.py
