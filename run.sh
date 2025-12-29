rm nohup.out

set -ex

nohup python3 -u run_all_steps_to_generate_the_catalog.py --skip-variation-cluster-annotations

