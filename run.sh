rm nohup.out

set -ex

nohup python3 -u run_all_steps_to_generate_the_catalog.py --timestamp 2026-02-01  # --timestamp 2026-01-12

