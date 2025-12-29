import argparse
import os
import pandas as pd

from annotation_utils.get_clinvar_table import export_clinvar_vcf  # import https://github.com/bw2/annotation-utils

p = argparse.ArgumentParser()
p.add_argument("-R", "--reference-fasta", default="~/hg38.fa")
p.add_argument("--trf-executable", default="trf")
args = p.parse_args()

if not os.path.isfile(os.path.expanduser(args.reference_fasta)):
    p.error(f"{args.reference_fasta} not found")

# Export all ClinVar variants to a VCF with the following info field keys: clinsig, stars, clinvarid, phenotypes.
# For example: clinsig=Pathogenic;stars=1;clinvarid=2994835;phenotypes=Peroxisome biogenesis disorder, complementation group 7
local_clinvar_vcf_path = export_clinvar_vcf(only_pathogenic=True, include_phenotypes=True)
#local_clinvar_vcf_path = "clinvar_2025_11_03.vcf.bgz"

def run(cmd):
    print(cmd)
    os.system(cmd)

local_clinvar_vcf_prefix = local_clinvar_vcf_path.replace('.vcf.bgz', '')
run(f"gunzip -c {local_clinvar_vcf_path} > {local_clinvar_vcf_prefix}.vcf")
run(f"rm {local_clinvar_vcf_path}.tbi")

run(f"bgzip -f {local_clinvar_vcf_prefix}.vcf")
run(f"tabix -f {local_clinvar_vcf_prefix}.vcf.gz")
    

# Determine which insertions and deletions in the ClinVar VCF represent tandem repeat expansions or contractions
cmd = f"python3 -m str_analysis.filter_vcf_to_tandem_repeats catalog --min-tandem-repeat-length 9 --min-repeats 3 --write-detailed --trf-min-repeats-in-reference 2 --trf-min-purity 0.66 "
cmd += f"--trf-executable {args.trf_executable} "  # --dont-run-trf
cmd += " ".join(f"--copy-info-field-keys-to-tsv {k}" for k in ["clinsig", "stars", "clinvarid", "phenotypes"])
cmd += f" -R {args.reference_fasta} "
cmd += f"{local_clinvar_vcf_prefix}.vcf.gz "
cmd += "--show-progress --verbose"
run(cmd)

cmd2 = f"python3 -m str_analysis.filter_vcf_to_tandem_repeats merge --write-detailed -R {args.reference_fasta} --show-progress --verbose {local_clinvar_vcf_prefix}.tandem_repeats.detailed.bed.gz"
run(cmd2)


run("rm -rf trf_working_dir")

