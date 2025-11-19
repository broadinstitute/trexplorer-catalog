import argparse
import os

from annotation_utils.get_clinvar_table import export_clinvar_vcf  # import https://github.com/bw2/annotation-utils

p = argparse.ArgumentParser()
p.add_argument("-R", "--reference-fasta", default="~/hg38.fa")
args = p.parse_args()

if not os.path.isfile(os.path.expanduser(args.reference_fasta)):
    p.error(f"{args.reference_fasta} not found")

# Export all ClinVar variants to a VCF with the following info field keys: clinsig, stars, clinvarid, phenotypes.
# For example: clinsig=Pathogenic;stars=1;clinvarid=2994835;phenotypes=Peroxisome biogenesis disorder, complementation group 7
local_path = export_clinvar_vcf(only_pathogenic=True, include_phenotypes=True)

# Determine which insertions and deletions in the ClinVar VCF represent tandem repeat expansions or contractions
cmd = f"python3 -m str_analysis.filter_vcf_to_tandem_repeats  catalog --min-tandem-repeat-length 9 --min-repeats 3 "
cmd += " ".join(f"--copy-info-field-keys-to-tsv {k}" for k in ["clinsig", "stars", "clinvarid", "phenotypes"])
cmd += f" -R {args.reference_fasta} "
cmd += f"{local_clinvar_vcf_path} "
cmd += "--show-progress --verbose"
print(cmd)
os.system(cmd)
