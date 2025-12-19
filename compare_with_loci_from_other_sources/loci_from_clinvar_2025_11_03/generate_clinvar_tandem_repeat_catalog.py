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

# Determine which insertions and deletions in the ClinVar VCF represent tandem repeat expansions or contractions
cmd = f"python3 -m str_analysis.filter_vcf_to_tandem_repeats  catalog --min-tandem-repeat-length 9 --min-repeats 3 --write-tsv "
cmd += f"--trf-executable {args.trf_executable} "  # --dont-run-trf
cmd += " ".join(f"--copy-info-field-keys-to-tsv {k}" for k in ["clinsig", "stars", "clinvarid", "phenotypes"])
cmd += f" -R {args.reference_fasta} "
cmd += f"{local_clinvar_vcf_path} "
cmd += "--show-progress --verbose"

print(cmd)
#os.system(cmd)

tsv_prefix = local_clinvar_vcf_path.replace(".vcf.bgz", "") 
df = pd.read_table(f"{tsv_prefix}.tandem_repeats.tsv.gz")
df["DetectionMode"] = df["DetectionMode"].str.replace("merged:", "")
df["DetectedUsingTRF"] = (df["DetectionMode"] == "trf")


output_bed_path = f"{tsv_prefix}.tandem_repeats.excluding_TRF.bed"
df_without_trf = df[~df["DetectedUsingTRF"]]
df_without_trf[["Chrom", "Start0Based", "End1Based", "Motif"]].to_csv(output_bed_path, sep="\t", header=False, index=False)
os.system(f"bgzip -f {output_bed_path}")
os.system(f"tabix -f {output_bed_path}.gz")
print(f"Wrote {len(df_without_trf):,d} rows to {output_bed_path}.gz")

for detected_using_trf in [False, True]:
    print("-" * 120)
    summary_df = df[df["DetectedUsingTRF"] == detected_using_trf]
    total = len(summary_df)
    summary_df = summary_df[["clinsig", "stars"]].value_counts().reset_index().sort_values(by=["clinsig", "stars"], ascending=False)
    summary_df.rename({"count": "TR loci"}, axis=1, inplace=True)
    if detected_using_trf:
        print(f"Expansion/contractions at {total:,d} interrupted TR loci that were identified by running TRF on variants in the ClinVar VCF:")
    else:
        print(f"Expansion/contractions at {total:,d} pure or mostly pure TR loci that were identified by scanning variants in the ClinVar VCF:")
    print(summary_df.to_string(index=False))


        
