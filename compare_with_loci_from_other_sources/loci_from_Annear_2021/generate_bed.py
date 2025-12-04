from liftover import get_lifter   # https://pypi.org/project/liftover/
import os
import pandas as pd
from pprint import pprint

os.chdir(os.path.dirname(__file__))


converter = get_lifter('hg19', 'hg38', one_based=False)

df = pd.read_table("./41598_2021_82050_MOESM2_ESM.txt")
print(f"Read {len(df):,d} loci from ./41598_2021_82050_MOESM2_ESM.txt")
#pprint(df.iloc[0].to_dict())

output_rows = []
for _, row in df.iterrows():
    chrom_hg19 = row["Chr"]
    start_hg19_0based = row["Start"] - 1

    liftover_result = converter[chrom_hg19][start_hg19_0based]
    if len(liftover_result) != 1 or len(liftover_result[0]) != 3:
        print(f"WARNING: Liftover result for {chrom_hg19}:{start_hg19_0based} has {len(liftover_result)} elements: {liftover_result}. Skipping...")
        continue
    chrom_hg38, start_hg38_0based, strand_hg38 = liftover_result[0]
    end_hg38_1based = start_hg38_0based + row["Reference.Repeat.Length"] * 3
    motif = "CCG"
    output_rows.append((chrom_hg38, start_hg38_0based, end_hg38_1based, motif))


output_bed_path = "Annear_2021_loci.bed"
with open(output_bed_path, "wt") as output_bed:
    for output_row in sorted(output_rows):
        chrom_hg38, start_hg38_0based, end_hg38_1based, motif = output_row
        output_bed.write(f"{chrom_hg38}\t{start_hg38_0based}\t{end_hg38_1based}\t{motif}\n")

os.system(f"bgzip -f {output_bed_path}")
os.system(f"tabix -f {output_bed_path}.gz")
print(f"Wrote {len(output_rows):,d} loci to {output_bed_path}.gz")
"""
Example row:

{'Chr': '1',
 'Start': 10057405,
 'Cohort.Polymorphism': 0.0,
 'Gene': 'RBP7',
 'Maximum.Repeat.Length': 5,
 'Median.Repeat.Length': 5,
 'Minimum.Repeat.Length': 5,
 'Reference.Repeat.Length': 4,
 'Region': 'intronic',
 'SD': 0.0,
 'Variance': 0.0}
"""

