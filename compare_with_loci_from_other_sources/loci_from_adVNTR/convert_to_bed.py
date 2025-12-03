#The TSV export code below is equivalent to:
#   sqlite3 -header -separator $'\t' p_vntrs_g_vntrs_v1.db "SELECT * FROM vntrs;" > vntrs.tsv

import collections
import os
import pandas as pd
import sqlite3

os.chdir(os.path.dirname(__file__))

 
# Source 1:
# One of the adVNTR papers
# (Variable number tandem repeats mediate the expression of proximal genes: Supplementary Data 1):
# https://www.nature.com/articles/s41467-021-22206-z#Sec29
# lists 10,264 VNTR loci originally identified by TRF (with unclear parameters) as described in the above paper.

df = pd.read_excel("41467_2021_22206_MOESM3_ESM.xlsx")
df[["chromosome", "start", "end"]] = df["VNTR Coordinates (GRCh38)"].str.split(r"[:-]", expand=True)
df["start_0based"] = df["start"].astype(int)
df["end_1based"] = df["end"].astype(int)
df.rename(columns={"First Repeating Pattern": "motif"}, inplace=True)

loci = set()
output_rows = []
for _, row in df.iterrows():
  loci.add((row["chromosome"], row["start_0based"], row["end_1based"]))
  output_rows.append((row["chromosome"], row["start_0based"], row["end_1based"], row["motif"]))


# Source 2:
# This SQLite database of gene-proximal and phenotype-associated VNTRs was downloaded from https://github.com/mehrdadbakhtiari/adVNTR
conn = sqlite3.connect("p_vntrs_g_vntrs_v1.db")
cursor = conn.cursor()

cursor.execute("SELECT * FROM vntrs;")
rows = cursor.fetchall()

# write to TSV
tsv_path = "vntrs.tsv"
with open(tsv_path, "wt") as f:
    f.write("\t".join(map(lambda x: x[0], cursor.description)) + "\n")
    for row in rows:
        f.write("\t".join(map(str, row)) + "\n")

df = pd.read_table(tsv_path)
#print(df.head())

"""
  id  nonoverlapping chromosome  ref_start  ...                                      left_flanking                                     right_flanking                                            repeats scaled_score
0   3            True       chr1      11166  ...  CGCCGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAG...  TTATAGGGAAACACCCGGAGCATATGCTGTTTGGTCTCAGTAGACT...  CGCCCCTTGCTTGCAGCCGGGCACTACAGGACCCGCTTGCTCACGG...          0.0
1  28            True       chr1     120967  ...  AGGCAATTTATCAAAGTCCCCTAATCCTCCAAAATCGCTATTTTTT...  TCATTTCCAAATTCCCCAGCGTTCATATTTGTCAGTGCAAGTAAAG...  TTATATA,TTATA,TATCTA,TTATATA,TAATATA,TATCTA,TT...          0.0
2  35            True       chr1     136192  ...  GCTTAGGGAAGTTGTGGGCCTACCAGGGCCGCTGGGAGCTGGGCAG...  CATGAGTTGGGCATCAACAGGCCACCGTGAGGGAGGAGCTGGGCCG...  GGCCTGTTGAGGCAGGGGGTCACGCTGACCTCTGTCCGCGTGGGAG...          0.0
3  39            True       chr1     138160  ...  TCGCCGGGAGGCCCAACCTTGGCGTGGAGGAGCCCACCGACCGGAG...  GACACCATCTGGGTCTGGAGGGTCCACTGTGAGGCAGAGGCTGACC...  GAGGCAGAGGCTGGGCCTGTGCAGGCCTTCGG,GAGGCAGGAGGCT...          0.0
4  41            True       chr1     138821  ...  GGCCTGGAGAGGCTGCCGAAAGGCAGGAGCTTCACCTGAGGATGCC...  CACCGTGAGGCATAAGCTGGATGTAGAGAGGCCAGTGTGAGGCAAG...  TGAGGCAGCAGTTGTGCCTGTAGACCCAGCCA,AGAGGAAGAGGTG...          0.0
"""

# convert to BED
df = df[["chromosome", "ref_start", "repeats"]]
for _, row in df.iterrows():
    chrom = row["chromosome"]
    start_0based = row["ref_start"]
    repeats = row["repeats"].split(",")
    # get the motif that has the most common length
    motif_length_counts = collections.Counter([len(r) for r in repeats])
    most_common_length = motif_length_counts.most_common(1)[0][0]
    most_common_motif = collections.Counter([r for r in repeats if len(r) == most_common_length]).most_common(1)[0][0]
    repeat_sequence = "".join(repeats)
    end_1based = start_0based + len(repeat_sequence)

    if (chrom, start_0based, end_1based) in loci:
      print(f"{len(loci):,d}: Duplicate locus: {chrom}:{start_0based}-{end_1based} {most_common_motif}")
      continue

    #print(f"{len(loci):,d}: New locus: {chrom}:{start_0based}-{end_1based} {most_common_motif}")
    #print([len(r) for r in repeats], most_common_length, most_common_motif)
    #print(f"{chrom}:{start}-{end}  {most_common_motif}  {repeat_sequence}")
    loci.add((chrom, start_0based, end_1based))
    output_rows.append((chrom, start_0based, end_1based, most_common_motif))

# write to BED
output_bed_filename = "adVNTR_loci.bed"
df_bed = pd.DataFrame(output_rows, columns=["chrom", "start_0based", "end_1based", "motif"])
df_bed.sort_values(by=["chrom", "start_0based", "end_1based"], inplace=True)
df_bed.to_csv(output_bed_filename, sep="\t", index=False, header=False)
os.system(f"bgzip -f {output_bed_filename}")
os.system(f"tabix -f {output_bed_filename}.gz")
print(f"Wrote {len(df_bed):,d} loci to {output_bed_filename}.gz")


# Source 1:
# The Adotto catalog lists adVNTR as one of its 9 source catalogs
# (Analysis and benchmarking of small and large genomic variants across tandem repeats: Supp. Table 1):
# https://www.nature.com/articles/s41587-024-02225-z

# Also, its Supp. Table 14 lists known pathogenic STR and VNTR loci (119 VNTRs with motifs)
df = pd.read_excel("41587_2024_2225_MOESM3_ESM.xlsx", 3, skiprows=1)
df = df.rename(columns={"Motifs 1": "motif"})
df = df[df.motif.notna() & (df.motif.str.len() > 6)]  # leaves 119 rows

output_bed_filename = "adotto_known_pathogenic_loci.bed"
df = df[["chrom", "start", "end", "motif"]]
df.sort_values(by=["chrom", "start", "end"], inplace=True)
df.to_csv(output_bed_filename, sep="\t", index=False, header=False)
os.system(f"bgzip -f {output_bed_filename}")
os.system(f"tabix -f {output_bed_filename}.gz")
print(f"Wrote {len(df):,d} loci to {output_bed_filename}.gz")

