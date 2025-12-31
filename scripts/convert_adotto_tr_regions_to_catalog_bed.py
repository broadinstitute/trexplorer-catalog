"""
0	chr	chromosome of tandem repeat region
1	start	start position of tandem repeat region
2	end	end position of tandem repeat region
3	ovl_flag	overlap categories of annotations inside the region
4	up_buff	number of bases upstream of the first annotation's start that are non-TR sequence
5	dn_buff	number of bases downstream of the last annotation's end that are non-TR sequence
6	hom_pct	percent of region’s range annotatable as homopolymer
7	n_filtered	number of annotations removed from the region
8	n_annos	number of annotations remaining in the region
9	n_subregions	number of subregions in the region
10	mu_purity	average purity of annotations in region
11	pct_annotated	percent of the region's range (minus buffer) annotated
12	interspersed	name of interspersed repeat class found within region by RepeatMasker
13	patho	name of gene affected by a pathogenic tandem repeat in region
14	codis	name of codis site contained in region
15	gene_flag	gene features intersecting region (Ensembl v105)
16	biotype	comma separated gene biotypes intersecting region (Ensembl v105)
17	annos	JSON of TRF annotations in the region (list of dicts with keys: motif, entropy, ovl_flag, etc)
"""
import json
import gzip
import os
import tqdm

if not os.path.abspath(os.getcwd()).endswith("other"):
    os.chdir("./ref/other/")

def convert_catalog_version(catalog_version):
    input_filename = f"adotto_TRregions_v{catalog_version}.bed.gz"
    output_filename = f"adotto_tr_catalog_v{catalog_version}.bed"
    print(f"Parsing {input_filename}")
    f = gzip.open(input_filename, "rt")
    
    all_lines = f.readlines()
    output_rows = []
    for i, line in enumerate(tqdm.tqdm(all_lines, unit=" rows", unit_scale=True)):
        fields = line.strip().split("\t")
        start = int(fields[1])
        end = int(fields[2])
        #n_filtered = int(fields[7])
        n_annos = int(fields[8])
        #n_subregions = int(fields[9])
        data = json.loads(fields[17].strip('"').replace('""', '"'))
        for d in data:
            if d["start"] < start or d["end"] > end:
                print("ERROR:", i, d, start, end)
                continue
            output_rows.append((
                fields[0],
                d["start"],
                d["end"] + 1,
                d["motif"],
                ".",
            ))

    f.close()

    #%%
    print(f"Parsed {len(output_rows):,d} rows from {input_filename}")
    
    #%%
    
    with open(output_filename, "wt") as output_bed:
        for row in sorted(output_rows):
            output_bed.write("\t".join(map(str, row)) + "\n")
            
    print(f"Wrote {len(output_rows):,d} rows to {output_filename}.gz")
            
    os.system(f"bgzip -f {output_filename}")
    os.system(f"tabix -f {output_filename}.gz")
            
    print(f"Done with catalog version {catalog_version}")
#%%

# v1.0 and v1.2 appear to be identical except for some minor format changes
convert_catalog_version("1.2")
convert_catalog_version("1.0")


