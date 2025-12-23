import collections
import ijson
import intervaltree
import os
import gzip
import tqdm

from str_analysis.utils.misc_utils import parse_interval

IS_POLYMORPHIC_THRESHOLD = 0.2

stats_lookup = collections.defaultdict(dict)

counters = collections.Counter()

hipstr_catalog_path = ("./compare_catalogs/hg38.hipstr_reference.catalog.bed.gz", "HipSTR")
gangstr_catalog_path = ("./compare_catalogs/hg38_ver17.bed.gz", "GangSTR")
illumina174k_catalog_path = ("./compare_catalogs/illumina_variant_catalog.sorted.bed.gz", "Illumina174k")
trexplorer_catalog_json_path = "~/code/tandem-repeat-explorer/downloads/TR_catalog.5591917_loci.20251209_110637.json.gz"

catalog_interval_trees = collections.defaultdict(intervaltree.IntervalTree)
for path, source in hipstr_catalog_path, gangstr_catalog_path, illumina174k_catalog_path:
    print(f"Parsing {path} to interval trees")
    with gzip.open(path, "rt") as f:
        for line in tqdm.tqdm(f, unit=" lines", unit_scale=True):
            fields = line.strip().split("\t")
            chrom = fields[0].replace("chr", "")
            start_0based = int(fields[1])
            end_1based = int(fields[2])
            #motif = fields[3]
            catalog_interval_trees[chrom].add(intervaltree.Interval(start_0based, end_1based))

    print(f"Merging overlaps in {source}")
    for chrom in catalog_interval_trees:
        catalog_interval_trees[chrom].merge_overlaps(strict=False)

with gzip.open(os.path.expanduser(trexplorer_catalog_json_path), "rt") as f:  
    for i, record in tqdm.tqdm(enumerate(ijson.items(f, "item")), total=5_700_000, unit=" records", unit_scale=True):
        if "TRExplorerV1" not in record["Source"]:
            continue

        counters["total"] += 1
        motif_size = len(record["CanonicalMotif"])

        if motif_size <= 6:
            counters["total_STRs"] += 1
        else:
            counters["total_VNTRs"] += 1

        locus_id = record["LocusId"]
        chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
        chrom = chrom.replace("chr", "")
        start_0based = int(start_0based)
        end_1based = int(end_1based)

        is_polymorphic_in_HPRC256 = False
        if "HPRC256_AlleleHistogram" in record:
            is_polymorphic_in_HPRC256 = float(record["HPRC256_Stdev"]) > IS_POLYMORPHIC_THRESHOLD
            
        is_polymorphic_in_AoU1027 = False
        if "AoU1027_Stdev" in record:
            is_polymorphic_in_AoU1027 = float(record["AoU1027_Stdev"]) > IS_POLYMORPHIC_THRESHOLD


        is_in_other_catalogs = catalog_interval_trees[chrom].overlap(start_0based, end_1based)

        if is_polymorphic_in_HPRC256 or is_polymorphic_in_AoU1027: 
            counters["polymorphic"] += 1
            if motif_size <= 6:
                counters["polymorphic_STRs"] += 1
            else:
                counters["polymorphic_VNTRs"] += 1

            if not is_in_other_catalogs:
                counters["polymorphic_and_not_in_other_catalogs"] += 1
                if motif_size <= 6:
                    counters["polymorphic_STRs_and_not_in_other_catalogs"] += 1
                else:
                    counters["polymorphic_VNTRs_and_not_in_other_catalogs"] += 1

        #if i > 10_000:
        #    break

for key, count in sorted(counters.items(), key=lambda x: x[1], reverse=True):
    print(f"{count:10,d} {key}")
