import collections
import json
import intervaltree
import gzip
import os

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

os.chdir(os.path.dirname(__file__))

SOURCES = [
    ("Annear2021", "loci_from_Annear_2021/loci_from_Annar_2021_CCG_loci.absent_from_TRExplorer_v2.bed.gz"),
    #("KnownVNTRs", "loci_from_crowd_sourced_VNTRs/functional_vntrs.loci_to_include_in_catalog.bed.gz"),
    ("ClinVar:2025-11-03", "loci_from_clinvar_2025_11_03/clinvar_2025_11_03.tandem_repeats.loci_to_include_in_catalog.bed.gz"),
    ("Mukamel2021", "loci_from_Mukamel_2021/vntrs_in_ST1.loci_to_include_in_catalog.bed.gz"),
    ("Tanudisastro2024:sc-eTRs", "loci_from_Tanudisastro_2024_sc-eTRs/TableS1v0.1.loci_to_include_in_catalog.bed.gz"),
    ("KnownDiseaseAssociatedLoci_v2", "loci_from_known_disease_associated_loci_v2/known_disease_associated_loci_v2.loci_to_include_in_catalog.bed.gz"),
    #("Danzi2025", "loci_from_Danzi_2025/Danzi_2025.loci_to_include_in_catalog.bed.gz"),
]


result_catalog = collections.defaultdict(intervaltree.IntervalTree)
for source_label, source_catalog_path in SOURCES:
    print(f"Processing {source_label} from {source_catalog_path}")
    total_counter = 0
    added_counter = 0
    with gzip.open(source_catalog_path, "rt") as f:
        # parse the bed file and add to the output catalog as long as it doesn't match an existing locus (boundaries differ by < 2 repeats)
        for line in f:
            total_counter += 1
            fields = line.strip().split("\t")
            chrom = fields[0].replace("chr", "")
            start_0based = int(fields[1])
            end_1based = int(fields[2])
            motif = fields[3]
            canonical_motif = compute_canonical_motif(motif)

            # make sure the definition is not the same as any existing locus (ie. boundaries differ by < 2 repeats)
            overlapping_loci = result_catalog[chrom].overlap(start_0based, end_1based)
            matches_existing_locus = False
            for overlapping_interval in overlapping_loci:
                intersection = min(overlapping_interval.end, end_1based) - max(overlapping_interval.begin, start_0based)
                union = max(overlapping_interval.end, end_1based) - min(overlapping_interval.begin, start_0based)
                overlapping_canonical_motif = overlapping_interval.data["canonical_motif"]
                if union - intersection <= 2*len(motif) and canonical_motif == overlapping_canonical_motif:
                    matches_existing_locus = True
                    break

            if matches_existing_locus:
                print(f"Skipping {source_label} locus {chrom}:{start_0based}-{end_1based} {motif} because it matches an existing locus definition")
                continue
            
            added_counter += 1
            result_catalog[chrom].add(intervaltree.Interval(start_0based, end_1based, data={
                "motif": motif,
                "canonical_motif": canonical_motif,
                "source": source_label,
            }))
    print(f"Added {added_counter:,d} out of {total_counter:,d} loci from {source_label}")

total = 0
output_path = "combined_catalog_of_selected_loci_from_other_sources.json.gz"
output_records = []
for chrom in result_catalog:
    for interval in result_catalog[chrom]:
        total += 1
        output_records.append({
            "LocusId": f"{chrom}-{interval.begin}-{interval.end}-{interval.data['motif']}",
            "ReferenceRegion": f"chr{chrom}:{interval.begin}-{interval.end}",
            "LocusStructure": "("+interval.data["motif"]+")*",
            "VariantType": "Repeat",
            "CanonicalMotif": interval.data["canonical_motif"],
            "Source": interval.data["source"],
        })
with gzip.open(output_path, "wt") as f:
    json.dump(output_records, f, indent=4)

print(f"Wrote {len(output_records):,d} loci to {output_path}")
