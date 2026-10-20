"""Generate the CHM13 (T2T-CHM13v2.0) version of the TRExplorer repeat catalog.

Following https://github.com/broadinstitute/trexplorer-catalog/issues/4, each source only adds loci that aren't
redundant with the loci of the sources before it:
1. known disease-associated loci: the str-analysis catalog that gnomAD uses, lifted over from hg38, and then
   STRchive's CHM13 catalog, keeping only the STRchive loci whose disease association is confident (Definitive,
   Strong, Moderate or Limited evidence) or newly reported
2. perfect repeats in CHM13 that span at least 9bp and 3 or more repeats (same as for the hg38 catalog)
3. the TRExplorer v2.1 hg38 catalog, lifted over to CHM13 with scripts/lift_over_catalog.py
4. not yet included: TR expansions and contractions found by running dipcall on T2T assemblies aligned to CHM13

None of the str-analysis catalog's known disease-associated loci are deduplicated: both a locus's primary definition
and its alternative definitions (like EP400 and EP400__12-132111361-132111451-CAG) are kept, even though they overlap.
The liftover replaces the hg38 coordinates in the alternative definitions' LocusIds with their CHM13 coordinates.

The KnownDiseaseAssociatedLocus labels come from the lifted str-analysis catalog, plus STRchive's CHM13 loci for the
known disease-associated loci that didn't lift over.

The hg38ToHs1 chain maps part of hg38's chrX PAR1 onto CHM13's chrY PAR1, which is all Ns in chm13v2.0.maskedY.fasta,
so the liftover drops those hg38 loci (about 400 of the hg38 catalog's 11,000 chrX PAR loci). The CHM13 perfect
repeats source still covers the repeats of CHM13's chrX PARs directly.

The hg38-only annotations (allele frequencies, TenK10K, HPRC256, AoU and variation clusters) aren't added. Mappability
and MANE annotations are skipped for now, since there's no CHM13 mappability track or MANE annotation yet. The Gencode
gene annotations come from the T2T consortium's CAT/Liftoff mapping of GENCODE v35 genes onto CHM13, whose gene and
transcript ids are LOFF_G... and LOFF_T... rather than ENSG... and ENST... ids. That GTF has no UTR features, so
GencodeGeneRegion reports a repeat in a 5' or 3' UTR as 'exon'; the RefseqGeneRegion field does distinguish UTRs.

The CentromereAndSatelliteRegion field holds the class of the centromeric satellite region (from the Miga lab's CenSat
annotation) that each locus overlaps, such as HSat3 or ct. About 1.3 million loci lie in CHM13's satellite and rDNA
arrays, most of them in HSat2 and HSat3 arrays that hg38 lacks, which is why this catalog has more loci than the hg38
one.
"""

import argparse
import datetime
import functools
import json
import os
import time
from types import SimpleNamespace as NS

# install str-analysis python package. This has to run before the str_analysis imports below, or a missing package
# fails the import before it can be installed.
if os.system("""
if ! python3 -m pip show str-analysis > /dev/null 2>&1
then
    python3 -m pip install --upgrade git+https://github.com/broadinstitute/str-analysis.git
fi
""") != 0:
    raise RuntimeError("Failed to install the str-analysis python package")

from str_analysis.utils.file_utils import tee_stdout_and_stderr_to_log_file

# imported after the install block above, since str-analysis is what installs the requests package this module uses
import run_all_steps_helpers
from run_all_steps_helpers import chdir, is_bgzipped, tabix_coordinate_flags


VERSION = "2.1"
MIN_MOTIF_SIZE = 1
MAX_MOTIF_SIZE = 1000

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--chm13-reference-fasta", default="chm13v2.0.maskedY.fasta",
                    help="Local path of the CHM13 reference genome FASTA file, which must have a .fai index. Use one "
                         "with hard-masked chrY PARs, such as a local copy of "
                         "gs://gcp-public-data--broad-references/t2t/v2/chm13v2.0.maskedY.fasta, to match the CHM13 "
                         "perfect repeats source, which excludes loci in the chrY PARs")
parser.add_argument("--hg38-reference-fasta", default="hg38.fa",
                    help="Local path of the hg38 reference genome FASTA file, which must have a .fai index. The "
                         "liftover steps use it to compare each locus's purity before and after liftover")
parser.add_argument("--hg38-catalog",
                    default=f"https://storage.googleapis.com/tandem-repeat-catalog/v{VERSION}/"
                            f"TRExplorer.repeat_catalog_v{VERSION}.hg38.1_to_1000bp_motifs.EH.json.gz",
                    help="Local path, gs:// path or URL of the TRExplorer hg38 catalog to lift over to CHM13")
parser.add_argument("--output-prefix", default=f"TRExplorer.repeat_catalog_v{VERSION}.chm13")
parser.add_argument("--only-step", type=int, help="Only run this one step")
parser.add_argument("--start-with-step", type=int, help="Start with a specific step number")
parser.add_argument("--end-with-step", type=int, help="End with a specific step number")
parser.add_argument("--timestamp", default=datetime.datetime.now().strftime('%Y-%m-%d'),
                    help="Timestamp to use in the output directory name")
parser.add_argument("--dry-run", action="store_true", help="Print commands without running them")

args = parser.parse_args()

run = functools.partial(run_all_steps_helpers.run, args=args)
download_unless_local_copy_is_up_to_date = functools.partial(
    run_all_steps_helpers.download_unless_local_copy_is_up_to_date, args=args)

print("TIMESTAMP:", args.timestamp)

base_dir = os.path.abspath(".")

# Remote inputs are downloaded here and shared by every run, so a file is only downloaded again after its version on
# the server changes. This is a subfolder of the hg38 pipeline's downloaded_data folder so that the two pipelines
# never overwrite each other's copy of a file that has the same name in both, like the known disease loci catalog.
downloaded_data_dir = os.path.join(base_dir, "downloaded_data", "chm13")
run(f"mkdir -p {downloaded_data_dir}")

for key in "chm13_reference_fasta", "hg38_reference_fasta":
    if not os.path.isfile(getattr(args, key)):
        parser.error(f"{key} must be a local file, and {getattr(args, key)} isn't one")
    setattr(args, key, os.path.abspath(getattr(args, key)))

# a local --hg38-catalog path is relative to the launch directory, so resolve it before changing to the results folder
if os.path.isfile(args.hg38_catalog):
    args.hg38_catalog = os.path.abspath(args.hg38_catalog)

# The working directory name includes the genome so that an hg38 and a CHM13 build started on the same day don't
# write to the same folder.
working_dir = os.path.abspath(f"results__chm13__{args.timestamp}")

run(f"mkdir -p {working_dir}")
chdir(working_dir)

release_draft_folder = os.path.abspath(f"release_draft_{args.timestamp}")
run(f"mkdir -p {release_draft_folder}")

start_datetime = datetime.datetime.now()

# Name the log after the run's start time so that a dry run or a partial rerun doesn't overwrite the
# log of an earlier full run.
log_file_path = f"{release_draft_folder}/{args.output_prefix}.all_steps.{start_datetime.strftime('%Y%m%d-%H%M%S')}.log"
print(f"Writing to log file: {log_file_path}")
tee_stdout_and_stderr_to_log_file(log_file_path)

script_start_time = time.time()
print(f"Started at: {start_datetime.strftime('%Y-%m-%d %H:%M:%S')}")

# Remote files that the steps read. They're downloaded here so that none of the steps need internet access.
remote_inputs = [
    # The same version of the str-analysis catalog that the v2.1 hg38 catalog's annotation step used. It's both a
    # source of locus definitions and the catalog that the KnownDiseaseAssociatedLocus annotations come from.
    NS(name="KnownDiseaseAssociatedLociOnHg38", url="https://raw.githubusercontent.com/broadinstitute/str-analysis/063508c1decae2689ddce02413eb85e2ca6520f6/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json"),
    # STRchive release v2.28.0 (commit e964fa88)
    NS(name="STRchiveKnownDiseaseAssociatedLoci", url="https://raw.githubusercontent.com/dashnowlab/STRchive/e964fa88df3b2189ee478cbd6eaad15bb28e3a0a/data/catalogs/STRchive-disease-loci.T2T-chm13.expansionhunter.json"),
    # STRchive's locus metadata at the same commit, which has each locus's disease, evidence and CHM13 coordinates
    NS(name="STRchiveLociMetadata", url="https://raw.githubusercontent.com/dashnowlab/STRchive/e964fa88df3b2189ee478cbd6eaad15bb28e3a0a/data/STRchive-loci.json"),
    # Generated by https://github.com/broadinstitute/colab-repeat-finder/blob/main/hail_batch_pipeline/run_hail_batch_pipeline.py
    # with --reference-genome chm13, after which the loci in the chrY PARs were removed since chm13v2.0.maskedY.fasta
    # has Ns there.
    NS(name="PerfectRepeatsInCHM13", url="https://storage.googleapis.com/str-truth-set/chm13/ref/other/colab-repeat-finder/chm13_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp/chm13_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp.bed.gz"),
    NS(name="TRExplorerHg38Catalog", url=args.hg38_catalog),
    # The UCSC chain file recommended by Genovese et al. 2024 (https://doi.org/10.1093/bioinformatics/btae038) rather
    # than hg38-chm13v2.over.chain.gz, which has erroneous mappings caused by a chaintools bug.
    NS(name="Hg38ToCHM13ChainFile", url="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz"),
    NS(name="GencodeCatLiftoffGTF", url="https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/genes/catLiftOffGenesV1.gtf.gz"),
    NS(name="RefseqGTF", url="https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz"),
    # The Miga lab's centromeric satellite (CenSat) annotation of CHM13, made with the same workflow as the hg38 one
    # that the hg38 pipeline uses, so that both catalogs' CentromereAndSatelliteRegion fields have the same classes.
    # It differs from the T2T consortium's chm13v2.0_censat_v2.1.bed, whose class names (hsat3, hor, bsat) don't match
    # the hg38 file's and whose ct regions also include 5 Mb of each chromosome arm next to the centromere.
    NS(name="CentromereAndSatelliteRegionsBed", url="https://raw.githubusercontent.com/hloucks/CenSatData/fb50a32b8ce443b9b0ea718e3c0c47b91507e59c/CHM13/chm13v2.0.cenSatv2.1.bed"),
]

input_paths = {}
for remote_input in remote_inputs:
    if os.path.isfile(remote_input.url):
        input_paths[remote_input.name] = os.path.abspath(remote_input.url)
    else:
        input_paths[remote_input.name] = download_unless_local_copy_is_up_to_date(remote_input.url, downloaded_data_dir)

print("Done downloading. The remaining steps don't need internet access.")

output_prefix = os.path.abspath(f"{args.output_prefix}.{MIN_MOTIF_SIZE}_to_{MAX_MOTIF_SIZE}bp_motifs")

lift_over_args = (f"--chain-file {input_paths['Hg38ToCHM13ChainFile']} "
                  f"--source-reference-fasta {args.hg38_reference_fasta} "
                  f"--target-reference-fasta {args.chm13_reference_fasta}")

# Keep only the STRchive loci whose disease association is confident: those with Definitive, Strong, Moderate or
# Limited evidence, plus those that STRchive's website tags as New, which for this 2026 STRchive version means first
# reported in 2024 or later. This filters both STRchive's locus metadata and its CHM13 catalog. These commands, like
# all of the following commands up to step 2, share step 1's number so that rerunning step 1 reapplies them.
strchive_loci_metadata_path = os.path.abspath("STRchive-loci.confident_disease_associations.json")
strchive_catalog_path = os.path.abspath("STRchive-disease-loci.T2T-chm13.confident_disease_associations.json")
for input_name, output_path in [
    ("STRchiveLociMetadata", strchive_loci_metadata_path),
    ("STRchiveKnownDiseaseAssociatedLoci", strchive_catalog_path),
]:
    run(f"""python3 -u {base_dir}/scripts/select_STRchive_loci_with_confident_disease_associations.py \
        --strchive-loci-json {input_paths['STRchiveLociMetadata']} \
        --min-year-first-reported-for-new-loci 2024 \
        -o {output_path} \
        {input_paths[input_name]}""", step_number=1)

# Lift the known disease-associated loci over to CHM13, keeping all their fields, and replacing the hg38 coordinates in
# the LocusIds of alternative definitions (like EP400__12-132062524-132062611-CAG) with their CHM13 coordinates. Then,
# for the catalog that the KnownDiseaseAssociatedLocus annotations below come from, add STRchive's CHM13 definitions of
# the loci that didn't lift over.
known_disease_associated_loci_path = os.path.abspath("known_disease_associated_loci.lifted_over_to_chm13.json")
run(f"""python3 -u {base_dir}/scripts/lift_over_catalog.py {lift_over_args} \
    --update-coordinates-in-locus-ids \
    --output-dropped-loci-tsv {known_disease_associated_loci_path.replace('.json', '.dropped_loci.tsv')} \
    -o {known_disease_associated_loci_path} \
    {input_paths['KnownDiseaseAssociatedLociOnHg38']}""", step_number=1)

known_disease_associated_loci_for_annotation_path = known_disease_associated_loci_path.replace(
    ".json", ".with_STRchive_loci.json")
run(f"""python3 -u {base_dir}/scripts/add_STRchive_CHM13_loci_to_known_disease_associated_loci_catalog.py \
    --strchive-loci-json {strchive_loci_metadata_path} \
    --hg38-known-disease-associated-loci {input_paths['KnownDiseaseAssociatedLociOnHg38']} \
    -o {known_disease_associated_loci_for_annotation_path} \
    {known_disease_associated_loci_path}""", step_number=1)

# As a source of locus definitions, and for the validation step, use the loci that have diseases. This includes the
# alternative definitions (like EP400__12-132111361-132111451-CAG). The merges keep only one of a locus's overlapping
# definitions, and step 11 then appends the others again.
known_disease_associated_loci_with_diseases_path = known_disease_associated_loci_path.replace(
    ".json", ".with_diseases.json")
if not args.dry_run:
    with open(known_disease_associated_loci_path) as f:
        known_disease_associated_loci = json.load(f)
    with open(known_disease_associated_loci_with_diseases_path, "wt") as f:
        json.dump([r for r in known_disease_associated_loci if r.get("Diseases")], f, indent=4)

# As for the hg38 catalog: split compound definitions, change the RFC1 motif from AARRG to AAAAG since the catalog
# doesn't currently support IUPAC codes, and compress
known_disease_associated_loci_source_path = known_disease_associated_loci_with_diseases_path.replace(
    ".json", ".split.json")
run(f"python3 -u -m str_analysis.split_adjacent_loci_in_expansion_hunter_catalog "
    f"-o {known_disease_associated_loci_source_path} {known_disease_associated_loci_with_diseases_path}", step_number=2)
run(f"sed -i 's/AARRG/AAAAG/g' {known_disease_associated_loci_source_path}", step_number=3)
run(f"gzip -f {known_disease_associated_loci_source_path}", step_number=4)
known_disease_associated_loci_source_path += ".gz"

# STRchive's catalog defines some loci as several adjacent repeats without the VariantIds that
# str_analysis.split_adjacent_loci_in_expansion_hunter_catalog requires, so split them by converting to BED instead
strchive_source_path = os.path.abspath("STRchive-disease-loci.T2T-chm13.split_adjacent_repeats.bed.gz")
run(f"python3 -u -m str_analysis.convert_expansion_hunter_catalog_to_bed --split-adjacent-repeats "
    f"-o {strchive_source_path} {strchive_catalog_path}", step_number=5)

lifted_hg38_catalog_path = os.path.abspath(f"TRExplorer.repeat_catalog_v{VERSION}.hg38.lifted_over_to_chm13.json.gz")
run(f"""python3 -u {base_dir}/scripts/lift_over_catalog.py {lift_over_args} \
    --keep-only-core-fields \
    --original-locus-id-field Hg38LocusId \
    --output-dropped-loci-tsv {lifted_hg38_catalog_path.replace('.json.gz', '.dropped_loci.tsv')} \
    -o {lifted_hg38_catalog_path} \
    {input_paths['TRExplorerHg38Catalog']}""", step_number=6)

# List of source catalogs, in order. If a locus is defined in more than one catalog (ie. overlapping boundaries,
# same motif), then the definition in the catalog that's earlier in the list takes precedence. The first merge combines
# the sources that the hg38 catalog's first merge used (known disease-associated loci and perfect repeats), and the
# second merge adds the lifted over hg38 catalog.
source_catalogs_in_order = [
    NS(name="TRExplorerCHM13:KnownDiseaseAssociatedLoci",         merge="v1", path=known_disease_associated_loci_source_path),
    NS(name="TRExplorerCHM13:STRchiveKnownDiseaseAssociatedLoci", merge="v1", path=strchive_source_path),
    NS(name="TRExplorerCHM13:PerfectRepeatsInReference",          merge="v1", path=input_paths["PerfectRepeatsInCHM13"]),
    NS(name="TRExplorerCHM13:LiftedOverFromTRExplorerHg38",      merge="v2", path=lifted_hg38_catalog_path),
]

filtered_source_catalog_paths = {}
for source_catalog in source_catalogs_in_order:
    if source_catalog.path.endswith(".json.gz"):
        filtered_catalog_path = source_catalog.path.replace(".json.gz", ".filtered.json.gz")
    elif source_catalog.path.endswith(".bed.gz"):
        filtered_catalog_path = source_catalog.path.replace(".bed.gz", ".filtered.json.gz")
    else:
        raise ValueError(f"Unexpected file extension for {source_catalog.path}")

    filtered_catalog_path = os.path.abspath(os.path.basename(filtered_catalog_path))
    filtered_source_catalog_paths[source_catalog.name] = filtered_catalog_path

    # As the hg38 pipeline does for its TRExplorerV2 sources, drop STRchive loci whose motif has IUPAC codes (like its
    # (GCN)* definitions of HOXD13, ZIC3 and PABPN1), since the merge treats GCN as a different motif from the GCG or
    # GCC of the same locus's definition in the first source, and would keep both. Also drop STRchive regions shorter
    # than one motif, which splitting its compound loci turns into loci of their own (like chr17:17755051-17755053
    # (TTTCA)*, a placeholder for the pathogenic motif next to RAI1's TTTTA repeat).
    extra_args = ""
    if source_catalog.name == "TRExplorerCHM13:STRchiveKnownDiseaseAssociatedLoci":
        extra_args = "--discard-loci-with-non-ACGT-bases-in-motif --min-repeats-in-reference 1"

    run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog \
        --verbose \
        --reference-fasta {args.chm13_reference_fasta} \
        {extra_args} \
        --min-motif-size {MIN_MOTIF_SIZE} \
        --max-motif-size {MAX_MOTIF_SIZE} \
        --min-interval-size-bp 1 \
        --skip-gene-annotations \
        --skip-mappability-annotations \
        --skip-disease-loci-annotations \
        --set-locus-id \
        --discard-loci-with-non-ACGT-bases-in-reference \
        --discard-loci-with-non-ACGTN-bases-in-motif \
        --output-path {filtered_catalog_path} \
        {source_catalog.path}""", step_number=7)

    print(f"Stats for {source_catalog.path}")
    run(f"python3 -m str_analysis.compute_catalog_stats --reference-fasta {args.chm13_reference_fasta} --verbose "
        f"{filtered_catalog_path}", step_number=8)

# NOTE: like the hg38 catalog, this doesn't use the --merge-adjacent-loci-with-same-motif option for
# str_analysis.merge_loci, so that locus definitions are preserved as they appear in the individual source catalogs.
source_catalog_paths_for_merge_command_v1 = " ".join([
    f"{source_catalog.name}:{filtered_source_catalog_paths[source_catalog.name]}"
    for source_catalog in source_catalogs_in_order if source_catalog.merge == "v1"
])

run(f"""python3 -u -m str_analysis.merge_loci --verbose \
    --add-found-in-fields \
    --output-format JSON \
    --discard-extra-fields-from-input-catalogs \
    --overlap-fraction 0.66 \
    --overlapping-loci-action keep-first \
    --write-merge-stats-tsv \
    --output-prefix {output_prefix}.merged_v1 \
    {source_catalog_paths_for_merge_command_v1}""", step_number=9)

run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog --verbose \
    --reference-fasta {args.chm13_reference_fasta} \
    --known-disease-associated-loci {known_disease_associated_loci_for_annotation_path} \
    --skip-gene-annotations \
    --skip-mappability-annotations \
    --min-motif-size {MIN_MOTIF_SIZE} \
    --max-motif-size {MAX_MOTIF_SIZE} \
    --min-interval-size-bp 1 \
    --discard-overlapping-intervals-with-similar-motifs \
    --output-path {output_prefix}.merged_v1.annotated.json.gz \
    {output_prefix}.merged_v1.json.gz""", step_number=10)

source_catalog_paths_for_merge_command_v2 = " ".join([
    f"merged:{output_prefix}.merged_v1.annotated.json.gz",
] + [
    f"{source_catalog.name}:{filtered_source_catalog_paths[source_catalog.name]}"
    for source_catalog in source_catalogs_in_order if source_catalog.merge == "v2"
])

# --only-compare-loci-from-different-catalogs keeps overlapping loci of the lifted hg38 catalog that the hg38 catalog
# deliberately has both of, like a locus and the wider definition that grew out of it (see the comment above step 15
# of the hg38 pipeline), while still dropping lifted loci that are redundant with the first merge's loci.
run(f"""python3 -u -m str_analysis.merge_loci --verbose \
    --add-found-in-fields \
    --output-format JSON \
    --overlapping-loci-action keep-first \
    --write-merge-stats-tsv \
    --overlap-fraction 0.66 \
    --motif-length-match-sufficient-for-VNTRs \
    --min-jaccard-similarity 0.66 \
    --only-compare-loci-from-different-catalogs \
    --output-prefix {output_prefix}.final_merged \
    {source_catalog_paths_for_merge_command_v2}""", step_number=11)

# Append the known disease-associated loci that the merges dropped, so that none of them are deduplicated. Step 9's
# merge keeps only the first of two overlapping definitions with the same motif, even from the same source, so it
# drops one of each pair of a primary and an alternative definition. The lifted disease loci are sorted by position, so
# it keeps whichever starts first, eg. EP400__12-132111361-132111451-CAG rather than EP400 at chr12:132111385-132111451.
# --discard-extra-fields drops the known disease-associated loci catalog's other fields (like Diseases), as step 9's
# merge does. This and the next command share step 11's number so that rerunning step 11 reapplies them.
run(f"""python3 -u {base_dir}/scripts/append_loci_not_already_in_catalog.py \
    --catalog-of-loci-to-append {filtered_source_catalog_paths['TRExplorerCHM13:KnownDiseaseAssociatedLoci']} \
    --source-name TRExplorerCHM13:KnownDiseaseAssociatedLoci \
    --discard-extra-fields \
    --output-catalog-json-path {output_prefix}.final_merged.with_all_known_disease_associated_loci.json.gz \
    {output_prefix}.final_merged.json.gz""", step_number=11)

# Remove the loci that are too close to a chromosome end for genotyping tools, like the telomeric repeats at CHM13's
# chromosome ends (hg38's chromosome ends are Ns).
final_merged_catalog_path = f"{output_prefix}.final_merged.without_loci_near_chromosome_ends.json.gz"
run(f"""python3 -u {base_dir}/scripts/remove_loci_near_chromosome_ends.py \
    -R {args.chm13_reference_fasta} \
    -o {final_merged_catalog_path} \
    {output_prefix}.final_merged.with_all_known_disease_associated_loci.json.gz""", step_number=11)

annotated_catalog_path = f"{output_prefix}.EH.with_annotations.json.gz"
run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog --verbose \
    --reference-fasta {args.chm13_reference_fasta} \
    --gene-models-source gencode \
    --gencode-gtf {input_paths['GencodeCatLiftoffGTF']} \
    --gene-models-source refseq \
    --refseq-gtf {input_paths['RefseqGTF']} \
    --skip-mappability-annotations \
    --known-disease-associated-loci {known_disease_associated_loci_for_annotation_path} \
    --min-motif-size {MIN_MOTIF_SIZE} \
    --max-motif-size {MAX_MOTIF_SIZE} \
    --min-interval-size-bp 1 \
    --output-path {annotated_catalog_path} \
    {final_merged_catalog_path}""", step_number=12)

# Give each locus the Hg38LocusId of its hg38 counterpart. When a lifted hg38 locus matches a known disease-associated
# locus or a perfect repeat, step 11's merge keeps the earlier locus and discards the lifted one along with its
# Hg38LocusId, so this finds the counterpart again among the lifted loci, using step 11's thresholds. Known
# disease-associated loci that still have no match (like NOTCH2NLC and XYLT1, whose hg38 loci don't lift over) get the
# hg38 LocusId that the known disease-associated loci catalog's TRExplorerV2 field gives for their label. It also adds
# an Hg38LocusMatchType field (SameBoundaries, Hg38LocusIsWider, Hg38LocusIsNarrower, Hg38LocusIsShifted or
# SameKnownDiseaseAssociatedLocus). These commands share step 12's number so that rerunning step 12 reapplies them.
run(f"""python3 -u {base_dir}/scripts/add_hg38_locus_ids_using_lifted_over_hg38_catalog.py \
    --lifted-over-hg38-catalog {filtered_source_catalog_paths['TRExplorerCHM13:LiftedOverFromTRExplorerHg38']} \
    --min-overlap-fraction 0.66 \
    --min-jaccard-similarity 0.66 \
    --motif-length-match-sufficient-for-VNTRs \
    --known-disease-associated-loci {known_disease_associated_loci_for_annotation_path} \
    --output-catalog-json-path {annotated_catalog_path}.with_hg38_locus_ids.json.gz \
    {annotated_catalog_path}""", step_number=12)
run(f"mv {annotated_catalog_path}.with_hg38_locus_ids.json.gz {annotated_catalog_path}", step_number=12)

# As in the hg38 pipeline, label the loci that the exact-boundary match above missed using the known disease-associated
# loci catalog's TRExplorerV2 field, which holds hg38 LocusIds. The step above put each locus's hg38 LocusId in the
# Hg38LocusId field. Since step 11 appends every known disease-associated locus's exact definition, this is a safety
# net that labels nothing with the current inputs. These commands share step 12's number so that rerunning step 12
# reapplies them.
run(f"""python3 {base_dir}/scripts/add_known_disease_associated_locus_labels_using_TRExplorerV2_ids.py \
    --known-disease-associated-loci {known_disease_associated_loci_for_annotation_path} \
    --catalog-locus-id-field Hg38LocusId \
    --output-catalog-json-path {annotated_catalog_path}.with_TRExplorerV2_disease_locus_labels.json.gz \
    {annotated_catalog_path}""", step_number=12)
run(f"mv {annotated_catalog_path}.with_TRExplorerV2_disease_locus_labels.json.gz {annotated_catalog_path}", step_number=12)

# Label each locus with the class of the centromeric satellite region it overlaps, if any (eg. HSat3, active_hor, or ct
# for a centromeric transition region). These commands also share step 12's number.
run(f"""python3 {base_dir}/scripts/add_centromere_and_satellite_region_annotations.py \
    --centromere-and-satellite-regions-bed {input_paths['CentromereAndSatelliteRegionsBed']} \
    --output-catalog-json-path {annotated_catalog_path}.with_centromere_and_satellite_regions.json.gz \
    {annotated_catalog_path}""", step_number=12)
run(f"mv {annotated_catalog_path}.with_centromere_and_satellite_regions.json.gz {annotated_catalog_path}", step_number=12)

# create a version of the ExpansionHunter catalog without extra annotations
run(f"""python3 {base_dir}/scripts/write_expansion_hunter_catalog_without_annotations.py \
    --output-catalog-json-path {output_prefix}.EH.json.gz \
    {annotated_catalog_path}""", step_number=13)

run(f"python3 -m str_analysis.filter_out_loci_with_Ns_in_flanks "
    f"-R {args.chm13_reference_fasta} "
    f"-o {output_prefix}.EH.without_loci_with_Ns_in_flanks.json.gz "
    f"--output-list-of-filtered-loci {output_prefix}.loci_with_Ns_in_flanks.txt "
    f"{output_prefix}.EH.json.gz", step_number=14)
run(f"mv {output_prefix}.EH.without_loci_with_Ns_in_flanks.json.gz {output_prefix}.EH.json.gz", step_number=15)

release_files = [
    f"{output_prefix}.bed.gz",
    f"{output_prefix}.bed.gz.tbi",
    f"{annotated_catalog_path}",
    f"{output_prefix}.EH.json.gz",
    f"{output_prefix}.TRGT.bed",
    f"{output_prefix}.LongTR.bed",
    f"{output_prefix}.HipSTR.bed",
    f"{output_prefix}.GangSTR.bed",
    f"{output_prefix}.ATaRVa.bed.gz",
    f"{output_prefix}.ATaRVa.bed.gz.tbi",
]

# annotate overlapping loci and output merged tandem repeat regions
run(f"""python3 {base_dir}/scripts/annotate_overlapping_loci_and_output_merged_tandem_repeat_regions.py \
    --verbose \
    -o {annotated_catalog_path}.with_overlap_annotations.json.gz \
    {annotated_catalog_path}""", step_number=16)
run(f"mv {annotated_catalog_path}.with_overlap_annotations.json.gz {annotated_catalog_path}", step_number=17)

# convert to BED, and annotate with "TRsInRegion" based on adjacent loci
run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_bed --split-adjacent-repeats "
    f"{annotated_catalog_path}  --output-file {output_prefix}.bed.gz", step_number=18)

run(f"python3 -m str_analysis.add_adjacent_loci_to_expansion_hunter_catalog "
    f"--ref-fasta {args.chm13_reference_fasta} "
    f"--source-of-adjacent-loci {output_prefix}.bed.gz "
    f"--add-extra-field TRsInRegion "
    f"--only-add-extra-fields "
    f"-o {annotated_catalog_path}.with_adjacent_loci_annotation.json.gz "
    f"{annotated_catalog_path}", step_number=19)
run(f"mv {annotated_catalog_path}.with_adjacent_loci_annotation.json.gz {annotated_catalog_path}", step_number=20)

# convert to TSV, leaving out the core columns of annotations that this catalog doesn't have
run(f"""python3 {base_dir}/scripts/convert_catalog_json_to_tsv.py \
    --omit-core-columns-missing-from-catalog \
    --output-tsv-path {annotated_catalog_path.replace('.json.gz', '') + '.tsv.gz'} \
    {annotated_catalog_path}""", step_number=21)

# convert the catalog from ExpansionHunter catalog format to TRGT, LongTR, HipSTR, and GangSTR formats
run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_trgt_catalog -R {args.chm13_reference_fasta} --split-adjacent-repeats {annotated_catalog_path}  --output-file {output_prefix}.TRGT.bed", step_number=22)
run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_longtr_format  {annotated_catalog_path}  --output-file {output_prefix}.LongTR.bed", step_number=23)
run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_hipstr_format  {annotated_catalog_path}  --output-file {output_prefix}.HipSTR.bed", step_number=24)
run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_gangstr_spec   {annotated_catalog_path}  --output-file {output_prefix}.GangSTR.bed", step_number=25)

# ATaRVa (https://github.com/SowpatiLab/ATaRVa) reads a 5-column BED: chromosome, start, end,
# motif, motif length. It reads one motif per row, hence --split-adjacent-repeats, and requires
# the file to be bgzipped and tabix-indexed, which the converter does for a .gz output path.
run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_bed --split-adjacent-repeats "
    f"--motif-size-column {annotated_catalog_path}  --output-file {output_prefix}.ATaRVa.bed.gz",
    step_number=26)

# Confirm that the TRGT catalog passes 'trgt validate'. trgt exits 0 even when it prints "Validation failed",
# so stop the pipeline unless its output contains "Validation successful".
run(f"trgt validate --genome {args.chm13_reference_fasta}  --repeats {output_prefix}.TRGT.bed 2>&1 "
    f"| tee {output_prefix}.trgt_validate.log  && grep -q 'Validation successful' {output_prefix}.trgt_validate.log",
    step_number=27)

# Print missing values
run(f"python3 {base_dir}/scripts/print_missing_values_percentages_in_json_or_tsv.py  {annotated_catalog_path}", step_number=28)

# Perform basic internal consistency checks on the JSON catalog. --check-for-presence-of-annotations isn't used since
# it requires the mappability and MANE annotations that this catalog doesn't have yet.
run(f"python3 {base_dir}/scripts/validate_catalog.py "
    f"--known-pathogenic-loci-json-path {known_disease_associated_loci_source_path} "
    f"--check-for-presence-of-all-known-loci "
    f"{annotated_catalog_path}", step_number=29)

# copy files to the release_draft folder and compute catalog stats
updated_release_files = []
for path in release_files:
    if path.endswith(".bed"):
        updated_release_files.append(f"{path}.gz")
    else:
        if path.endswith(".json") or path.endswith(".json.gz") and ".EH." in path:
            run(f"python3 {base_dir}/scripts/validate_json.py -k LocusId -k LocusStructure -k ReferenceRegion -k VariantType {path}", step_number=30)
        updated_release_files.append(path)

# Compress and tabix-index the release BEDs so users can query them by region. Compression and
# indexing are one step on purpose: as separate steps, resuming with --only-step could index an
# archive without first rebuilding it from the .bed the conversion steps just wrote, and publish
# stale contents. The plain BED and the ATaRVa BED arrive already compressed and indexed.
for path in list(updated_release_files):
    if not path.endswith(".bed.gz") or f"{path}.tbi" in updated_release_files:
        continue
    uncompressed = path[:-len(".gz")]
    if os.path.isfile(uncompressed):
        # Always rebuild from the .bed an earlier step wrote, so the published archive can
        # never be a leftover from a previous run into this output directory.
        prepare = f"bgzip -f {uncompressed} && "
    elif os.path.isfile(path) and not is_bgzipped(path):
        # Only an older run's archive survives, compressed with gzip rather than bgzip, which
        # tabix cannot index. "bgzip -d" reads plain gzip too, and unlike "gzcat | bgzip" it
        # fails loudly on a truncated input instead of quietly writing a short file.
        prepare = f"bgzip -f -d {path} && bgzip -f {uncompressed} && "
    else:
        prepare = ""
    run(f"{prepare}tabix -f {tabix_coordinate_flags(path)} {path}", step_number=31)
    updated_release_files.append(f"{path}.tbi")

for path in updated_release_files:
    run(f"cp {path} {release_draft_folder}", step_number=32)

run(f"python3 -m str_analysis.compute_catalog_stats --reference-fasta {args.chm13_reference_fasta} --verbose {annotated_catalog_path}", step_number=33)

# Print source statistics table
run(f"python3 {base_dir}/scripts/generate_catalog_sources_stats_table.py {release_draft_folder}/{os.path.basename(annotated_catalog_path)}", step_number=34)

# Print final timing summary
end_datetime = datetime.datetime.now()
total_elapsed_time = time.time() - script_start_time
print("\n" + "="*80)
print(f"Started at:  {start_datetime.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Finished at: {end_datetime.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Total time:  {total_elapsed_time//3600:.0f}h, {(total_elapsed_time%3600)//60:.0f}m, {total_elapsed_time%60:.0f}s")
print("="*80)
