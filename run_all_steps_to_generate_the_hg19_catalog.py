"""Generate the hg19 (GRCh37) version of the TRExplorer repeat catalog.

As for the hg38 and CHM13 catalogs (https://github.com/broadinstitute/trexplorer-catalog/issues/4), each source only adds
loci that aren't redundant with the loci of the sources before it:
1. known disease-associated loci: the GRCh37 version of the str-analysis catalog that gnomAD uses, keeping the loci that
   have diseases, as the CHM13 pipeline does
2. the hg19 version of Illumina's catalog of 174k polymorphic TRs (https://github.com/Illumina/RepeatCatalogs), with its
   compound loci split into their repeats, like the hg38 catalog's TRExplorerV1:Illumina174kPolymorphicTRs source
3. perfect repeats in hg19 that span at least 9bp and 3 or more repeats, found by colab-repeat-finder's
   modal_pipeline/find_perfect_repeats_on_modal.py with the same settings as the hg38 and CHM13 perfect repeats
4. the TRExplorer v2.1 hg38 catalog, lifted over to hg19 with scripts/lift_over_catalog.py and UCSC's hg38ToHg19 chain

None of the known disease-associated loci are deduplicated: both a locus's primary definition and its alternative
definitions (like EP400 and EP400__12-132547069-132547156-CAG) are kept, even though they overlap.

The reference is the Broad's GRCh37 build (Homo_sapiens_assembly19, also called b37) with chr-prefixed contig names. Its
chr1-22 and chrX have the same sequence as UCSC's hg19, as does its chrY outside the pseudoautosomal regions (PARs),
which, like the hg38 and CHM13 references, it masks with Ns. Its chrM is the revised Cambridge Reference Sequence, like hg38's, rather than
UCSC hg19's older chrM, which is what the hg38ToHg19 chain maps hg38's chrM to. So the catalog leaves out chrM. None of
the sources has chrM loci by default: the perfect repeats were found without chrM, and the hg38 catalog that's lifted
over by default (the EH.json.gz release) has none, since the hg38 pipeline leaves them out of it. Step 6 still
excludes chrM, in case --hg38-catalog is an hg38 catalog that has chrM loci, like the annotated one, whose 14 chrM loci
would otherwise land at positions that don't match this chrM.

The annotations are:
- Gencode genes from GENCODE's GRCh37 mapping of release 50 (gencode.v50lift37), the release that the hg38 catalog
  uses, so that the GencodeGeneId and GencodeTranscriptId fields hold the same Ensembl ids as in the hg38 catalog
- RefSeq genes from UCSC's hg19.ncbiRefSeq.gtf.gz, which UCSC last updated on 11/02/2021, so it has an older set of
  RefSeq transcripts than the hg38 catalog's
- MANE genes from the gencode.v50lift37 transcripts tagged MANE_Select or MANE_Plus_Clinical, since MANE only
  publishes GTF files for GRCh38
- mappability from UCSC's ENCODE/CRG 36-mer alignability track for hg19, which allows up to 2 mismatches. It's the
  closest public counterpart of the hg38 catalog's k36_m2 track, though it wasn't necessarily made the same way.
- KnownDiseaseAssociatedLocus from the GRCh37 known disease-associated loci catalog
- Hg38LocusId and Hg38LocusMatchType, from the lifted hg38 catalog, as for the CHM13 catalog. For a known
  disease-associated locus that matches no lifted hg38 locus, the hg38 LocusId comes from the known disease-associated
  loci catalog's TRExplorerV2 field. The GRCh37 version of that catalog has no TRExplorerV2 field, so step 12 copies it
  from the GRCh38 version by scripts/add_TRExplorerV2_ids_to_GRCh37_known_disease_associated_loci.py.
- CentromereAndSatelliteRegion from the Miga lab's hg38 CenSat annotation, lifted to hg19 by
  scripts/lift_over_bed_regions.py since there's no hg19 CenSat annotation. hg19's centromeres are gaps, so only the
  satellite arrays and centromeric transition regions outside them get labels. The script keeps only the lifted pieces
  that map back into the hg38 CenSat regions with the hg19ToHg38 chain, which drops pieces that come from alignments to
  paralogous copies (like GRCh38's falsely duplicated copies of 21q22, which would otherwise label 21q22 loci 'ct').
- PseudoautosomalRegion, and UCSC_SegDupMaxSequenceIdentity from UCSC's hg19 genomicSuperDups table. There's no hg19
  version of the Vollger et al. 2022 segdups.
- the hg38 catalog's variation clusters, lifted over by scripts/lift_over_variation_clusters_tsv.py as for the CHM13
  catalog. Loci with no hg38 counterpart aren't in a variation cluster.

The hg38-only annotations (allele frequencies, TenK10K, HPRC256 and AoU) aren't added, and like the CHM13 catalog,
there's no hg19 version of the superseded TRGT catalog that only has the repeats outside variation clusters.
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

import str_analysis.annotate_and_filter_str_catalog
from str_analysis.utils.file_utils import tee_stdout_and_stderr_to_log_file

# UCSC's hg19 mappability track has no values over assembly gaps, and a str-analysis from before
# compute_mean_mappability crashes on the loci next to a gap rather than giving their flanks a mappability of 0. The
# install block above doesn't upgrade an existing install, so check that the installed str-analysis has it.
if not hasattr(str_analysis.annotate_and_filter_str_catalog, "compute_mean_mappability"):
    raise RuntimeError("The installed str-analysis package can't annotate the loci next to the assembly gaps of the "
                       "hg19 mappability track. Upgrade it to a version whose annotate_and_filter_str_catalog module "
                       "has compute_mean_mappability.")

# imported after the install block above, since str-analysis is what installs the requests package this module uses
import run_all_steps_helpers
from run_all_steps_helpers import chdir

# Step 8 merges the known disease-associated loci, the Illumina catalog and the perfect repeats with keep-first, which
# needs a merge_loci that doesn't let the loci it discarded make later loci redundant
run_all_steps_helpers.check_that_merge_loci_ignores_discarded_loci()


VERSION = "2.1"
MIN_MOTIF_SIZE = 1
MAX_MOTIF_SIZE = 1000

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--hg19-reference-fasta", default="Homo_sapiens_assembly19.chr_prefixed.fasta",
                    help="Local path of the hg19 reference genome FASTA file, which must have a .fai index. Use the "
                         "Broad's GRCh37 build (Homo_sapiens_assembly19.fasta) with chr-prefixed contig names, which "
                         "masks the chrY PARs with Ns, as the hg19 perfect repeats source was made from it")
parser.add_argument("--hg38-reference-fasta", default="hg38.fa",
                    help="Local path of the hg38 reference genome FASTA file, which must have a .fai index. The "
                         "liftover step uses it to compare each locus's purity before and after liftover")
parser.add_argument("--hg38-catalog",
                    default=f"https://storage.googleapis.com/tandem-repeat-catalog/v{VERSION}/"
                            f"TRExplorer.repeat_catalog_v{VERSION}.hg38.1_to_1000bp_motifs.EH.json.gz",
                    help="Local path, gs:// path or URL of the TRExplorer hg38 catalog to lift over to hg19")
parser.add_argument("--perfect-repeats-in-hg19",
                    default="https://storage.googleapis.com/str-truth-set/hg19/ref/other/colab-repeat-finder/"
                            "hg19_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp/"
                            "hg19_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp.bed.gz",
                    help="Local path, gs:// path or URL of the BED file of perfect repeats in hg19 made by "
                         "colab-repeat-finder's modal_pipeline/find_perfect_repeats_on_modal.py. The default is where "
                         "that BED is to be published, next to the hg38 and CHM13 perfect repeats, and it works only "
                         "once the BED has been uploaded there. Until then, pass the path of the local BED.")
parser.add_argument("--output-prefix", default=f"TRExplorer.repeat_catalog_v{VERSION}.hg19")
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
# the server changes. This is a subfolder of the hg38 pipeline's downloaded_data folder so that the pipelines never
# overwrite each other's copy of a file that has the same name in both, like genomicSuperDups.txt.gz.
downloaded_data_dir = os.path.join(base_dir, "downloaded_data", "hg19")
run(f"mkdir -p {downloaded_data_dir}")

for key in "hg19_reference_fasta", "hg38_reference_fasta":
    if not os.path.isfile(getattr(args, key)):
        parser.error(f"{key} must be a local file, and {getattr(args, key)} isn't one")
    setattr(args, key, os.path.abspath(getattr(args, key)))

# local paths are relative to the launch directory, so resolve them before changing to the results folder
for key in "hg38_catalog", "perfect_repeats_in_hg19":
    if os.path.isfile(getattr(args, key)):
        setattr(args, key, os.path.abspath(getattr(args, key)))

# The working directory name includes the genome so that builds of different genomes started on the same day don't
# write to the same folder.
working_dir = os.path.abspath(f"results__hg19__{args.timestamp}")

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
    # The GRCh37 version of the str-analysis known disease-associated loci catalog, at the commit that made it
    # concordant with the GRCh38 version. It's both a source of locus definitions and the catalog that the
    # KnownDiseaseAssociatedLocus annotations come from.
    NS(name="KnownDiseaseAssociatedLociOnHg19", url="https://raw.githubusercontent.com/broadinstitute/str-analysis/7476acc51eaa0e502674901b536cea42dca5650c/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh37.json"),
    # The same version of the GRCh38 catalog that the CHM13 pipeline uses, which is the one the v2.1 hg38 catalog's
    # annotation step used. It's only read for its TRExplorerV2 field, which holds the hg38 catalog's LocusId of each
    # known disease-associated locus, and which step 12 copies to the GRCh37 catalog.
    NS(name="KnownDiseaseAssociatedLociOnHg38", url="https://raw.githubusercontent.com/broadinstitute/str-analysis/063508c1decae2689ddce02413eb85e2ca6520f6/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json"),
    NS(name="Illumina174kPolymorphicTRsOnHg19", url="https://raw.githubusercontent.com/Illumina/RepeatCatalogs/dd0dc170170ad442d232c6e3426e41446405d659/hg19/variant_catalog.json"),
    NS(name="PerfectRepeatsInHg19", url=args.perfect_repeats_in_hg19),
    NS(name="TRExplorerHg38Catalog", url=args.hg38_catalog),
    # The variation clusters of the hg38 catalog (the hg38 pipeline's --variation-clusters-tsv), which are lifted over
    NS(name="VariationClustersTsvOnHg38", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/genome_clusters-v2.tsv.gz"),
    NS(name="Hg38ToHg19ChainFile", url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"),
    # The chain file back from hg19 to hg38, which lift_over_bed_regions.py uses to drop the pieces of the lifted CenSat
    # regions that come from alignments to paralogous copies
    NS(name="Hg19ToHg38ChainFile", url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"),
    NS(name="GencodeGTF", url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh37_mapping/gencode.v50lift37.annotation.gtf.gz"),
    NS(name="RefseqGTF", url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.ncbiRefSeq.gtf.gz"),
    NS(name="MappabilityTrack", url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig"),
    # The Miga lab's centromeric satellite (CenSat) annotation of hg38, which the hg38 pipeline uses, lifted over below
    NS(name="CentromereAndSatelliteRegionsBedOnHg38", url="https://raw.githubusercontent.com/hloucks/CenSatData/fb50a32b8ce443b9b0ea718e3c0c47b91507e59c/HG38/hg38.cenSatv1.0.bed"),
    # UCSC's segmental duplications table, for the UCSC_SegDupMaxSequenceIdentity field
    NS(name="UcscSegDups", url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz"),
]

input_paths = {}
for remote_input in remote_inputs:
    if os.path.isfile(remote_input.url):
        input_paths[remote_input.name] = os.path.abspath(remote_input.url)
    else:
        input_paths[remote_input.name] = download_unless_local_copy_is_up_to_date(remote_input.url, downloaded_data_dir)

print("Done downloading. The remaining steps don't need internet access.")

output_prefix = os.path.abspath(f"{args.output_prefix}.{MIN_MOTIF_SIZE}_to_{MAX_MOTIF_SIZE}bp_motifs")

# As a source of locus definitions, and for the validation step, use the loci that have diseases, as the CHM13 pipeline
# does. This includes the alternative definitions (like EP400__12-132547069-132547156-CAG). The merges keep only one of a
# locus's overlapping definitions, and step 11 then appends the others again.
known_disease_associated_loci_with_diseases_path = os.path.abspath(
    "variant_catalog_without_offtargets.GRCh37.with_diseases.json")
if not args.dry_run:
    with open(input_paths["KnownDiseaseAssociatedLociOnHg19"]) as f:
        known_disease_associated_loci = json.load(f)
    with open(known_disease_associated_loci_with_diseases_path, "wt") as f:
        json.dump([r for r in known_disease_associated_loci if r.get("Diseases")], f, indent=4)

# As for the hg38 catalog: split compound definitions, change the RFC1 motif from AARRG to AAAAG since the catalog
# doesn't currently support IUPAC codes, and compress
known_disease_associated_loci_source_path = known_disease_associated_loci_with_diseases_path.replace(
    ".json", ".split.json")
run(f"python3 -u -m str_analysis.split_adjacent_loci_in_expansion_hunter_catalog "
    f"-o {known_disease_associated_loci_source_path} {known_disease_associated_loci_with_diseases_path}", step_number=1)
run(f"sed -i 's/AARRG/AAAAG/g' {known_disease_associated_loci_source_path}", step_number=2)
run(f"gzip -f {known_disease_associated_loci_source_path}", step_number=3)
known_disease_associated_loci_source_path += ".gz"

# Split the Illumina catalog's compound loci (like the HTT locus's CAG and CCG repeats) into their repeats by converting
# it to BED, as for the hg38 catalog's Illumina source
illumina_source_path = os.path.abspath("Illumina174kPolymorphicTRs.hg19.split_adjacent_repeats.bed.gz")
run(f"python3 -u -m str_analysis.convert_expansion_hunter_catalog_to_bed --split-adjacent-repeats "
    f"-o {illumina_source_path} {input_paths['Illumina174kPolymorphicTRsOnHg19']}", step_number=4)

lift_over_args = (f"--chain-file {input_paths['Hg38ToHg19ChainFile']} "
                  f"--source-reference-fasta {args.hg38_reference_fasta} "
                  f"--target-reference-fasta {args.hg19_reference_fasta}")

lifted_hg38_catalog_path = os.path.abspath(f"TRExplorer.repeat_catalog_v{VERSION}.hg38.lifted_over_to_hg19.json.gz")
run(f"""python3 -u {base_dir}/scripts/lift_over_catalog.py {lift_over_args} \
    --keep-only-core-fields \
    --original-locus-id-field Hg38LocusId \
    --output-dropped-loci-tsv {lifted_hg38_catalog_path.replace('.json.gz', '.dropped_loci.tsv')} \
    -o {lifted_hg38_catalog_path} \
    {input_paths['TRExplorerHg38Catalog']}""", step_number=5)

# List of source catalogs, in order. If a locus is defined in more than one catalog (ie. overlapping boundaries,
# same motif), then the definition in the catalog that's earlier in the list takes precedence. The first merge combines
# the sources that the hg38 catalog's first merge used (known disease-associated loci, the Illumina catalog and perfect
# repeats), and the second merge adds the lifted over hg38 catalog, as in the CHM13 pipeline.
source_catalogs_in_order = [
    NS(name="TRExplorerHg19:KnownDiseaseAssociatedLoci",      merge="v1", path=known_disease_associated_loci_source_path),
    NS(name="TRExplorerHg19:Illumina174kPolymorphicTRs",      merge="v1", path=illumina_source_path),
    NS(name="TRExplorerHg19:PerfectRepeatsInReference",       merge="v1", path=input_paths["PerfectRepeatsInHg19"]),
    NS(name="TRExplorerHg19:LiftedOverFromTRExplorerHg38",    merge="v2", path=lifted_hg38_catalog_path),
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

    # --exclude-chrom chrM leaves out the chrM loci of a lifted hg38 catalog that has them (see the docstring)
    run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog \
        --verbose \
        --reference-fasta {args.hg19_reference_fasta} \
        --min-motif-size {MIN_MOTIF_SIZE} \
        --max-motif-size {MAX_MOTIF_SIZE} \
        --min-interval-size-bp 1 \
        --skip-gene-annotations \
        --skip-mappability-annotations \
        --skip-disease-loci-annotations \
        --set-locus-id \
        --exclude-chrom chrM \
        --discard-loci-with-non-ACGT-bases-in-reference \
        --discard-loci-with-non-ACGTN-bases-in-motif \
        --output-path {filtered_catalog_path} \
        {source_catalog.path}""", step_number=6)

    print(f"Stats for {source_catalog.path}")
    run(f"python3 -m str_analysis.compute_catalog_stats --reference-fasta {args.hg19_reference_fasta} --verbose "
        f"{filtered_catalog_path}", step_number=7)

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
    {source_catalog_paths_for_merge_command_v1}""", step_number=8)

run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog --verbose \
    --reference-fasta {args.hg19_reference_fasta} \
    --known-disease-associated-loci {input_paths['KnownDiseaseAssociatedLociOnHg19']} \
    --skip-gene-annotations \
    --skip-mappability-annotations \
    --min-motif-size {MIN_MOTIF_SIZE} \
    --max-motif-size {MAX_MOTIF_SIZE} \
    --min-interval-size-bp 1 \
    --discard-overlapping-intervals-with-similar-motifs \
    --output-path {output_prefix}.merged_v1.annotated.json.gz \
    {output_prefix}.merged_v1.json.gz""", step_number=9)

source_catalog_paths_for_merge_command_v2 = " ".join([
    f"merged:{output_prefix}.merged_v1.annotated.json.gz",
] + [
    f"{source_catalog.name}:{filtered_source_catalog_paths[source_catalog.name]}"
    for source_catalog in source_catalogs_in_order if source_catalog.merge == "v2"
])

# The same settings as the CHM13 pipeline's step 11. --only-compare-loci-from-different-catalogs keeps overlapping loci
# of the lifted hg38 catalog that the hg38 catalog deliberately has both of, like a locus and the wider definition that
# grew out of it (see the comment above step 15 of the hg38 pipeline), while still dropping lifted loci that are
# redundant with the first merge's loci.
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
    {source_catalog_paths_for_merge_command_v2}""", step_number=10)

# Append the known disease-associated loci that the merges dropped, so that none of them are deduplicated. Step 8's
# merge keeps only the first of two overlapping definitions with the same motif, even from the same source, so it
# drops one of each pair of a primary and an alternative definition. --discard-extra-fields drops the known
# disease-associated loci catalog's other fields (like Diseases), as step 8's merge does. This and the next command
# share step 11's number so that rerunning step 11 reapplies them.
run(f"""python3 -u {base_dir}/scripts/append_loci_not_already_in_catalog.py \
    --catalog-of-loci-to-append {filtered_source_catalog_paths['TRExplorerHg19:KnownDiseaseAssociatedLoci']} \
    --source-name TRExplorerHg19:KnownDiseaseAssociatedLoci \
    --discard-extra-fields \
    --output-catalog-json-path {output_prefix}.final_merged.with_all_known_disease_associated_loci.json.gz \
    {output_prefix}.final_merged.json.gz""", step_number=11)

# Remove the loci that are too close to a chromosome end for genotyping tools. All of hg19's chromosome ends are Ns
# except chr17's, which has sequence right up to both of its ends, with telomeric repeats at its end.
final_merged_catalog_path = f"{output_prefix}.final_merged.without_loci_near_chromosome_ends.json.gz"
run(f"""python3 -u {base_dir}/scripts/remove_loci_near_chromosome_ends.py \
    -R {args.hg19_reference_fasta} \
    -o {final_merged_catalog_path} \
    {output_prefix}.final_merged.with_all_known_disease_associated_loci.json.gz""", step_number=11)

# Make the MANE gene annotations from the GENCODE transcripts tagged as MANE transcripts, and lift the CenSat regions
# over to hg19. These commands share step 12's number so that rerunning step 12 reapplies them.
mane_gtf_path = os.path.abspath("gencode.v50lift37.MANE_transcripts.gtf.gz")
run(f"""python3 -u {base_dir}/scripts/extract_MANE_transcripts_from_GENCODE_gtf.py \
    -o {mane_gtf_path} \
    {input_paths['GencodeGTF']}""", step_number=12)

centromere_and_satellite_regions_bed_path = os.path.abspath("hg38.cenSatv1.0.lifted_over_to_hg19.bed")
run(f"""python3 -u {base_dir}/scripts/lift_over_bed_regions.py \
    --chain-file {input_paths['Hg38ToHg19ChainFile']} \
    --reverse-chain-file {input_paths['Hg19ToHg38ChainFile']} \
    -o {centromere_and_satellite_regions_bed_path} \
    {input_paths['CentromereAndSatelliteRegionsBedOnHg38']}""", step_number=12)

annotated_catalog_path = f"{output_prefix}.EH.with_annotations.json.gz"
run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog --verbose \
    --reference-fasta {args.hg19_reference_fasta} \
    --gene-models-source gencode \
    --gencode-gtf {input_paths['GencodeGTF']} \
    --gene-models-source refseq \
    --refseq-gtf {input_paths['RefseqGTF']} \
    --gene-models-source mane \
    --mane-gtf {mane_gtf_path} \
    --mappability-track-bigwig {input_paths['MappabilityTrack']} \
    --known-disease-associated-loci {input_paths['KnownDiseaseAssociatedLociOnHg19']} \
    --min-motif-size {MIN_MOTIF_SIZE} \
    --max-motif-size {MAX_MOTIF_SIZE} \
    --min-interval-size-bp 1 \
    --output-path {annotated_catalog_path} \
    {final_merged_catalog_path}""", step_number=12)

# Give the GRCh37 known disease-associated loci catalog the TRExplorerV2 field of the GRCh38 one, which holds the hg38
# LocusId of each known disease-associated locus. This shares step 12's number so that rerunning step 12 reapplies it.
known_disease_associated_loci_with_TRExplorerV2_ids_path = os.path.abspath(
    "variant_catalog_without_offtargets.GRCh37.with_TRExplorerV2_ids.json")
run(f"""python3 -u {base_dir}/scripts/add_TRExplorerV2_ids_to_GRCh37_known_disease_associated_loci.py \
    --grch38-known-disease-associated-loci {input_paths['KnownDiseaseAssociatedLociOnHg38']} \
    -o {known_disease_associated_loci_with_TRExplorerV2_ids_path} \
    {input_paths['KnownDiseaseAssociatedLociOnHg19']}""", step_number=12)

# Give each locus the Hg38LocusId of its hg38 counterpart, as the CHM13 pipeline does. When a lifted hg38 locus matches
# a known disease-associated locus, an Illumina locus or a perfect repeat, step 10's merge keeps the earlier locus and
# discards the lifted one along with its Hg38LocusId, so this finds the counterpart again among the lifted loci, using
# step 10's thresholds. A known disease-associated locus that still has no match (like PABPN1's alternative definition
# PABPN1__14-23790681-23790712-GCG, which is much wider than its lifted hg38 definition) gets the hg38 LocusId that the
# TRExplorerV2 field gives for its label. It also adds an Hg38LocusMatchType field. These commands share step 12's
# number so that rerunning step 12 reapplies them.
run(f"""python3 -u {base_dir}/scripts/add_hg38_locus_ids_using_lifted_over_hg38_catalog.py \
    --lifted-over-hg38-catalog {filtered_source_catalog_paths['TRExplorerHg19:LiftedOverFromTRExplorerHg38']} \
    --min-overlap-fraction 0.66 \
    --min-jaccard-similarity 0.66 \
    --motif-length-match-sufficient-for-VNTRs \
    --known-disease-associated-loci {known_disease_associated_loci_with_TRExplorerV2_ids_path} \
    --output-catalog-json-path {annotated_catalog_path}.with_hg38_locus_ids.json.gz \
    {annotated_catalog_path}""", step_number=12)
run(f"mv {annotated_catalog_path}.with_hg38_locus_ids.json.gz {annotated_catalog_path}", step_number=12)

# As in the hg38 and CHM13 pipelines, label the loci that the exact-boundary match of the annotation above missed, using
# the known disease-associated loci catalog's TRExplorerV2 field, which holds hg38 LocusIds. The step above put each
# locus's hg38 LocusId in the Hg38LocusId field. Since step 11 appends every known disease-associated locus's exact
# definition, this is a safety net. These commands share step 12's number so that rerunning step 12 reapplies them.
run(f"""python3 {base_dir}/scripts/add_known_disease_associated_locus_labels_using_TRExplorerV2_ids.py \
    --known-disease-associated-loci {known_disease_associated_loci_with_TRExplorerV2_ids_path} \
    --catalog-locus-id-field Hg38LocusId \
    --output-catalog-json-path {annotated_catalog_path}.with_TRExplorerV2_disease_locus_labels.json.gz \
    {annotated_catalog_path}""", step_number=12)
run(f"mv {annotated_catalog_path}.with_TRExplorerV2_disease_locus_labels.json.gz {annotated_catalog_path}", step_number=12)

# Label each locus with the class of the centromeric satellite region it overlaps, if any (eg. HSat3, active_hor, or ct
# for a centromeric transition region). These commands also share step 12's number.
run(f"""python3 {base_dir}/scripts/add_centromere_and_satellite_region_annotations.py \
    --centromere-and-satellite-regions-bed {centromere_and_satellite_regions_bed_path} \
    --output-catalog-json-path {annotated_catalog_path}.with_centromere_and_satellite_regions.json.gz \
    {annotated_catalog_path}""", step_number=12)
run(f"mv {annotated_catalog_path}.with_centromere_and_satellite_regions.json.gz {annotated_catalog_path}", step_number=12)

# Label each locus that lies inside a pseudoautosomal region with PAR1 or PAR2. These commands also share step 12's
# number.
run(f"""python3 {base_dir}/scripts/add_pseudoautosomal_region_annotations.py \
    --genome-version hg19 \
    --output-catalog-json-path {annotated_catalog_path}.with_pseudoautosomal_regions.json.gz \
    {annotated_catalog_path}""", step_number=12)
run(f"mv {annotated_catalog_path}.with_pseudoautosomal_regions.json.gz {annotated_catalog_path}", step_number=12)

# Label each locus with the highest sequence identity among UCSC's segmental duplications that overlap it, if any.
# These commands also share step 12's number.
run(f"""python3 {base_dir}/scripts/add_segdup_max_sequence_identity_annotations.py \
    --segdups-file {input_paths['UcscSegDups']} \
    --segdups-source ucsc \
    --genome-version hg19 \
    --output-catalog-json-path {annotated_catalog_path}.with_segdups.json.gz \
    {annotated_catalog_path}""", step_number=12)
run(f"mv {annotated_catalog_path}.with_segdups.json.gz {annotated_catalog_path}", step_number=12)

# Lift the hg38 catalog's variation clusters over to hg19, using the Hg38LocusId field added above, and annotate each
# locus with its variation cluster as the hg38 pipeline does. These commands also share step 12's number.
lifted_variation_clusters_tsv_path = os.path.abspath("genome_clusters-v2.lifted_over_to_hg19.tsv.gz")
run(f"""python3 -u {base_dir}/scripts/lift_over_variation_clusters_tsv.py \
    --chain-file {input_paths['Hg38ToHg19ChainFile']} \
    --lifted-catalog {annotated_catalog_path} \
    -o {lifted_variation_clusters_tsv_path} \
    {input_paths['VariationClustersTsvOnHg38']}""", step_number=12)
run(f"""python3 {base_dir}/scripts/add_variation_cluster_annotations_to_catalog.py \
    --verbose \
    --output-catalog-json-path {annotated_catalog_path}.with_variation_clusters.json.gz \
    {lifted_variation_clusters_tsv_path} \
    {annotated_catalog_path}""", step_number=12)
run(f"mv {annotated_catalog_path}.with_variation_clusters.json.gz {annotated_catalog_path}", step_number=12)

# create a version of the ExpansionHunter catalog without extra annotations
run(f"""python3 {base_dir}/scripts/write_expansion_hunter_catalog_without_annotations.py \
    --output-catalog-json-path {output_prefix}.EH.json.gz \
    {annotated_catalog_path}""", step_number=13)

run(f"python3 -m str_analysis.filter_out_loci_with_Ns_in_flanks "
    f"-R {args.hg19_reference_fasta} "
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
    f"--ref-fasta {args.hg19_reference_fasta} "
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
run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_trgt_catalog -R {args.hg19_reference_fasta} --split-adjacent-repeats {annotated_catalog_path}  --output-file {output_prefix}.TRGT.bed", step_number=22)
run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_longtr_format  {annotated_catalog_path}  --output-file {output_prefix}.LongTR.bed", step_number=23)
run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_hipstr_format  {annotated_catalog_path}  --output-file {output_prefix}.HipSTR.bed", step_number=24)
run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_gangstr_spec   {annotated_catalog_path}  --output-file {output_prefix}.GangSTR.bed", step_number=25)

# ATaRVa (https://github.com/SowpatiLab/ATaRVa) reads a 5-column BED: chromosome, start, end,
# motif, motif length. It reads one motif per row, hence --split-adjacent-repeats, and requires
# the file to be bgzipped and tabix-indexed, which the converter does for a .gz output path.
run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_bed --split-adjacent-repeats "
    f"--motif-size-column {annotated_catalog_path}  --output-file {output_prefix}.ATaRVa.bed.gz",
    step_number=26)

# As in the hg38 pipeline, write the TRGT catalog with a row for every repeat plus a row for every variation cluster, and
# its LongTR version. --name-clusters-by-trid since every repeat also has a row of its own, so a single-repeat cluster
# named after its member would collide with it. These commands share step 26's number.
variation_clusters_and_all_repeats_path = os.path.abspath(
    f"TRExplorer.variation_clusters_and_all_repeats_v{VERSION}.hg19.TRGT.bed.gz")
run(f"""python3 {base_dir}/scripts/generate_TRGT_catalog_with_variation_clusters_and_all_repeats.py \
    -o {variation_clusters_and_all_repeats_path} \
    {lifted_variation_clusters_tsv_path} \
    {annotated_catalog_path}""", step_number=26)
run(f"python3 {base_dir}/scripts/convert_trgt_catalog_to_longtr_format.py --name-clusters-by-trid "
    f"{variation_clusters_and_all_repeats_path}", step_number=26)
release_files += [
    variation_clusters_and_all_repeats_path,
    variation_clusters_and_all_repeats_path.replace(".TRGT.bed.gz", ".LongTR.bed.gz"),
]

# Confirm that the TRGT catalogs pass 'trgt validate'. trgt exits 0 even when it prints "Validation failed",
# so stop the pipeline unless its output contains "Validation successful".
for trgt_catalog_path, trgt_validate_log_path in [
    (f"{output_prefix}.TRGT.bed", f"{output_prefix}.trgt_validate.log"),
    (variation_clusters_and_all_repeats_path, f"{variation_clusters_and_all_repeats_path}.trgt_validate.log"),
]:
    run(f"trgt validate --genome {args.hg19_reference_fasta}  --repeats {trgt_catalog_path} 2>&1 "
        f"| tee {trgt_validate_log_path}  && grep -q 'Validation successful' {trgt_validate_log_path}",
        step_number=27)

# Check that the repeat rows of the TRGT catalog with variation clusters reproduce the repeat catalog's own TRGT file.
# --expect-filtered-locus replaces the script's default, an hg38 locus in NOTCH2NLA that the liftover drops (it and the
# NOTCH2NLC disease locus map onto the same hg19 repeat), with a perfect repeat in the PHRF1 coding sequence that both
# the hg38 and the lifted variation clusters tables mark DEPTH.
run(f"python3 {base_dir}/scripts/validate_variation_clusters_catalog.py "
    f"--expect-filtered-locus 11-587346-587357-TGA "
    f"{variation_clusters_and_all_repeats_path} {output_prefix}.TRGT.bed", step_number=27)

# Print missing values
run(f"python3 {base_dir}/scripts/print_missing_values_percentages_in_json_or_tsv.py  {annotated_catalog_path}", step_number=28)

# Perform basic internal consistency checks on the JSON catalog. Unlike the CHM13 catalog, this one has the mappability
# and MANE annotations that --check-for-presence-of-annotations requires.
run(f"python3 {base_dir}/scripts/validate_catalog.py "
    f"--known-pathogenic-loci-json-path {known_disease_associated_loci_source_path} "
    f"--check-for-presence-of-annotations "
    f"--check-for-presence-of-all-known-loci "
    f"{annotated_catalog_path}", step_number=29)

# copy files to the release_draft folder and compute catalog stats
updated_release_files = run_all_steps_helpers.validate_compress_and_index_release_files(
    release_files, base_dir, run, validate_json_step_number=30, tabix_step_number=31)

for path in updated_release_files:
    run(f"cp {path} {release_draft_folder}", step_number=32)

run(f"python3 -m str_analysis.compute_catalog_stats --reference-fasta {args.hg19_reference_fasta} --verbose {annotated_catalog_path}", step_number=33)

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
