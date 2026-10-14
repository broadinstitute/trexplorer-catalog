import argparse
import datetime
import gzip
import json
import os
import re
import requests
import subprocess
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

from str_analysis.utils.file_utils import file_exists, tee_stdout_and_stderr_to_log_file


VERSION = "2.1"

def run(command, step_number=None):
    command = re.sub("[ \\t]{2,}", "  ", command)  # remove extra spaces
    if args.dry_run:
        if step_number is not None:
            print(f"STEP #{step_number}: {command}")
        else:
            print(command)

    if not args.dry_run or command.startswith("mkdir"):
        if step_number is None or (
            not (args.only_step is not None and step_number != args.only_step) and
            not (args.start_with_step is not None and step_number < args.start_with_step) and
            not (args.end_with_step is not None and step_number > args.end_with_step)
        ):
            if step_number is not None:
                print(f"\nSTEP #{step_number}: {command}")
            else:
                print(f"\n--- $ {command}")
            subprocess.run(command, shell=True, check=True)

def chdir(d):
    print(f"cd {d}")
    os.chdir(d)

def is_bgzipped(path):
    """Returns True if path is bgzipped, and so can be tabix-indexed, rather than plain gzipped.

    A bgzipped file starts with a gzip header that sets the FEXTRA flag and carries a "BC"
    extra subfield. Checking the bytes means this works for a file that a previous run
    already compressed, where the branch that compressed it is no longer known.
    """
    with open(path, "rb") as f:
        header = f.read(14)
    return header[:4] == b"\x1f\x8b\x08\x04" and header[12:14] == b"BC"

# The LongTR, HipSTR and GangSTR converters in str-analysis write start_0based + 1, so those three
# formats are 1-based rather than true BED. Indexing them with "-p bed", which means "-s1 -b2 -e3 -0",
# shifts every record one base to the right and makes tabix miss it at the left edge of a query.
ONE_BASED_BED_SUFFIXES = (".LongTR.bed.gz", ".HipSTR.bed.gz", ".GangSTR.bed.gz")

def tabix_coordinate_flags(path):
    """Returns the tabix flags describing path's coordinate convention."""
    return "-s1 -b2 -e3" if path.endswith(ONE_BASED_BED_SUFFIXES) else "-p bed"

def get_remote_file_etag(url):
    """Returns the ETag that the server reports for url. It changes whenever the file's contents change."""
    if url.startswith("gs://"):
        return subprocess.run(["gcloud", "storage", "objects", "describe", url, "--format=value(etag)"],
                              capture_output=True, text=True, check=True).stdout.strip()

    response = requests.head(url, allow_redirects=True, timeout=60)
    response.raise_for_status()
    if not response.headers.get("ETag"):
        raise ValueError(f"{url} has no ETag, so there's no way to tell whether a local copy of it is up to date")
    return response.headers["ETag"]

def download_unless_local_copy_is_up_to_date(url, local_dir):
    """Downloads url into local_dir, unless the copy already there has the file's current ETag.

    The ETag of each download is saved next to it in a .etag file, which is what later runs compare against.

    Args:
        url (str): gs:// path or URL of the file
        local_dir (str): directory to download the file into

    Returns:
        str: path of the local copy
    """
    local_path = os.path.join(local_dir, os.path.basename(url))
    if args.dry_run:
        print(f"Would download {url} to {local_path} unless the local copy is up to date")
        return local_path

    etag_path = f"{local_path}.etag"
    remote_etag = get_remote_file_etag(url)
    if os.path.isfile(local_path) and os.path.isfile(etag_path):
        with open(etag_path) as f:
            if f.read() == remote_etag:
                print(f"Local copy of {url} is up to date: {local_path}")
                return local_path

    if url.startswith("gs://"):
        run(f"gcloud storage cp {url} {local_path}.tmp")
    else:
        run(f"wget -q -O {local_path}.tmp {url}")
    run(f"mv {local_path}.tmp {local_path}")
    with open(etag_path, "w") as f:
        f.write(remote_etag)

    return local_path



parser = argparse.ArgumentParser()
parser.add_argument("--hg38-reference-fasta", default="hg38.fa", help="Path of hg38 reference genome FASTA file")
parser.add_argument("--gencode-gtf", default="gs://str-truth-set/hg38/ref/other/gencode.v50.annotation.gtf.gz",
                    help="Gencode GTF file used to fill in the catalog's Gencode gene annotation fields")
parser.add_argument("--output-prefix", default=f"TRExplorer.repeat_catalog_v{VERSION}.hg38")
parser.add_argument("--variation-clusters-output-prefix", default=f"TRExplorer.variation_clusters_and_all_repeats_v{VERSION}.hg38")
parser.add_argument("--legacy-variation-clusters-output-prefix",
                    default=f"TRExplorer.variation_clusters_and_isolated_TRs_v{VERSION}.hg38",
                    help="Output prefix for the superseded catalog that holds a row only for "
                         "repeats outside a variation cluster. It ships alongside its "
                         "replacement for one release, so remove this option and the two steps "
                         "that use it once consumers have migrated.")

parser.add_argument("--only-step", type=int, help="Only run this one step")
parser.add_argument("--start-with-step", type=int, help="Start with a specific step number")
parser.add_argument("--end-with-step", type=int, help="End with a specific step number")
parser.add_argument("--tenk10k-annotations", default="gs://tandem-repeat-catalog/v2.0/tenk10k_str_mt_rows.reformatted.tsv.gz",
                    help="Path of the TenK10K annotations data")
parser.add_argument("--skip-tenk10k-annotations", action="store_true",
                    help="Skip adding TenK10K annotations to the catalog")
parser.add_argument("--hprc256-annotations", default="gs://tandem-repeat-catalog/v2.0/hprc_lps.2025_12.per_locus_and_motif.256_samples.tsv.gz",
                    help="Path of the HPRC256 annotations data")
parser.add_argument("--skip-hprc256-annotations", action="store_true",
                    help="Skip adding HPRC256 annotations to the catalog")
parser.add_argument("--aou1027-annotations", default="gs://tandem-repeat-catalog/v2.0/AoULR_phase1_TRGT_Weisburd_v1_combined.txt.gz",
                    help="Path of the AoU1027 annotations data")
parser.add_argument("--skip-aou1027-annotations", action="store_true",
                    help="Skip adding AoU1027 annotations to the catalog")
parser.add_argument("--variation-clusters-tsv", default="gs://tandem-repeat-catalog/v2.0/genome_clusters-v2.tsv.gz",
                    help="Variation clusters TSV file shared by Egor Dolzhenko")
parser.add_argument("--skip-variation-cluster-annotations", action="store_true",
                    help="Skip adding variation cluster annotations to the catalog")
parser.add_argument("--extended-loci-catalog",
                    default="gs://tandem-repeat-catalog/v2.1/TRExplorer.extended_loci_v2.1.hg38.bed.gz",
                    help="Catalog of wider locus definitions produced by running the gap-purity "
                         "extension rule over the v2 catalog. Generated once, by hand, by "
                         "extended_loci/generate_extended_loci_catalog.sh, and published alongside "
                         "the release so this pipeline is reproducible without regenerating it")
parser.add_argument("--timestamp", default=datetime.datetime.now().strftime('%Y-%m-%d'),
                    help="Timestamp to use in the output directory name")
parser.add_argument("--dry-run", action="store_true", help="Print commands without running them")

args = parser.parse_args()

print("TIMESTAMP:", args.timestamp)

base_dir = os.path.abspath(".")

# Remote inputs that aren't specific to one run are downloaded here. Every run shares this folder, so a file is
# only downloaded again after its version on the server changes.
downloaded_data_dir = os.path.join(base_dir, "downloaded_data")
run(f"mkdir -p {downloaded_data_dir}")

for key in "hg38_reference_fasta", "gencode_gtf", "extended_loci_catalog", "variation_clusters_tsv", "tenk10k_annotations", "hprc256_annotations", "aou1027_annotations":
    if (key == "variation_clusters_tsv" and args.skip_variation_cluster_annotations
        ) or (key == "tenk10k_annotations" and args.skip_tenk10k_annotations
        ) or (key == "hprc256_annotations" and args.skip_hprc256_annotations
        ) or (key == "aou1027_annotations" and args.skip_aou1027_annotations
    ):
        setattr(args, key, None)
        continue
    path = getattr(args, key)
    if not file_exists(path):
        parser.error(f"{key} file not found {path}")

    if path.startswith("gs://"):
        path = download_unless_local_copy_is_up_to_date(path, downloaded_data_dir)

    setattr(args, key, os.path.abspath(path))

working_dir = os.path.abspath(f"results__{args.timestamp}")

run(f"mkdir -p {working_dir}")
chdir(working_dir)

# create a release draft folder
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

# STEP #1  (already done for hg38, and results are publicly available. Will just download them in step #2)
#git clone git@github.com:broadinstitute/colab-repeat-finder.git
#cd colab-repeat-finder/python
#python3 perfect_repeat_finder.py  --min-repeats 3  --min-span 9  --min-motif-size 1  --max-motif-size 1000  --output-prefix perfect_repeats.hg38  --show-progress-bar  ${REFERENCE_FASTA_PATH}


# STEP #2:  generate catalog

# list of source catalogs, in order. If a locus is defined in more than one catalog (ie. overlapping boundaries,
# same motif), then the definition in the catalog that's earlier in the list will take precedence over definitions in
# subsequent catalogs.
source_catalogs_in_order = [
    NS(name="TRExplorerV1:KnownDiseaseAssociatedLoci",      merge="v1", url="https://raw.githubusercontent.com/broadinstitute/str-analysis/69dd90ecbc1dcbb23d5ca84ab4022850a283114f/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json"),
    NS(name="TRExplorerV1:Illumina174kPolymorphicTRs",      merge="v1", url="https://storage.googleapis.com/str-truth-set/hg38/ref/other/illumina_variant_catalog.sorted.bed.gz"),
    NS(name="TRExplorerV1:PerfectRepeatsInReference",       merge="v1", url="https://storage.googleapis.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp.bed.gz"),
    NS(name="TRExplorerV1:PolymorphicTRsInT2TAssemblies",   merge="v1", url="https://storage.googleapis.com/str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers_keeping_loci_that_have_overlapping_variants/combined/merged_expansion_hunter_catalog.78_samples.json.gz"),

    NS(name="TRExplorerV2:KnownDiseaseAssociatedLociV2",    merge="v2a", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/known_disease_associated_loci_v2.loci_to_include_in_catalog.bed.gz"),
    NS(name="TRExplorerV2:PolymorphicTRsInT2TAssembliesV2", merge="v2a", url="https://storage.googleapis.com/str-truth-set-v2/filter_vcf_v2__2025_12_29/combined.321_catalogs.tandem_repeats.bed.gz"),
    NS(name="TRExplorerV2:KnownFunctionalVNTRs",            merge="v2a", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/known_functional_VNTRs.loci_to_include_in_catalog.bed.gz"),

    NS(name="TRExplorerV2:HipSTRCatalog"        ,           merge="v2b", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/hg38.hipstr_reference.loci_to_include_in_catalog.bed.gz"),
    NS(name="TRExplorerV2:AdottoTRsFromDanzi2025",          merge="v2b", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/Danzi_2025.adotto.loci_to_include_in_catalog.bed.gz"),
    NS(name="TRExplorerV2:ClinvarIndelsThatAreTRs2025",     merge="v2b", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/clinvar_2025_11_03.loci_to_include_in_catalog.bed.gz"),
    NS(name="TRExplorerV2:Tanudisastro2025",                merge="v2b", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/Tanudisastro_2025.loci_to_include_in_catalog.bed.gz"),
    NS(name="TRExplorerV2:Manigbas2024",                    merge="v2b", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/Manigbas_2024_genotyped_TRs.bed.gz"),
    NS(name="TRExplorerV2:Sulovari2021",                    merge="v2b", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/Sulovari_2019_human_specific_STRs.bed.gz"),
    NS(name="TRExplorerV2:Garg2021",                        merge="v2b", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/Garg_2021_eVNTRs_mVNTRs.bed.gz"),
    NS(name="TRExplorerV2:Mukamel2021",                     merge="v2b", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/Mukamel_2021.loci_to_include_in_catalog.bed.gz"),
    NS(name="TRExplorerV2:Annear2021",                      merge="v2b", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/Annear_2021.loci_to_include_in_catalog.bed.gz"),
    NS(name="TRExplorerV2:Gymrek2016",                      merge="v2b", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/Gymrek_2016_eSTRs.tandem_repeats.bed.gz"),
    NS(name="TRExplorerV2:Hause2016",                       merge="v2b", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/Hause_2016_cancer_MSI_loci.bed.gz"),
    NS(name="TRExplorerV2:VamosV3",                         merge="v2b", url="https://storage.googleapis.com/tandem-repeat-catalog/v2.0/vamosGenomicTR_v3.0.loci_to_include_in_catalog.bed.gz"),    
]

source_catalog_paths = {}
for source_catalog in source_catalogs_in_order:
    source_catalog_paths[source_catalog.name] = download_unless_local_copy_is_up_to_date(
        source_catalog.url, os.path.abspath("."))

# Remote files other than the source catalogs that later steps read. They're downloaded here, together with
# everything else the pipeline needs, so that none of the steps need internet access.
additional_data_urls = [
    NS(name="MappabilityTrack",              url="gs://tgg-viewer/ref/GRCh38/mappability/GRCh38_no_alt_analysis_set_GCA_000001405.15-k36_m2.bw"),
    NS(name="KnownDiseaseAssociatedLoci",    url="https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json"),
    NS(name="RefseqGTF",                     url="gs://str-truth-set/hg38/ref/other/hg38.ncbiRefSeq.gtf.gz"),
    NS(name="ManeGTF",                       url="gs://str-truth-set/hg38/ref/other/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz"),
    NS(name="Illumina174kAlleleHistograms",  url="https://github.com/Illumina/RepeatCatalogs/raw/master/hg38/genotype/1000genomes/1kg.gt.hist.tsv.gz"),
    NS(name="T2TAssemblyGenotypes",          url="gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers_keeping_loci_that_have_overlapping_variants/combined/joined.78_samples.variants.tsv.gz"),
    NS(name="TRExplorerCatalogV1",           url="https://github.com/broadinstitute/tandem-repeat-catalog/releases/download/v1.0/repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.json.gz"),
    NS(name="GangSTRCatalogV17",             url="https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver17.bed.gz"),
]

additional_data_paths = {}
for additional_data in additional_data_urls:
    additional_data_paths[additional_data.name] = download_unless_local_copy_is_up_to_date(
        additional_data.url, downloaded_data_dir)

print("Done downloading. The remaining steps don't need internet access.")

# preprocess catalog of known disease-associated loci: split compound definitions
run(f"python3 -u -m str_analysis.split_adjacent_loci_in_expansion_hunter_catalog {source_catalog_paths['TRExplorerV1:KnownDiseaseAssociatedLoci']}", step_number=1)
source_catalog_paths['TRExplorerV1:KnownDiseaseAssociatedLoci'] = source_catalog_paths['TRExplorerV1:KnownDiseaseAssociatedLoci'].replace(".json", ".split.json")
# change motif definition for the RFC1 locus from AARRG => AAAAG since our catalog doesn't currently support IUPAC codes
run(f"sed -i 's/AARRG/AAAAG/g' {source_catalog_paths['TRExplorerV1:KnownDiseaseAssociatedLoci']}", step_number=2)
run(f"gzip -f {source_catalog_paths['TRExplorerV1:KnownDiseaseAssociatedLoci']}", step_number=3)

source_catalog_paths['TRExplorerV1:KnownDiseaseAssociatedLoci'] += ".gz"

if not args.dry_run:
    # compute stats for primary disease-associated loci
    with gzip.open(source_catalog_paths["TRExplorerV1:KnownDiseaseAssociatedLoci"]) as f:
        known_disease_associated_loci = json.load(f)

    primary_disease_associated_loci = [
        x for x in known_disease_associated_loci if x["Diseases"] and (
            x["LocusId"].startswith("HOXA") or x["LocusId"].startswith("ARX") or "_" not in x["LocusId"]
        )
    ]


    # check the number of primary disease-associated loci, not counting adjacent repeats and historic candidate loci that
    # are not currently considered monogenic
    if len(primary_disease_associated_loci) != 63:
        raise ValueError(f"Expected 63 primary disease-associated loci, found {len(primary_disease_associated_loci)}")

primary_disease_associated_loci_path = source_catalog_paths["TRExplorerV1:KnownDiseaseAssociatedLoci"].replace(
    ".json.gz", ".primary_disease_associated_loci.json.gz")

if not args.dry_run:
    with gzip.open(primary_disease_associated_loci_path, "wt") as f:
        json.dump(primary_disease_associated_loci, f, indent=4)

    run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog \
        --verbose \
        --reference-fasta {args.hg38_reference_fasta} \
        --min-interval-size-bp 1 \
        --skip-gene-annotations \
        --skip-mappability-annotations \
        --skip-disease-loci-annotations \
        --discard-loci-with-non-ACGT-bases-in-reference \
        --discard-loci-with-non-ACGTN-bases-in-motif \
        --output-path {primary_disease_associated_loci_path} \
        {primary_disease_associated_loci_path}""", step_number=4)

    run(f"python3 -m str_analysis.compute_catalog_stats --reference-fasta {args.hg38_reference_fasta} --verbose {primary_disease_associated_loci_path}", step_number=5)

adjacent_repeats_source_bed = None
for motif_size_label, min_motif_size, max_motif_size, release_tar_gz_path in [
    ("1_to_1000bp_motifs",  1, 1000, None),
    #("2_to_1000bp_motifs",  2, 1000, f"{args.output_prefix}.subset.all_loci_except_homopolymers.tar.gz"),
    #("homopolymers",        1, 1,    f"{args.output_prefix}.subset.only_homopolymer_loci.tar.gz"),
    #("2_to_6bp_motifs",     2, 6,    f"{args.output_prefix}.subset.only_loci_with_2_to_6bp_motifs.tar.gz"),
    #("7_to_1000bp_motifs",  7, 1000, f"{args.output_prefix}.subset.only_loci_with_7_to_1000bp_motifs.tar.gz"),
]:

    if (args.start_with_step and args.start_with_step > 2) and motif_size_label != "1_to_1000bp_motifs":
        continue

    print("="*200)
    chdir(working_dir)
    run(f"mkdir -p {motif_size_label}")
    chdir(motif_size_label)

    output_prefix = os.path.abspath(f"{args.output_prefix}.{motif_size_label}")

    filtered_source_catalog_paths = {}
    for i, source_catalog in enumerate(source_catalogs_in_order):
        catalog_path = source_catalog_paths[source_catalog.name]

        if catalog_path.endswith(".json.gz"):
            filtered_catalog_path = catalog_path.replace(".json.gz", ".filtered.json.gz")
        elif catalog_path.endswith(".bed.gz"):
            filtered_catalog_path = catalog_path.replace(".bed.gz", ".filtered.json.gz")
        else:
            raise ValueError(f"Unexpected file extension for {catalog_path}")

        filtered_catalog_path = os.path.abspath(os.path.basename(filtered_catalog_path))
        filtered_source_catalog_paths[source_catalog.name] = filtered_catalog_path

        extra_args = ""
        if source_catalog.name == "TRExplorerV2:PolymorphicTRsInT2TAssembliesV2":
            extra_args = "--min-repeats-in-reference 1 "
        if source_catalog.name.startswith("TRExplorerV2:"):
            extra_args += "--discard-loci-with-non-ACGT-bases-in-motif "

        run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog \
            --verbose \
            --known-disease-associated-loci {primary_disease_associated_loci_path} \
            --reference-fasta {args.hg38_reference_fasta} \
            {extra_args} \
            --min-motif-size {min_motif_size} \
            --max-motif-size {max_motif_size} \
            --min-interval-size-bp 1 \
            --skip-gene-annotations \
            --skip-mappability-annotations \
            --skip-disease-loci-annotations \
            --set-locus-id \
            --discard-loci-with-non-ACGT-bases-in-reference \
            --discard-loci-with-non-ACGTN-bases-in-motif \
            --output-path {filtered_catalog_path} \
            {catalog_path}""",
            step_number=6)

        print(f"Stats for {catalog_path}")
        run(f"python3 -m str_analysis.compute_catalog_stats "
            f"--reference-fasta {args.hg38_reference_fasta} "
            f"--verbose "
            f"{filtered_catalog_path}",
            step_number=7)

    # NOTE: we don't use the --merge-adjacent-loci-with-same-motif option for str_analysis.merge_loci because
    # it's important to preserve locus definitions as they appear in the individual source catalogs. If loci
    # are merged together, this can create incompatibility with loci in source catalogs (ie. the illumina catalog) or
    # between future versions of the overall Simple Repeat Catalog since newly-added loci could merge with previously
    # added loci, causing loss of those loci from new versions of the catalog. Variation cluster analysis is a better
    # way to merge adjacent loci where needed.
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
        --write-outer-join-table \
        --write-bed-files-with-unique-loci \
        --outer-join-overlap-table-min-sources 1 \
        --output-prefix {output_prefix}.merged_v1 \
        {source_catalog_paths_for_merge_command_v1}""", step_number=8)

    run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog --verbose \
        --reference-fasta {args.hg38_reference_fasta} \
        --known-disease-associated-loci {additional_data_paths['KnownDiseaseAssociatedLoci']} \
        --mappability-track-bigwig {additional_data_paths['MappabilityTrack']} \
        --min-motif-size {min_motif_size} \
        --max-motif-size {max_motif_size} \
        --min-interval-size-bp 1 \
        --discard-overlapping-intervals-with-similar-motifs \
        --output-path {output_prefix}.merged_v1.annotated.json.gz \
        {output_prefix}.merged_v1.json.gz""", step_number=9)

    source_catalog_paths_for_merge_command_v2a = " ".join([
        f"merged:{output_prefix}.merged_v1.annotated.json.gz",
    ] + [
        f"{source_catalog.name}:{filtered_source_catalog_paths[source_catalog.name]}"
        for source_catalog in source_catalogs_in_order if source_catalog.merge == "v2a"
    ])

    run(f"""python3 -u -m str_analysis.merge_loci --verbose \
        --add-found-in-fields \
        --output-format JSON \
        --overlapping-loci-action keep-first \
        --write-merge-stats-tsv \
        --overlap-fraction 0.66 \
        --motif-length-match-sufficient-for-VNTRs \
        --min-jaccard-similarity 0.66 \
        --output-prefix {output_prefix}.merged_v2a \
        {source_catalog_paths_for_merge_command_v2a}""", step_number=10)

    run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog --verbose \
        --reference-fasta {args.hg38_reference_fasta} \
        --known-disease-associated-loci {additional_data_paths['KnownDiseaseAssociatedLoci']} \
        --mappability-track-bigwig {additional_data_paths['MappabilityTrack']} \
        --min-motif-size {min_motif_size} \
        --max-motif-size {max_motif_size} \
        --min-interval-size-bp 1 \
        --output-path {output_prefix}.merged_v2a.annotated.json.gz \
        {output_prefix}.merged_v2a.json.gz""", step_number=11)

    source_catalog_paths_for_merge_command_v2b = " ".join([
       f"merged:{output_prefix}.merged_v2a.annotated.json.gz",
    ] + [
       f"{source_catalog.name}:{filtered_source_catalog_paths[source_catalog.name]}"
        for source_catalog in source_catalogs_in_order if source_catalog.merge == "v2b"
    ])

    run(f"""python3 -u -m str_analysis.merge_loci --verbose \
        --add-found-in-fields \
        --output-format JSON \
        --overlapping-loci-action keep-first \
        --write-merge-stats-tsv \
        --overlap-fraction 0.66 \
        --motif-length-match-sufficient-for-VNTRs \
        --min-jaccard-similarity 0.2 \
        --output-prefix {output_prefix}.merged_v2b \
        {source_catalog_paths_for_merge_command_v2b}""", step_number=12)

    run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog --verbose \
        --reference-fasta {args.hg38_reference_fasta} \
        --known-disease-associated-loci {additional_data_paths['KnownDiseaseAssociatedLoci']} \
        --mappability-track-bigwig {additional_data_paths['MappabilityTrack']} \
        --min-motif-size {min_motif_size} \
        --max-motif-size {max_motif_size} \
        --min-interval-size-bp 1 \
        --output-path {output_prefix}.merged_v2b.annotated.json.gz \
        {output_prefix}.merged_v2b.json.gz""", step_number=13)

    source_catalog_paths_for_final_merge_command = " ".join([
        f"merged_v1:{output_prefix}.merged_v1.annotated.json.gz",
        f"merged_v2a:{output_prefix}.merged_v2a.annotated.json.gz",
        f"merged_v2b:{output_prefix}.merged_v2b.annotated.json.gz",
    ])

    run(f"""python3 -u -m str_analysis.merge_loci --verbose \
        --add-found-in-fields \
        --output-format JSON \
        --overlap-fraction 0.66 \
        --min-jaccard-similarity 0.66 \
        --motif-length-match-sufficient-for-VNTRs \
        --only-compare-loci-from-different-catalogs \
        --overlapping-loci-action keep-first \
        --write-merge-stats-tsv \
        --output-prefix {output_prefix}.final_merged \
        {source_catalog_paths_for_final_merge_command}""", step_number=14)

    # Add the extended-boundary loci. These are a separate source of loci rather than something this
    # pipeline computes: extended_loci/generate_extended_loci_catalog.sh runs the gap-purity rule over
    # the finished v2 catalog once, by hand, and keeps only the wider definitions it produced. They
    # are appended rather than merged because a locus and the wider definition that grew out of it
    # necessarily overlap, and merge_loci would discard one of the two. The rule already dropped any
    # definition that v2 described, so nothing appended here duplicates what is already present.
    # Sorting matters: the annotation steps below group overlapping loci and reject a catalog whose
    # coordinates ever decrease, and an extended definition can start before the locus it grew from.
    #
    # Known trade-off, accepted rather than overlooked. These definitions carry locus ids generated
    # from their own coordinates, so none of them appears in the variation clusters TSV, which was
    # computed against v2's ids. They therefore get no VariationCluster annotations at the step
    # below, and the TRGT variation-cluster catalog writes each as its own row. Of the 58,135, only
    # 7,844 sit wholly inside a cluster row that already covers them; 22,089 overlap a cluster but
    # reach past its boundary, and 28,202 fall outside every cluster. This is the same path every
    # locus takes when its source postdates the cluster table: all 56,846 TRExplorerV2:VamosV3 loci
    # in v2 are already in exactly this position. Fixing it properly means either regenerating the
    # variation clusters against this catalog, or teaching the TRGT catalog script to skip loci
    # whose Source is TRExplorerV2.1:ExtendedLocusBoundaries when restoring what the TSV lacks.
    run(f"""python3 -u << EOF
from str_analysis.extend_str_catalog_boundaries import sort_catalog_records
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator
from str_analysis.utils.export_json import export_json

records = list(get_variant_catalog_iterator("{output_prefix}.final_merged.json.gz"))
print(f"Read {{len(records):,d}} loci from {output_prefix}.final_merged.json.gz")

extended_records = list(get_variant_catalog_iterator("{args.extended_loci_catalog}"))
print(f"Read {{len(extended_records):,d}} extended loci from {args.extended_loci_catalog}")

already_in_catalog = {{
    (r["ReferenceRegion"], r["LocusStructure"]) for r in records
}}
new_records = [
    r for r in extended_records
    if (r["ReferenceRegion"], r["LocusStructure"]) not in already_in_catalog
]
if len(new_records) < len(extended_records):
    print(f"Skipped {{len(extended_records) - len(new_records):,d}} extended loci already defined in "
          f"the catalog")

# Name them as their own source. Without this they reach the sources-stats table as "Unknown", since
# it reads record["Source"], and the catalog would not say where a fifth of its widest definitions
# came from.
for record in new_records:
    record["Source"] = "TRExplorerV2.1:ExtendedLocusBoundaries"
    record["FoundInTRExplorerV2.1:ExtendedLocusBoundaries"] = True

export_json(sort_catalog_records(records + new_records),
            "{output_prefix}.final_merged.with_extended_loci.json.gz")
print(f"Wrote {{len(records) + len(new_records):,d}} loci")
EOF
""", step_number=15)

    annotated_catalog_path = f"{output_prefix}.EH.with_annotations.json.gz"
    run(f"""python3 -u -m str_analysis.annotate_and_filter_str_catalog --verbose \
        --reference-fasta {args.hg38_reference_fasta} \
        --gene-models-source gencode \
        --gencode-gtf {args.gencode_gtf} \
        --gene-models-source refseq \
        --refseq-gtf {additional_data_paths['RefseqGTF']} \
        --gene-models-source mane \
        --mane-gtf {additional_data_paths['ManeGTF']} \
        --mappability-track-bigwig {additional_data_paths['MappabilityTrack']} \
        --known-disease-associated-loci {additional_data_paths['KnownDiseaseAssociatedLoci']} \
        --min-motif-size {min_motif_size} \
        --max-motif-size {max_motif_size} \
        --min-interval-size-bp 1 \
        --output-path {annotated_catalog_path} \
        {output_prefix}.final_merged.with_extended_loci.json.gz""", step_number=16)

    # The annotation above labels a known disease-associated locus only when the catalog locus has exactly the same
    # boundaries, which misses loci like PRE-MIR7-2 (chr15:88569433-88569449, while the catalog locus also includes
    # the trailing partial repeat: 15-88569433-88569452-TTTG). Label those using the TRExplorerV2 field of the known
    # disease-associated loci catalog. These commands share step 16's number so that rerunning step 16 reapplies them.
    run(f"""python3 {base_dir}/scripts/add_known_disease_associated_locus_labels_using_TRExplorerV2_ids.py \
        --known-disease-associated-loci {additional_data_paths['KnownDiseaseAssociatedLoci']} \
        --output-catalog-json-path {annotated_catalog_path}.with_TRExplorerV2_disease_locus_labels.json.gz \
        {annotated_catalog_path}""", step_number=16)
    run(f"mv {annotated_catalog_path}.with_TRExplorerV2_disease_locus_labels.json.gz {annotated_catalog_path}", step_number=16)

    # create a version of the ExpansionHunter catalog without extra annotations
    run(f"""python3 << EOF
import gzip, ijson, json

f = gzip.open("{annotated_catalog_path}", "rt")
out = gzip.open("{output_prefix}.EH.json.gz", "wt")
i = 0
total = 0
out.write("[")
for record in ijson.items(f, "item", use_float=True):
    total += 1
    # skip chrM loci because ExpansionHunter prints an error like: 'Unable to extract chrM:-793-207 from hg38.fa'
    is_chrM = record["LocusId"].startswith("M-") or record["LocusId"].startswith("chrM-")
    if is_chrM: print(f"Skipping chrM locus: {{record['LocusId']}}")
    if is_chrM: continue
    if record["NsInFlanks"] > 5: continue
    if i > 0: out.write(", ")
    i += 1
    out.write(json.dumps({{ 
        k: v for k, v in record.items() if k in {{"LocusId", "ReferenceRegion", "VariantType", "LocusStructure"}} 
    }}, indent=4))
out.write("]")
print(f"Wrote {{i:,d}} out of {{total:,d}} ({{i/total:0.1%}}%) loci to {output_prefix}.EH.json.gz")
EOF
""", step_number=17)

    run(f"python3 -m str_analysis.filter_out_loci_with_Ns_in_flanks "
        f"-R {args.hg38_reference_fasta} "
        f"-o {output_prefix}.EH.without_loci_with_Ns_in_flanks.json.gz "
        f"--output-list-of-filtered-loci {output_prefix}.loci_with_Ns_in_flanks.txt "
        f"{output_prefix}.EH.json.gz", step_number=18)
    run(f"mv {output_prefix}.EH.without_loci_with_Ns_in_flanks.json.gz {output_prefix}.EH.json.gz", step_number=19)

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


    # add variation cluster annotations to the catalog
    if args.variation_clusters_tsv and not args.skip_variation_cluster_annotations:
        variation_clusters_and_all_repeats_release_filename = f"{args.variation_clusters_output_prefix}.TRGT.bed.gz"

        run(f"""python3 {base_dir}/scripts/add_variation_cluster_annotations_to_catalog.py \
            --verbose \
            --output-catalog-json-path {output_prefix}.EH.with_annotations.with_variation_clusters.json.gz \
            {args.variation_clusters_tsv} \
            {annotated_catalog_path}""", step_number=20)

        run(f"mv {output_prefix}.EH.with_annotations.with_variation_clusters.json.gz {annotated_catalog_path}", step_number=21)

        run(f"""python3 {base_dir}/scripts/generate_TRGT_catalog_with_variation_clusters_and_all_repeats.py \
            -o {variation_clusters_and_all_repeats_release_filename} \
            {args.variation_clusters_tsv} \
            {annotated_catalog_path}""", step_number=22)

        release_files.append(variation_clusters_and_all_repeats_release_filename)

        # --name-clusters-by-trid only for this catalog: it is the one that also gives every repeat
        # a row of its own, so a single-repeat cluster named after its member would collide with it.
        run(f"python3 {base_dir}/scripts/convert_trgt_catalog_to_longtr_format.py --name-clusters-by-trid {variation_clusters_and_all_repeats_release_filename}",
            step_number=23)

        release_files.append(variation_clusters_and_all_repeats_release_filename.replace(".TRGT.bed.gz", ".LongTR.bed.gz"))

        # trgt reads a bgzipped BED directly, so there is no need to decompress it first. This
        # catalog holds a row per repeat as well as a row per variation cluster, which means rows
        # deliberately overlap, and this step is what confirms trgt still accepts them.
        run(f"trgt validate --genome {args.hg38_reference_fasta} --repeats {variation_clusters_and_all_repeats_release_filename}",
            step_number=24)

        # The superseded catalog, which holds a row only for repeats that are not inside a variation
        # cluster. It ships alongside its replacement for one release so consumers have time to
        # migrate. Delete this block, and the --legacy-variation-clusters-output-prefix option, after
        # that release.
        legacy_variation_clusters_release_filename = f"{args.legacy_variation_clusters_output_prefix}.TRGT.bed.gz"

        run(f"""python3 {base_dir}/scripts/generate_TRGT_catalog_with_isolated_repeats_and_variation_clusters.py \
            -o {legacy_variation_clusters_release_filename} \
            {args.variation_clusters_tsv} \
            {annotated_catalog_path}""", step_number=25)

        release_files.append(legacy_variation_clusters_release_filename)

        run(f"python3 {base_dir}/scripts/convert_trgt_catalog_to_longtr_format.py {legacy_variation_clusters_release_filename}",
            step_number=26)

        release_files.append(legacy_variation_clusters_release_filename.replace(".TRGT.bed.gz", ".LongTR.bed.gz"))

        # The two catalogs above are built by separate scripts from the same TSV and repeat catalog,
        # so nothing else keeps them in sync. This is what would have caught
        # https://github.com/PacificBiosciences/trgt-lps/issues/5#issuecomment-5683540748, where the
        # legacy catalog was missing repeats that the all_repeats catalog (and the HPRC TRGT VCF) had.
        run(f"""python3 {base_dir}/scripts/validate_variation_cluster_catalogs_are_consistent.py \
            {variation_clusters_and_all_repeats_release_filename} \
            {legacy_variation_clusters_release_filename}""", step_number=27)

    # annotate overlapping loci and output merged tandem repeat regions
    run(f"""python3 {base_dir}/scripts/annotate_overlapping_loci_and_output_merged_tandem_repeat_regions.py \
        --verbose \
        -o {annotated_catalog_path}.with_overlap_annotations.json.gz \
        {annotated_catalog_path}""", step_number=28)
    run(f"mv {annotated_catalog_path}.with_overlap_annotations.json.gz {annotated_catalog_path}", step_number=29)

    # add allele frequencies to the catalog
    run(f"""python3 -u {base_dir}/scripts/add_allele_frequency_annotations.py \
            --illumina174k-allele-histograms-tsv {additional_data_paths['Illumina174kAlleleHistograms']} \
            --t2t-assembly-genotypes-tsv {additional_data_paths['T2TAssemblyGenotypes']} \
            --add-t2t-assembly-frequencies-to-overlapping-loci \
            -o {annotated_catalog_path}.with_allele_frequencies.json.gz  {annotated_catalog_path}""", step_number=30)

    run(f"mv {annotated_catalog_path}.with_allele_frequencies.json.gz {annotated_catalog_path}", step_number=31)

    # Add TenK10K population statistics annotations
    if args.tenk10k_annotations:
        run(f"""python3 {base_dir}/scripts/add_TenK10K_annotations_to_catalog.py \
            --tsv-path {args.tenk10k_annotations} \
            -o {annotated_catalog_path}.with_TenK10K_annotations.json.gz \
            {annotated_catalog_path}""", step_number=32)
        run(f"mv {annotated_catalog_path}.with_TenK10K_annotations.json.gz {annotated_catalog_path}", step_number=33)

    # Add HPRC256 population statistics annotations
    if args.hprc256_annotations:
        run(f"""python3 {base_dir}/scripts/add_HPRC256_annotations_to_catalog.py \
            --tsv-path {args.hprc256_annotations} \
            -o {annotated_catalog_path}.with_HPRC256_annotations.json.gz \
            {annotated_catalog_path}""", step_number=34)
        run(f"mv {annotated_catalog_path}.with_HPRC256_annotations.json.gz {annotated_catalog_path}", step_number=35)

    # Add AoU1027 population statistics annotations
    if args.aou1027_annotations:
        run(f"""python3 {base_dir}/scripts/add_AoU_annotations_to_catalog.py \
            --tsv-path {args.aou1027_annotations} \
            -o {annotated_catalog_path}.with_AoU1027_annotations.json.gz \
            {annotated_catalog_path}""", step_number=36)
        run(f"mv {annotated_catalog_path}.with_AoU1027_annotations.json.gz {annotated_catalog_path}", step_number=37)

    # annotate with "TRsInRegion" based on adjacent loci
    if motif_size_label == "1_to_1000bp_motifs":
        adjacent_repeats_source_bed = f"{output_prefix}.bed.gz"

    # convert to BED
    run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_bed --split-adjacent-repeats "
        f"{annotated_catalog_path}  --output-file {output_prefix}.bed.gz", step_number=38)

    run(f"python3 -m str_analysis.add_adjacent_loci_to_expansion_hunter_catalog "
        f"--ref-fasta {args.hg38_reference_fasta} "
        f"--source-of-adjacent-loci {adjacent_repeats_source_bed} "
        f"--add-extra-field TRsInRegion "
        f"--only-add-extra-fields "
        f"-o {annotated_catalog_path}.with_adjacent_loci_annotation.json.gz "
        f"{annotated_catalog_path}", step_number=39)

    run(f"mv {annotated_catalog_path}.with_adjacent_loci_annotation.json.gz {annotated_catalog_path}", step_number=40)


    # convert to TSV
    run(f"""python3 << EOF
import gzip
import json
import pandas as pd
from pprint import pformat
with gzip.open("{annotated_catalog_path}", "rt") as f:
    data = json.load(f)
print(f"Loaded {{len(data):,d}} records from {annotated_catalog_path}")
print("Writing it to TSV file...")
df = pd.DataFrame(data)
core_columns = [
    'LocusId', 'ReferenceRegion', 'LocusStructure', 'CanonicalMotif', 'TRsInRegion',
    'Source', 'GencodeGeneRegion', 'GencodeGeneId', 'GencodeGeneName', 'GencodeTranscriptId',
    'RefseqGeneRegion', 'RefseqGeneId', 'RefseqGeneName', 'RefseqTranscriptId', 
    'ManeGeneRegion', 'ManeGeneId', 'ManeGeneName', 'ManeTranscriptId',
    'KnownDiseaseAssociatedMotif',  'KnownDiseaseAssociatedLocus', 'NsInFlanks',
    'LeftFlankMappability', 'FlanksAndLocusMappability', 'RightFlankMappability', 
    'FoundInTRExplorerV1:KnownDiseaseAssociatedLoci', 'FoundInTRExplorerV1:Illumina174kPolymorphicTRs',
    'FoundInTRExplorerV1:PerfectRepeatsInReference', 'FoundInTRExplorerV1:PolymorphicTRsInT2TAssemblies',
    'NumRepeatsInReference', 'ReferenceRepeatPurity', 
    'AlleleFrequenciesFromIllumina174k', 'StdevFromIllumina174k',
    'AlleleFrequenciesFromT2TAssemblies', 'StdevFromT2TAssemblies',
    'VariationCluster', 'VariationClusterMotifs', 'VariationClusterSizeDiff',
    #'AlleleFrequenciesFromHPRC256', 'StdevFromHPRC256',
    'AoU1027_Stdev', 'AoU1027_Median', 'AoU1027_MaxAllele',
]

drop_columns = ['VariantType', ]
for c in set(core_columns)  - set(df.columns): df[c] = None
df = df[core_columns + [c for c in df.columns if c not in (core_columns + drop_columns)]]

output_tsv_path = "{annotated_catalog_path.replace('.json.gz', '') + '.tsv.gz'}"
df.to_csv(output_tsv_path, sep="\\t", index=False)
print(f"Wrote {{len(df):,d}} rows to {{output_tsv_path}} with columns: {{pformat(list(df.columns))}}")
EOF
""", step_number=41)

    # convert the catalog from ExpansionHunter catalog format to TRGT, LongTR, HipSTR, and GangSTR formats
    run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_trgt_catalog -R {args.hg38_reference_fasta} --split-adjacent-repeats {annotated_catalog_path}  --output-file {output_prefix}.TRGT.bed", step_number=42)
    run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_longtr_format  {annotated_catalog_path}  --output-file {output_prefix}.LongTR.bed", step_number=43)
    run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_hipstr_format  {annotated_catalog_path}  --output-file {output_prefix}.HipSTR.bed", step_number=44)
    run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_gangstr_spec   {annotated_catalog_path}  --output-file {output_prefix}.GangSTR.bed", step_number=45)

    # ATaRVa (https://github.com/SowpatiLab/ATaRVa) reads a 5-column BED: chromosome, start, end,
    # motif, motif length. It reads one motif per row, hence --split-adjacent-repeats, and requires
    # the file to be bgzipped and tabix-indexed, which the converter does for a .gz output path.
    run(f"python3 -m str_analysis.convert_expansion_hunter_catalog_to_bed --split-adjacent-repeats "
        f"--motif-size-column {annotated_catalog_path}  --output-file {output_prefix}.ATaRVa.bed.gz",
        step_number=46)

    # Confirm that the TRGT catalog passes 'trgt validate'
    run(f"trgt validate --genome {args.hg38_reference_fasta}  --repeats {output_prefix}.TRGT.bed", step_number=47)

    # Check that the combined catalog's repeat rows still reproduce the repeat catalog. This runs
    # here rather than beside the steps that build the combined catalog because it compares against
    # {output_prefix}.TRGT.bed, which step 42 writes well after those steps. The guard repeats the
    # one on the block that defines the filename.
    if args.variation_clusters_tsv and not args.skip_variation_cluster_annotations:
        run(f"""python3 {base_dir}/scripts/validate_variation_clusters_catalog.py \
            {variation_clusters_and_all_repeats_release_filename} \
            {output_prefix}.TRGT.bed""", step_number=48)

    # Print missing values
    run(f"python3 {base_dir}/scripts/print_missing_values_percentages_in_json_or_tsv.py  {annotated_catalog_path}", step_number=49)

    # Perform basic internal consistency checks on the JSON catalog
    run(f"python3 {base_dir}/scripts/validate_catalog.py " +
        f"--known-pathogenic-loci-json-path {source_catalog_paths['TRExplorerV1:KnownDiseaseAssociatedLoci']} " +
        f"--catalog-v1-json-path {additional_data_paths['TRExplorerCatalogV1']} " +
        ("--check-for-presence-of-annotations --check-for-presence-of-all-known-loci --check-for-presence-of-all-loci-from-v1 " if motif_size_label == "1_to_1000bp_motifs" else "") +
        f"{annotated_catalog_path}", step_number=50)

    # copy files to the release_draft folder and compute catalog stats
    updated_release_files = []
    for path in release_files:
        if path.endswith(".bed"):
            updated_release_files.append(f"{path}.gz")
        else:
            if path.endswith(".json") or path.endswith(".json.gz") and ".EH." in path:
                run(f"python3 {base_dir}/scripts/validate_json.py -k LocusId -k LocusStructure -k ReferenceRegion -k VariantType {path}", step_number=51)
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
        run(f"{prepare}tabix -f {tabix_coordinate_flags(path)} {path}", step_number=52)
        updated_release_files.append(f"{path}.tbi")

    if release_tar_gz_path is None:
        for path in updated_release_files:
            run(f"cp {path} {release_draft_folder}", step_number=53)
    else:
        run(f"tar czf {release_tar_gz_path} -C {os.path.dirname(output_prefix)} " + " ".join([os.path.basename(p) for p in updated_release_files]), step_number=54)
        run(f"cp {release_tar_gz_path} {release_draft_folder}", step_number=55)

    run(f"python3 -m str_analysis.compute_catalog_stats --reference-fasta {args.hg38_reference_fasta} --verbose {annotated_catalog_path}", step_number=56)

    # Print source statistics table
    run(f"python3 {base_dir}/scripts/generate_catalog_sources_stats_table.py {release_draft_folder}/{os.path.basename(annotated_catalog_path)}", step_number=57)

    # report hours, minutes, seconds relative to script_start_time
    diff = time.time() - script_start_time
    print(f"Done generating {output_prefix} catalog. Took {diff//3600:.0f}h, {(diff%3600)//60:.0f}m, {diff%60:.0f}s")

    if motif_size_label != "1_to_1000bp_motifs":
        continue

    # compare to the GangSTR_v17 catalog to make sure it's included
    comparison_catalog_paths = {
        "GangSTR_v17": additional_data_paths["GangSTRCatalogV17"],
        #"vamos_catalog_v2.1": "https://storage.googleapis.com/str-truth-set/hg38/ref/other/vamos_catalog.v2.1.bed.gz",
    }

    # write the converted catalog here rather than next to the download, which is shared by every run
    path_after_conversion = os.path.abspath(os.path.basename(comparison_catalog_paths["GangSTR_v17"]).replace(".bed.gz", ".json.gz"))
    run(f"python3 -u -m str_analysis.convert_gangstr_spec_to_expansion_hunter_catalog --verbose {comparison_catalog_paths['GangSTR_v17']} -o {path_after_conversion}", step_number=58)
    comparison_catalog_paths["GangSTR_v17"] = path_after_conversion

    # compare catalog to other catalogs
    start_time = time.time()

    for catalog_name, path in comparison_catalog_paths.items():
        filtered_comparison_catalog_path = re.sub("(.json|.bed)(.gz)?$", "", path) + ".filtered.json.gz"

        run(f"""python3 -m str_analysis.annotate_and_filter_str_catalog \
            --reference-fasta {args.hg38_reference_fasta} \
            --skip-gene-annotations \
            --skip-mappability-annotations \
            --skip-disease-loci-annotations \
            --min-motif-size {min_motif_size} \
            --max-motif-size {max_motif_size} \
            --output-path {filtered_comparison_catalog_path} \
            --verbose \
            {path}""", step_number=59)

        run(f"python3 -m str_analysis.compute_catalog_stats --reference-fasta {args.hg38_reference_fasta} --verbose {filtered_comparison_catalog_path}", step_number=60)

        run(f"""python3 -u -m str_analysis.merge_loci \
            --output-prefix {catalog_name} \
            --output-format JSON \
            --overlapping-loci-action keep-first \
            --verbose \
            --write-merge-stats-tsv \
            {annotated_catalog_path} \
            {filtered_comparison_catalog_path}""", step_number=61)

    diff = time.time() - start_time
    print(f"Done with comparisons. Took {diff//3600:.0f}h, {(diff%3600)//60:.0f}m, {diff%60:.0f}s")

# Print final timing summary
end_datetime = datetime.datetime.now()
total_elapsed_time = time.time() - script_start_time
print("\n" + "="*80)
print(f"Started at:  {start_datetime.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Finished at: {end_datetime.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Total time:  {total_elapsed_time//3600:.0f}h, {(total_elapsed_time%3600)//60:.0f}m, {total_elapsed_time%60:.0f}s")
print("="*80)

