import hailtop.fs as hfs
import os
import re

from step_pipeline import pipeline, Backend, Localize

DOCKER_IMAGE = "us-central1-docker.pkg.dev/cmg-analysis/docker-repo/str-analysis-with-expansion-hunter@sha256:af2bece9abb119210e0f45a9e4963b2fee47c55ccb0c22b3855a678082c8c92a"

REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"

INPUT_CRAM_PATH = "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.cram"
INPUT_CRAI_PATH = "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.cram.crai"
VARIANT_CATALOG = "gs://tandem-repeat-catalog/v1.0/repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.json.gz"

SUBSET_CATALOGS = {
    "subset1": "gs://tandem-repeat-catalog/v1.0/subsets/TRExplorer_v1_subset1__excluding_1bp_and_2bp_motifs_and_intergenic.EH.json.gz",
    "subset2": "gs://tandem-repeat-catalog/v1.0/subsets/TRExplorer_v1_subset2__polymorphic_in_HPRC256.EH.json.gz",
    "subset3": "gs://tandem-repeat-catalog/v1.0/subsets/TRExplorer_v1_subset3__motifs_like_those_at_known_disease_associated_loci.EH.json.gz",
}

OUTPUT_BASE_DIR = "gs://bw2-delete-after-5-days/tool_results/expansion_hunter"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--subset", choices=SUBSET_CATALOGS.keys(), help="Use a subset catalog instead of the full catalog")
    parser.add_argument("--variant-catalog", default=VARIANT_CATALOG, help="Path of variant catalog json file to process")
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=INPUT_CRAM_PATH)
    parser.add_argument("--input-bai", default=INPUT_CRAI_PATH)
    parser.add_argument("--min-locus-coverage", type=int, help="Sets ExpansionHunter's --min-locus-coverage arg")
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    args = bp.parse_known_args()

    if args.subset:
        args.variant_catalog = SUBSET_CATALOGS[args.subset]

    catalog_name = os.path.basename(args.variant_catalog).replace(".json", "").replace(".gz", "")
    output_dir = os.path.join(args.output_dir, catalog_name)

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"Benchmarking: ExpansionHunter: {catalog_name}: {bam_path_ending}")

    create_expansion_hunter_step(
        bp,
        reference_fasta=args.reference_fasta,
        reference_fasta_fai=args.reference_fasta_fai,
        input_bam=args.input_bam,
        input_bai=args.input_bai,
        variant_catalog_path=args.variant_catalog,
        output_dir=output_dir,
        male_or_female="male",
        min_locus_coverage=args.min_locus_coverage)
    bp.run()


def create_expansion_hunter_step(bp, *, reference_fasta, reference_fasta_fai, input_bam, input_bai,
                                 variant_catalog_path, output_dir, male_or_female="male",
                                 min_locus_coverage=None):

    min_locus_coverage_arg = ""
    if min_locus_coverage is not None:
        min_locus_coverage_arg = f"--min-locus-coverage {min_locus_coverage}"

    hfs_ls_results = hfs.ls(input_bam)
    if len(hfs_ls_results) == 0:
        raise ValueError(f"No files found matching {input_bam}")
    input_bam_file_stats = hfs_ls_results[0]

    s1 = bp.new_step(
        f"Run ExpansionHunter on {os.path.basename(input_bam)} ({os.path.basename(variant_catalog_path)})",
        arg_suffix=f"run-expansion-hunter-step",
        step_number=1,
        image=DOCKER_IMAGE,
        cpu=2,
        memory="highmem",
        localize_by=Localize.GSUTIL_COPY,
        storage=f"{int(input_bam_file_stats.size/10**9) + 25}Gi",
        preemptible=False,
        output_dir=output_dir)

    local_fasta = s1.input(reference_fasta)
    if reference_fasta_fai:
        s1.input(reference_fasta_fai)

    local_bam = s1.input(input_bam)
    if input_bai:
        s1.input(input_bai)

    s1.command("set -ex")

    local_variant_catalog = s1.input(variant_catalog_path)

    output_prefix = re.sub(".json$", "", local_variant_catalog.filename)
    s1.command(f"echo Genotyping $(cat {local_variant_catalog} | grep LocusId | wc -l) loci")

    s1.command(f"""/usr/bin/time --verbose ExpansionHunter {min_locus_coverage_arg} \
        --reference {local_fasta} \
        --reads {local_bam} \
        --analysis-mode optimized-streaming \
        --sex {male_or_female} \
        --variant-catalog {local_variant_catalog} \
        --output-prefix {output_prefix}""")

    s1.command("ls -lhrt")

    s1.output(f"{output_prefix}.json", output_dir=os.path.join(output_dir, f"json"))

    return s1


if __name__ == "__main__":
    main()
