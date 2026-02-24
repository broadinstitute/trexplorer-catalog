import hailtop.fs as hfs
import os
import re

from step_pipeline import pipeline, Backend, Localize

DOCKER_IMAGE = "us-central1-docker.pkg.dev/cmg-analysis/docker-repo/str-analysis-with-longtr@sha256:0aa8f9eac56b4b52ce0ee190b86f8d1047a1a8b6f99c7057483638e0d8b63617"

REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"

INPUT_BAM_PATH = "gs://str-truth-set-v2/raw_data/HG002/pacbio/HG002.downsampled_to_30x.bam"
INPUT_BAI_PATH = "gs://str-truth-set-v2/raw_data/HG002/pacbio/HG002.downsampled_to_30x.bam.bai"

REGIONS_BED = "gs://tandem-repeat-catalog/v1.0/variation_clusters_and_isolated_TRs_v1.hg38.LongTR.bed.gz"

SUBSET_CATALOGS = {
    "subset1": "gs://tandem-repeat-catalog/v1.0/subsets/TRExplorer_v1_subset1__excluding_1bp_and_2bp_motifs_and_intergenic.LongTR.bed.gz",
    "subset2": "gs://tandem-repeat-catalog/v1.0/subsets/TRExplorer_v1_subset2__polymorphic_in_HPRC256.LongTR.bed.gz",
    "subset3": "gs://tandem-repeat-catalog/v1.0/subsets/TRExplorer_v1_subset3__motifs_like_those_at_known_disease_associated_loci.LongTR.bed.gz",
}

OUTPUT_BASE_DIR = "gs://bw2-delete-after-5-days/tool_results/longtr"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--subset", choices=SUBSET_CATALOGS.keys(), help="Use a subset catalog instead of the full catalog")
    parser.add_argument("--regions-bed", default=REGIONS_BED, help="Path of regions bed file to process")
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=INPUT_BAM_PATH)
    parser.add_argument("--input-bai", default=INPUT_BAI_PATH)
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    args = bp.parse_known_args()

    if args.subset:
        args.regions_bed = SUBSET_CATALOGS[args.subset]

    catalog_name = os.path.basename(args.regions_bed).replace(".bed", "").replace(".gz", "")

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"Benchmarking: LongTR: {catalog_name}: {bam_path_ending}")
    output_dir = os.path.join(args.output_dir, catalog_name)

    create_longtr_step(
        bp,
        reference_fasta=args.reference_fasta,
        reference_fasta_fai=args.reference_fasta_fai,
        input_bam=args.input_bam,
        input_bai=args.input_bai,
        regions_bed_path=args.regions_bed,
        output_dir=output_dir)
    bp.run()


def create_longtr_step(bp, *, reference_fasta, reference_fasta_fai, input_bam, input_bai, regions_bed_path, output_dir):
    hfs_ls_results = hfs.ls(input_bam)
    if len(hfs_ls_results) == 0:
        raise ValueError(f"No files found matching {input_bam}")
    input_bam_file_stats = hfs_ls_results[0]

    input_bam_filename_prefix = re.sub("(.bam|.cram)$", "", os.path.basename(input_bam))

    s1 = bp.new_step(f"Run LongTR on {os.path.basename(input_bam)} ({os.path.basename(regions_bed_path)})",
                     arg_suffix=f"run-longtr-step",
                     step_number=1,
                     image=DOCKER_IMAGE,
                     cpu=1,
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

    local_regions_bed = s1.input(regions_bed_path)
    output_prefix = re.sub(".bed(.gz)?$", "", local_regions_bed.filename)
    decompressed_bed = re.sub(".gz$", "", local_regions_bed.filename)
    if regions_bed_path.endswith(".gz"):
        s1.command(f"gunzip -c {local_regions_bed} > {decompressed_bed}")
    else:
        decompressed_bed = str(local_regions_bed)
    s1.command(f"echo Genotyping $(cat {decompressed_bed} | wc -l) loci")
    s1.command("set -ex")
    s1.command(f"""/usr/bin/time --verbose LongTR \
            --min-reads 2 \
            --bam-samps {input_bam_filename_prefix} \
            --bam-libs {input_bam_filename_prefix} \
            --bams {local_bam} \
            --fasta {local_fasta} \
            --regions {decompressed_bed} \
            --log {output_prefix}.log \
            --tr-vcf {output_prefix}.vcf.gz""")

    s1.command("ls -lhrt")

    s1.output(f"{output_prefix}.vcf.gz", output_dir=os.path.join(output_dir, f"vcf"))
    s1.output(f"{output_prefix}.log", output_dir=os.path.join(output_dir, f"log"))

    return s1

if __name__ == "__main__":
    main()
