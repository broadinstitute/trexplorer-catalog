import hailtop.fs as hfs
import os
import re

from step_pipeline import pipeline, Backend, Localize

DOCKER_IMAGE = "us-central1-docker.pkg.dev/cmg-analysis/docker-repo/str-analysis-with-trgt@sha256:1b3ee5020d94baa000b46d28871d154f7bfa4ea9ca174ba88431c0bcd4c65cad"

REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"

INPUT_BAM_PATH = "gs://str-truth-set-v2/raw_data/HG002/pacbio/HG002.downsampled_to_30x.bam"
INPUT_BAI_PATH = "gs://str-truth-set-v2/raw_data/HG002/pacbio/HG002.downsampled_to_30x.bam.bai"

TRGT_CATALOG_BED = "gs://tandem-repeat-catalog/v1.0/variation_clusters_and_isolated_TRs_v1.hg38.TRGT.bed.gz"

SUBSET_CATALOGS = {
    "subset1": "gs://tandem-repeat-catalog/v1.0/subsets/TRExplorer_v1_subset1__excluding_1bp_and_2bp_motifs_and_intergenic.TRGT.bed.gz",
    "subset2": "gs://tandem-repeat-catalog/v1.0/subsets/TRExplorer_v1_subset2__polymorphic_in_HPRC256.TRGT.bed.gz",
    "subset3": "gs://tandem-repeat-catalog/v1.0/subsets/TRExplorer_v1_subset3__motifs_like_those_at_known_disease_associated_loci.TRGT.bed.gz",
}

OUTPUT_BASE_DIR = "gs://bw2-delete-after-5-days/tool_results/trgt"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--subset", choices=SUBSET_CATALOGS.keys(), help="Use a subset catalog instead of the full catalog")
    parser.add_argument("--trgt-catalog-bed", default=TRGT_CATALOG_BED, help="Path of TRGT catalog bed file to process")
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=INPUT_BAM_PATH)
    parser.add_argument("--input-bai", default=INPUT_BAI_PATH)
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    args = bp.parse_known_args()

    if args.subset:
        args.trgt_catalog_bed = SUBSET_CATALOGS[args.subset]

    catalog_name = os.path.basename(args.trgt_catalog_bed).replace(".bed", "").replace(".gz", "")
    trgt_catalog_bed_paths = [x.path for x in hfs.ls(args.trgt_catalog_bed)]
    if len(trgt_catalog_bed_paths) == 0:
        raise ValueError(f"No files found matching {args.trgt_catalog_bed}")

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"Benchmarking: TRGT: {catalog_name}: {bam_path_ending}")
    output_dir = os.path.join(args.output_dir, catalog_name)

    create_trgt_step(
        bp,
        reference_fasta=args.reference_fasta,
        reference_fasta_fai=args.reference_fasta_fai,
        input_bam=args.input_bam,
        input_bai=args.input_bai,
        trgt_catalog_bed_path=trgt_catalog_bed_paths[0],
        output_dir=output_dir)
    bp.run()


def create_trgt_step(bp, *, reference_fasta, reference_fasta_fai, input_bam, input_bai, trgt_catalog_bed_path,
                     output_dir, male_or_female="male"):
    output_prefix = re.sub(".bed(.gz)?$", "", os.path.basename(trgt_catalog_bed_path))

    hfs_ls_results = hfs.ls(input_bam)
    if len(hfs_ls_results) == 0:
        raise ValueError(f"No files found matching {input_bam}")
    input_bam_file_stats = hfs_ls_results[0]

    s1 = bp.new_step(f"Run TRGT on {os.path.basename(input_bam)}  {os.path.basename(trgt_catalog_bed_path)}",
                     arg_suffix=f"run-trgt-step",
                     step_number=1,
                     image=DOCKER_IMAGE,
                     cpu=1,
                     localize_by=Localize.GSUTIL_COPY,
                     storage=f"{int(input_bam_file_stats.size/10**9) + 25}Gi",
                     preemptible=False,
                     output_dir=output_dir)
    s1.command("set -ex")
    local_fasta = s1.input(reference_fasta)
    if reference_fasta_fai:
        s1.input(reference_fasta_fai)
    else:
        s1.input(f"{reference_fasta}.fai")
    local_bam = s1.input(input_bam)
    if input_bai:
        s1.input(input_bai)
    local_trgt_catalog_bed = s1.input(trgt_catalog_bed_path)
    s1.command("df -kh")
    karyotype = "XX" if male_or_female == "female" else "XY"
    s1.command(f"echo Genotyping $(cat {local_trgt_catalog_bed} | wc -l) loci in {local_bam.filename}  karyotype={karyotype}")
    s1.command(f"""/usr/bin/time --verbose trgt genotype \
                                     --genome {local_fasta} \
                                     --reads {local_bam} \
                                     --karyotype {karyotype} \
                                     --repeats {local_trgt_catalog_bed} \
                                     --output-prefix {output_prefix}
    """)
    s1.command("ls -lhrt")

    s1.output(f"{output_prefix}.vcf.gz")
    s1.output(f"{output_prefix}.spanning.bam")

    return s1

if __name__ == "__main__":
    main()
