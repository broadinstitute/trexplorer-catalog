set -ex

gcloud storage cp results__2026-10-19/release_draft_2026-10-19/TRExplorer.*.gz*   gs://tandem-repeat-catalog/v2.1/

# The annotated TSV is written to the build folder rather than the release folder
gcloud storage cp results__2026-10-19/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.EH.with_annotations.tsv.gz   gs://tandem-repeat-catalog/v2.1/

# Each build names its log after its start time, so publish it under one fixed name that replaces the previous log
gcloud storage cp results__2026-10-19/release_draft_2026-10-19/TRExplorer.repeat_catalog_v2.1.hg38.all_steps.20261019-022932.log   gs://tandem-repeat-catalog/v2.1/TRExplorer.repeat_catalog_v2.1.hg38.all_steps.log

# The extended locus boundaries source catalog is generated separately by
# extended_loci/generate_extended_loci_catalog.sh and is the default input for the pipeline's
# --extended-loci-catalog argument, so it is published alongside the release it belongs to.
gcloud storage cp extended_loci/TRExplorer.extended_loci_v2.1.hg38.bed.gz*   gs://tandem-repeat-catalog/v2.1/
