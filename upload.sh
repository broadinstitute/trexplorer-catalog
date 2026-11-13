set -ex

gcloud storage cp results__2026-10-29/release_draft_2026-10-29/TRExplorer.*.gz*   gs://tandem-repeat-catalog/v2.1/

# The annotated TSV is written to the build folder rather than the release folder
gcloud storage cp results__2026-10-29/1_to_1000bp_motifs/TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.EH.with_annotations.tsv.gz   gs://tandem-repeat-catalog/v2.1/

# Each build names its log after its start time, so publish it under one fixed name that replaces the previous log.
# This build ran steps 1-11, then the rest (--start-with-step 12), then on 11/10/2026 steps 16 onward again to add the
# PseudoautosomalRegion and SegDup fields (--start-with-step 16), and then step 61 again (--only-step 61), so its four
# logs are joined.
cat results__2026-10-29/release_draft_2026-10-29/TRExplorer.repeat_catalog_v2.1.hg38.all_steps.20261029-233040.log \
    results__2026-10-29/release_draft_2026-10-29/TRExplorer.repeat_catalog_v2.1.hg38.all_steps.20261030-100102.log \
    results__2026-10-29/release_draft_2026-10-29/TRExplorer.repeat_catalog_v2.1.hg38.all_steps.20261110-183246.log \
    results__2026-10-29/release_draft_2026-10-29/TRExplorer.repeat_catalog_v2.1.hg38.all_steps.20261110-221125.log \
    | gcloud storage cp -   gs://tandem-repeat-catalog/v2.1/TRExplorer.repeat_catalog_v2.1.hg38.all_steps.log

# The extended locus boundaries source catalog is generated separately by
# extended_loci/generate_extended_loci_catalog.sh and is the default input for the pipeline's
# --extended-loci-catalog argument, so it is published alongside the release it belongs to.
gcloud storage cp extended_loci/TRExplorer.extended_loci_v2.1.hg38.bed.gz*   gs://tandem-repeat-catalog/v2.1/

# The CHM13 catalog (run_all_steps_to_generate_the_CHM13_catalog.py). Its file names include "chm13", so it shares the
# v2.1 folder with the hg38 catalog without overwriting anything.
gcloud storage cp results__chm13__2026-10-30/release_draft_2026-10-30/TRExplorer.*.gz*   gs://tandem-repeat-catalog/v2.1/
gcloud storage cp results__chm13__2026-10-30/TRExplorer.repeat_catalog_v2.1.chm13.1_to_1000bp_motifs.EH.with_annotations.tsv.gz   gs://tandem-repeat-catalog/v2.1/
# This build ran all steps, then on 11/10/2026 steps 12 onward again to add the PseudoautosomalRegion and SegDup fields
# (--start-with-step 12), so its two logs are joined.
cat results__chm13__2026-10-30/release_draft_2026-10-30/TRExplorer.repeat_catalog_v2.1.chm13.all_steps.20261030-142833.log \
    results__chm13__2026-10-30/release_draft_2026-10-30/TRExplorer.repeat_catalog_v2.1.chm13.all_steps.20261110-204705.log \
    | gcloud storage cp -   gs://tandem-repeat-catalog/v2.1/TRExplorer.repeat_catalog_v2.1.chm13.all_steps.log

# The hg19 catalog (run_all_steps_to_generate_the_hg19_catalog.py). Its file names include "hg19", so it shares the
# v2.1 folder with the hg38 and CHM13 catalogs without overwriting anything.
gcloud storage cp results__hg19__2026-11-12/release_draft_2026-11-12/TRExplorer.*.gz*   gs://tandem-repeat-catalog/v2.1/
gcloud storage cp results__hg19__2026-11-12/TRExplorer.repeat_catalog_v2.1.hg19.1_to_1000bp_motifs.EH.with_annotations.tsv.gz   gs://tandem-repeat-catalog/v2.1/
# This build ran steps 1-5, where step 5 stopped on a locus that lifts onto a contig the hg19 reference doesn't have,
# then step 5 again once lift_over_catalog.py dropped such loci, then steps 6-26 once the perfect repeats in hg19 were
# ready, and then steps 27-34, so its four logs are joined.
cat results__hg19__2026-11-12/release_draft_2026-11-12/TRExplorer.repeat_catalog_v2.1.hg19.all_steps.20261112-190359.log \
    results__hg19__2026-11-12/release_draft_2026-11-12/TRExplorer.repeat_catalog_v2.1.hg19.all_steps.20261112-190854.log \
    results__hg19__2026-11-12/release_draft_2026-11-12/TRExplorer.repeat_catalog_v2.1.hg19.all_steps.20261113-010313.log \
    results__hg19__2026-11-12/release_draft_2026-11-12/TRExplorer.repeat_catalog_v2.1.hg19.all_steps.20261113-024423.log \
    | gcloud storage cp -   gs://tandem-repeat-catalog/v2.1/TRExplorer.repeat_catalog_v2.1.hg19.all_steps.log
