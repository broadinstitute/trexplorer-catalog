# Mukamel 2021 (protein-coding VNTRs)

**Title:** Protein-coding repeat polymorphisms strongly shape diverse human phenotypes

**Year:** 2021

**URL:** https://www.science.org/doi/10.1126/science.abg8289

**Type:** coding VNTRs

**Status:** this source is already in the v2 catalog, but via `../loci_from_Mukamel_2021`, not via this folder

**Total Loci in Paper:** 118 (table S1)

**Loci Extracted Here:** 105

## What this folder is

`generate_bed.py` re-extracts table S1 from the Excel supplement
(`../loci_from_Mukamel_2021/science.abg8289_tables_s1_to_s8/science.abg8289_Tables_S1_to_S8.xlsx`,
sheet "ST1 Protein-coding VNTRs") and writes the paper's raw `Repeat unit` motifs to
`Mukamel_2021_coding_VNTRs.bed.gz`.

Those 105 loci are identical (same coordinates, same motifs) to the pre-existing
`../loci_from_Mukamel_2021/VNTRs_in_ST1.bed`, which had no generating script. So this folder documents
the provenance of that file rather than producing a new source catalog.

## Selection Criteria

118 rows in table S1, minus 13:

- 10 whose Notes column says "Period length probably incorrectly estimated."
- 3 with repeat units > 2,000 bp (10,553 bp, 5,544 bp, 4,059 bp)

Leaves 105 loci. No filtering against the v2 draft catalog is done here.

## Relationship to what v2 actually includes

`run_all_steps_to_generate_the_catalog.py` (source `TRExplorerV2:Mukamel2021`, merge tier v2b)
downloads `gs://tandem-repeat-catalog/v2.0/Mukamel_2021.loci_to_include_in_catalog.bed.gz`, which is
produced by `../loci_from_Mukamel_2021/select_loci_to_include_in_catalog.py`. That script starts from
the same 105 loci but runs `select_loci(..., adjust_motifs_to_maximize_purity=True,
min_repeats_in_reference=2, min_adjusted_motif_purity=0.2, keep_only_motifs_with_ACGT_bases=True)` and
writes the adjusted motifs.

Result: same 105 coordinates, but 104 of the 105 motifs differ from the bed in this folder. 34 are pure
rotations, 66 have base substitutions, and 4 collapse a 9 bp motif to its 3 bp period (for example
chr1:179534854-179534923 is `AAAGAAGAA` here and `AAG` in the shipped bed).

All 105 loci are passed to the pipeline; deduplication against the earlier source catalogs happens
inside the v2b merge. For reference, scoring the 105 against the v2 draft
(`../loci_from_Mukamel_2021/vntrs_in_ST1.overlap_with_TRExplorer_v2.tsv.gz`) gives 51 with "diff <= 2
repeats" and 13 with Jaccard > 0.66, leaving 41 that are new; those 41 are in
`../loci_from_Mukamel_2021/vntrs_in_ST1.loci_to_include_in_catalog.bed.gz`.

## Output

- `Mukamel_2021_coding_VNTRs.bed.gz`: 105 loci (hg38) with the raw table S1 motifs. Not uploaded to
  `gs://tandem-repeat-catalog/` and not read by any pipeline step.
