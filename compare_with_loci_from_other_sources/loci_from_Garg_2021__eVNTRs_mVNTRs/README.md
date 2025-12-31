# Garg 2021

**Title:** Pervasive cis effects of variation in copy number of large tandem repeats on local DNA methylation and gene expression

**Year:** 2021

**URL:** https://www.cell.com/ajhg/fulltext/S0002-9297(21)00100-2

**Type:** eVNTRs, mVNTRs

**Authors:** Garg P, Martin-Trujillo A, Rodriguez OL, Gies SJ, Halvardson J, Jazber A, Lauer E, Frey J, Huber N, Mazzarella G, Moreno LI, Pearlman A, Wolford BN, Muzny DM, Gibbs RA, Feuk L, Sharp AJ

**Status:** processed

**Total Loci Extracted:** 3,824 unique VNTRs

## Summary

Used read-depth data from Illumina WGS to perform association analysis between copy number of ~70,000 VNTRs (motif size >= 10 bp) with both gene expression (404 samples in 48 tissues) and DNA methylation (235 samples in peripheral blood). Identified thousands of VNTRs associated with local gene expression (eVNTRs) and DNA methylation levels (mVNTRs). Validated 73-80% of signals in independent cohort. Strong enrichments of eVNTRs and mVNTRs for regulatory features such as enhancers and promoters.

## Data Availability

Supplementary tables available at: https://www.cell.com/ajhg/fulltext/S0002-9297(21)00100-2#supplementaryMaterial

## Selection Criteria

Extracted unique VNTR loci from:
- Table S4 (mmc5.xlsx): 13,752 eVNTR associations → 2,949 unique loci
- Table S5 (mmc6.xlsx): 3,152 mVNTR associations → 875 additional unique loci

Filtering:
- Motif length <= 1000bp
- Valid hg38 coordinates

Note: All loci in this dataset have motif >= 10bp (VNTR definition per paper methods)

## Output

- `Garg_2021_eVNTRs_mVNTRs.bed.gz`: 3,824 unique VNTR loci (hg38)

