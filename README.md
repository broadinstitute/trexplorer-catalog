

### Tandem Repeat Catalog & Variation Clusters

This repo contains a genome-wide TR catalog that stratifies TRs into 2 groups:
1) isolated TRs suitable for traditional repeat copy number analysis using short-read or long-read data
2) TRs embedded within wider polymorphic regions (ie. variation clusters) that are best studied through sequence-level analysis

TRExplorer catalog [release v2.0](https://github.com/broadinstitute/tandem-repeat-catalogs/releases/tag/v2.0) is now available. 

The [http://trexplorer.broadinstitute.org](http://trexplorer.broadinstitute.org) companion portal supports interactive exploration, visualization, filtering, and export of the catalog.

See the `Using the TRExplorer Catalog` section of the manuscript below for practical considerations on how to use these resources.

----

### Publication

<b>Defining a tandem repeat catalog and variation clusters for genome-wide analyses and population databases</b><br />
Ben Weisburd, Egor Dolzhenko, Mark F. Bennett, Matt C. Danzi, Adam English, Laurel Hiatt, Hope Tanudisastro, Nehir Edibe Kurtas, Helyaneh Ziaei Jam, Harrison Brand, Fritz J. Sedlazeck, Melissa Gymrek, Harriet Dashnow, Michael A. Eberle, Heidi L. Rehm

[bioRxiv 2024.10.04.615514v5](https://www.biorxiv.org/content/10.1101/2024.10.04.615514v5)

----

### Background

Tandem repeats (TRs) are regions of the genome that consist of consecutive copies of some motif sequence. For example, `CAGCAGCAG` is a tandem repeat of the `CAG` motif. Many types of genomic studies require annotations of tandem repeats in the reference genome, called repeat catalogs, which specify the genomic start and end coordinates of each tandem repeat region, as well as the one or more motifs that repeat there. 

For example, if a hypothetical region at the beginning of `chrX` consisted of the following nucleotide sequence:  
`ATCAGTAGA ATATATATAT CAGACAGCAGCAG TGAGTGCGTAC...`  
it could be represented in a repeat catalog as two entries:  
`chrX:10-19 (AT)*`  
`chrX:20-32 (CAG)*`   
indicating that a repeat of the `AT` motif occurs between positions 10 and 19 (inclusive), and of the `CAG` motif between positions 20 and 32.
A genome-wide catalog would contain such entries for all repeat regions of interest found anywhere in the genome. 


----

### v2.0 Catalog Construction

See the [release notes](https://github.com/broadinstitute/tandem-repeat-catalogs/releases/tag/v2.0) 

### v1.0 Catalog Construction

The genome-wide TR catalog was created by combining 4 source catalogs in order:
1) [Known disease-associated loci](https://github.com/broadinstitute/str-analysis/blob/69dd90ecbc1dcbb23d5ca84ab4022850a283114f/str_analysis/variant_catalogs/variant_catalog_with_offtargets.GRCh37.json)
2) [Illumina catalog of 174k polymorphic repeats](https://github.com/Illumina/RepeatCatalogs?tab=readme-ov-file) 
3) All perfect repeats in hg38 that span ≥ 9bp and consist of at least 3 repeats of any motif between 1 and 1000 bp in size. These were identified using [ColabRepeatFinder](https://github.com/broadinstitute/colab-repeat-finder).
4) Catalog of polymorphic loci computed by applying methods described in [[Weisburd et al. 2023](https://www.biorxiv.org/content/10.1101/2023.05.05.539588v1)] to 78 haplotype-resolved T2T assemblies from the HPRC and HGSVC

The numbers (and %) of loci in the combined catalog that were added from each of the source catalogs were as follows:

```
           83 out of        83 (100.0%) TRs from source 1: known disease-associated loci as well as 20 adjacent or historical candidate loci 
      174,244 out of   174,286 (100.0%) TRs from source 2: Illumina catalog of 174k polymorphic loci
    4,391,197 out of 4,558,281 ( 96.3%) TRs from source 3: perfect repeats in hg38
      297,517 out of 1,937,805 ( 15.4%) TRs from source 4: polymorphic loci in 78 haplotype-resolved T2T assemblies
```
The merging procedure involved taking all loci from the 1st catalog, then all non-duplicate loci from the next catalog, then from the third catalog and so on, in the order listed above. A locus was considered a duplicate if it overlapped a previously-added locus by at least 66% and the two loci had the same motif after cyclic shift and/or reverse complement (ie. CAG, AGC, GCA, CTG, TGC, GCT were considered to be the same motif). 

### v1.0 Catalog Stats

The following catalog stats for v1.0 were computed using [str_analysis/compute_catalog_stats.py](https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/compute_catalog_stats.py):

```
Stats for repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz:
    4,863,041 total loci
   65,678,112 base pairs spanned by all loci (2.127% of the genome)
    3,210,115 out of  4,863,041 ( 66.0%) repeat interval size is an integer multiple of the motif size (aka. trimmed)
    1,567,337 out of  4,863,041 ( 32.2%) repeat intervals are homopolymers
       18,340 out of  4,863,041 (  0.4%) repeat intervals overlap each other by at least two motif lengths
           11 out of  4,863,041 (  0.0%) repeat intervals have non-ACGT motifs
Examples of overlapping repeats: chr1:82008141-82008152, chr3:78937990-78938032, chr4:1046750-1046794, chr5:52437153-52437201, chr6:34683425-34683453, chr4:107646310-107646327, chr18:52413438-52413450, chr6:150149299-150149323, chr7:40295582-40295597, chr9:35561915-35561931

Ranges:
   Motif size range:  1-833bp
   Locus size range:  1-2523bp
   Num repeats range: 1-300x repeats

   Max locus size =   2,523bp               @ chrX:71520430-71522953 (CCAGCACTTTGGGAGGCCGAGGCAGGCTGATCACTAGGTCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCCCGTCTCTACTAAAAATACAAAAATTACCTGGGTGTGGGGGTGGGCACCTGTAATCCCAGCTACTCGGGAGGCTGGGGAGGCAGGAGAATTGCCTGAACCTGAGAGGCAGAGGCTGCAGTGAGCTGAGATTGTGCCACTGCACTCCAGCCTGGGCGACAGAGTGAGACTCAGTCTCAAAACAAAAAAAAAAAAAGATTTTAGTAACTTTTATCCTGTTTTAATAATACTGACTCAGAAACTATAATGTGTACTTTATAATTTACTTCCTAGATGACACTTGATTTTCTTCAAGAGCAAGATAGCTGCCCTGTGCAGTTGGTCTCCTTGAAAACTATTTTAGTTCTATCATAATTTCCTGTGATAAATATTTTGACCTTCTAAAATTTCAGAATATTGCACCAAGTAGAAAGAAAATAGGTTTTTTCTCTTTTCTTCTTCTTCCTTTTTTTTTTCTGAGAAAGAGGGAATGAGAACTTTAGTGTTCTTTCAATAGCGTTCTTATTTGTAGAAATGCATAATAGTGTCCTAGTAAGGCTTGACAATAACTCTGGTCTTCATCATATTTTGTGATAAAACTTTTGATTTAAAAAAACCTCTGATCTATTTATCATGGCAAATGGATAGAGCTTTCCTGCCTGTTTTCTTTCTTTTCTTTTTTCTTTCTTTCCTTTTTTTTCCTTTGAGCTTAGATTTTTAGAAGCACATATTTAAAAATCAGGTATAAGACTGGATGCAGTGGCTCACGCCTGTAATC)
   Min reference repeat purity   =  0.43    @ chr3:112804380-112804514 (TCT)
   Min overall mappability       =  0.00    @ chrY:56887882-56887891 (TGA)
   Base-level   purity   median: 1.000,  mean: 0.999

          chrX:    244,191 out of  4,863,041 (  5.0%) repeat intervals
          chrY:     39,257 out of  4,863,041 (  0.8%) repeat intervals
          chrM:         14 out of  4,863,041 (  0.0%) repeat intervals
   alt contigs:          0 out of  4,863,041 (  0.0%) repeat intervals

Motif size distribution:
          1bp:  1,567,337 out of  4,863,041 ( 32.2%) repeat intervals
          2bp:    978,972 out of  4,863,041 ( 20.1%) repeat intervals
          3bp:  1,432,117 out of  4,863,041 ( 29.4%) repeat intervals
          4bp:    590,787 out of  4,863,041 ( 12.1%) repeat intervals
          5bp:    177,422 out of  4,863,041 (  3.6%) repeat intervals
          6bp:     56,731 out of  4,863,041 (  1.2%) repeat intervals
       7-24bp:     43,996 out of  4,863,041 (  0.9%) repeat intervals
        25+bp:     15,679 out of  4,863,041 (  0.3%) repeat intervals

Num repeats in reference:
           1x:     10,443 out of  4,863,041 (  0.2%) repeat intervals
           2x:     38,922 out of  4,863,041 (  0.8%) repeat intervals
           3x:  1,799,189 out of  4,863,041 ( 37.0%) repeat intervals
           4x:    650,397 out of  4,863,041 ( 13.4%) repeat intervals
           5x:    356,525 out of  4,863,041 (  7.3%) repeat intervals
           6x:    151,893 out of  4,863,041 (  3.1%) repeat intervals
           7x:     85,760 out of  4,863,041 (  1.8%) repeat intervals
           8x:    257,475 out of  4,863,041 (  5.3%) repeat intervals
           9x:    352,993 out of  4,863,041 (  7.3%) repeat intervals
       10-15x:    759,188 out of  4,863,041 ( 15.6%) repeat intervals
       16-25x:    348,837 out of  4,863,041 (  7.2%) repeat intervals
       26-35x:     45,478 out of  4,863,041 (  0.9%) repeat intervals
       36-50x:      5,610 out of  4,863,041 (  0.1%) repeat intervals
         51+x:        331 out of  4,863,041 (  0.0%) repeat intervals

Reference repeat purity distribution:
          0.0:          0 out of  4,863,041 (  0.0%) repeat intervals
          0.1:          0 out of  4,863,041 (  0.0%) repeat intervals
          0.2:          0 out of  4,863,041 (  0.0%) repeat intervals
          0.3:          0 out of  4,863,041 (  0.0%) repeat intervals
          0.4:          3 out of  4,863,041 (  0.0%) repeat intervals
          0.5:         14 out of  4,863,041 (  0.0%) repeat intervals
          0.6:         44 out of  4,863,041 (  0.0%) repeat intervals
          0.7:      1,570 out of  4,863,041 (  0.0%) repeat intervals
          0.8:     12,336 out of  4,863,041 (  0.3%) repeat intervals
          0.9:     21,126 out of  4,863,041 (  0.4%) repeat intervals
          1.0:  4,827,948 out of  4,863,041 ( 99.3%) repeat intervals

Mappability distribution:
          0.0:    154,279 out of  4,863,041 (  3.2%) loci
          0.1:    214,471 out of  4,863,041 (  4.4%) loci
          0.2:    246,877 out of  4,863,041 (  5.1%) loci
          0.3:    236,856 out of  4,863,041 (  4.9%) loci
          0.4:    391,388 out of  4,863,041 (  8.0%) loci
          0.5:    561,639 out of  4,863,041 ( 11.5%) loci
          0.6:    352,273 out of  4,863,041 (  7.2%) loci
          0.7:    306,208 out of  4,863,041 (  6.3%) loci
          0.8:    337,715 out of  4,863,041 (  6.9%) loci
          0.9:    626,047 out of  4,863,041 ( 12.9%) loci
          1.0:  1,435,288 out of  4,863,041 ( 29.5%) loci

Locus sizes at each motif size:
     1bp motifs: locus size range:      1 bp to      90 bp  (median:   11 bp) based on  1,567,337 loci. Mean base purity: 1.00.  Mean mappability: 0.66
     2bp motifs: locus size range:      2 bp to     600 bp  (median:   10 bp) based on    978,972 loci. Mean base purity: 1.00.  Mean mappability: 0.76
     3bp motifs: locus size range:      3 bp to     632 bp  (median:    9 bp) based on  1,432,117 loci. Mean base purity: 1.00.  Mean mappability: 0.75
     4bp motifs: locus size range:      4 bp to     533 bp  (median:   14 bp) based on    590,787 loci. Mean base purity: 1.00.  Mean mappability: 0.68
     5bp motifs: locus size range:      5 bp to     400 bp  (median:   18 bp) based on    177,422 loci. Mean base purity: 1.00.  Mean mappability: 0.61
     6bp motifs: locus size range:      6 bp to   1,103 bp  (median:   20 bp) based on     56,731 loci. Mean base purity: 1.00.  Mean mappability: 0.62
     7bp motifs: locus size range:      7 bp to     151 bp  (median:   22 bp) based on     15,083 loci. Mean base purity: 1.00.  Mean mappability: 0.58
     8bp motifs: locus size range:      8 bp to     312 bp  (median:   25 bp) based on      7,107 loci. Mean base purity: 1.00.  Mean mappability: 0.57
     9bp motifs: locus size range:      9 bp to     153 bp  (median:   28 bp) based on      3,231 loci. Mean base purity: 1.00.  Mean mappability: 0.51
    10bp motifs: locus size range:     10 bp to     150 bp  (median:   31 bp) based on      2,713 loci. Mean base purity: 1.00.  Mean mappability: 0.50
```

Additional stats can be found in the [[run log](https://raw.githubusercontent.com/broadinstitute/tandem-repeat-catalogs/main/all_steps.log)]



