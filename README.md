
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EnrichKit

<!-- badges: start -->

<!-- badges: end -->

The goal of EnrichKit is to **perform over-representation** test of a
given gene set pair (*SignificantGenes* and *TotolGenes*) based on
hypergeometric distribution (**Fisher’s exact test**). Gene sets could
possily be non-preserved co-expression modules, differentially expressed
gene (DEG), genes flagged by significant SNPs and etc..

Currently, **six pathway/annotation databases** are integrated:

  - [Gene
    Ontology](http://ensemblgenomes.org/info/access/biomart)
  - [KEGG](https://www.genome.jp/kegg/)
  - [Interpro](http://ensemblgenomes.org/info/access/biomart)
  - [MeSH](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C10&q=The+MeSH+translation+maintenance+system%3A+structure%2C+interface+design%2C+and+implementation.&btnG=)
  - [Reactome](https://reactome.org/download-data)
  - [Molecular
    Signatures](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/)

Also, gene identifiers could be *Ensembl Gene ID*, *EntrezID* or *HGNC
Gene Symbol*.

Note current release can only support **Bos Taurus**, other organism
might be included in future release.

## Example

Soppose we have identified 2 out 5 DEG in each of the two lactations.

``` r
require(EnrichKit,quietly = TRUE)
#> 
#> Warning: replacing previous import 'data.table::last' by 'dplyr::last' when
#> loading 'EnrichKit'
#> Warning: replacing previous import 'data.table::first' by 'dplyr::first' when
#> loading 'EnrichKit'
#> Warning: replacing previous import 'biomaRt::select' by 'dplyr::select' when
#> loading 'EnrichKit'
#> Warning: replacing previous import 'data.table::between' by 'dplyr::between'
#> when loading 'EnrichKit'
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
#> Warning: replacing previous import 'WGCNA::cor' by 'stats::cor' when loading
#> 'EnrichKit'
#> Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when
#> loading 'EnrichKit'
#> Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
#> 'EnrichKit'
Sig_lac1 =   c("ENSBTAG00000012594","ENSBTAG00000049850")
Sig_lac2 =   c("ENSBTAG00000009188","ENSBTAG00000001258")
Tot_lac1 = c("ENSBTAG00000012594","ENSBTAG00000049850","ENSBTAG00000018278","ENSBTAG00000021997","ENSBTAG00000008482")
Tot_lac2 = c("ENSBTAG00000009188","ENSBTAG00000001258","ENSBTAG00000021819","ENSBTAG00000019404","ENSBTAG00000015212")

GeneInfo = convertNformatID(GeneSetNames=c("lactation1","lactation2"),
                            SigGene_list = list(Sig_lac1,Sig_lac2),
                            TotalGene_list = list(Tot_lac1,Tot_lac2),
                            IDtype = "ens")
#> Maps last updated on: Thu Oct 24 12:31:05 2019
#> Warning in checkGeneSymbols(notsure): x contains non-approved gene symbols
#> Maps last updated on: Thu Oct 24 12:31:05 2019
GeneInfo
#> $lactation1
#>                 Gene  ENTREZID       SYMBOL SYMBOL_Suggested Sig
#> 1 ENSBTAG00000012594    615431        MRPS6            MRPS6   1
#> 2 ENSBTAG00000049850 112448353 LOC112448353     LOC112448353   1
#> 3 ENSBTAG00000018278    281640       ATP5PO           ATP5PO   0
#> 4 ENSBTAG00000021997    510879        ITSN1            ITSN1   0
#> 5 ENSBTAG00000008482    516462          SON              SON   0
#> 
#> $lactation2
#>                 Gene ENTREZID  SYMBOL SYMBOL_Suggested Sig
#> 1 ENSBTAG00000009188   281183    GART             GART   1
#> 2 ENSBTAG00000001258   615627 TMEM50B          TMEM50B   1
#> 3 ENSBTAG00000021819   282257  IFNAR1           IFNAR1   0
#> 4 ENSBTAG00000019404   767864  IL10RB           IL10RB   0
#> 5 ENSBTAG00000015212   282258  IFNAR2           IFNAR2   0
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!

## Installation

You can install the released version of EnrichKit from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("EnrichKit")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("liulihe954/EnrichKit")
```
