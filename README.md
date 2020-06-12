
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EnrichKit

Authors: Lihe Liu and Francisco Peñagaricano  
Maintainer: Lihe Liu (<lihe.liu@ufl.edu>)

The goal of EnrichKit is to **perform over-representation** test of
biological pathways (gene sets) within a given gene list pair
(*SignificantGenes* and *TotolGenes*) based on hypergeometric
distribution (**Fisher’s exact test**). Gene lists of interests could
possily be non-preserved co-expression modules, differentially expressed
gene (DEG), genes flagged by significant SNPs and
etc..

<div class="figure">

<img src="man/figures/README-Enrich_Illustration.png" alt=" " width="400%" />

<p class="caption">

</p>

</div>

Currently, **six pathway/annotation databases** are integrated in
current release:

  - [Gene
    Ontology](http://ensemblgenomes.org/info/access/biomart)
  - [KEGG](https://www.genome.jp/kegg/)
  - [Interpro](http://ensemblgenomes.org/info/access/biomart)
  - [MeSH](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C10&q=The+MeSH+translation+maintenance+system%3A+structure%2C+interface+design%2C+and+implementation.&btnG=)
  - [Reactome](https://reactome.org/download-data)
  - [Molecular
    Signatures](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/)

Note current release can only support **Bos Taurus**, other organism
might be included in future release.

Latest update 06-08-2020.

## Installation

Currently not published on [CRAN](https://CRAN.R-project.org).

Users could use the development version from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("liulihe954/EnrichKit") # Depends on R (>= 3.5.0)
```

## Example

Soppose we have identified 2 out 5 DEG in each of the two lactations.

``` r
library(EnrichKit)
# input format
Sig_lac1 =   c("ENSBTAG00000012594","ENSBTAG00000004139")
Sig_lac2 =   c("ENSBTAG00000009188","ENSBTAG00000001258")
Tot_lac1 = c("ENSBTAG00000012594","ENSBTAG00000004139","ENSBTAG00000018278","ENSBTAG00000021997","ENSBTAG00000008482")
Tot_lac2 = c("ENSBTAG00000009188","ENSBTAG00000001258","ENSBTAG00000021819","ENSBTAG00000019404","ENSBTAG00000015212")

# convert and orgnize
GeneInfo = convertNformatID(GeneSetNames=c("lactation1","lactation2"),
                            SigGene_list = list(Sig_lac1,Sig_lac2),
                            TotalGene_list = list(Tot_lac1,Tot_lac2),
                            IDtype = "ens") # Need to choose from c('ens','entrez','symbol')

# Resulting an integreted gene identifier object
GeneInfo
#> $lactation1
#>                 Gene ENTREZID SYMBOL SYMBOL_Suggested Sig
#> 1 ENSBTAG00000012594   615431  MRPS6            MRPS6   1
#> 2 ENSBTAG00000004139   531350  BACH1            BACH1   1
#> 3 ENSBTAG00000018278   281640 ATP5PO           ATP5PO   0
#> 4 ENSBTAG00000021997   510879  ITSN1            ITSN1   0
#> 5 ENSBTAG00000008482   516462    SON              SON   0
#> 
#> $lactation2
#>                 Gene ENTREZID  SYMBOL SYMBOL_Suggested Sig
#> 1 ENSBTAG00000009188   281183    GART             GART   1
#> 2 ENSBTAG00000001258   615627 TMEM50B          TMEM50B   1
#> 3 ENSBTAG00000021819   282257  IFNAR1           IFNAR1   0
#> 4 ENSBTAG00000019404   767864  IL10RB           IL10RB   0
#> 5 ENSBTAG00000015212   282258  IFNAR2           IFNAR2   0
```

With simply providing **significant/total gened** as **list** R objects,
***convertNformatID()*** automatically match and organize genes across
different identifiers, namely, coordinating **Ensembl Gene ID**,
**EntrezID**, **Gene Symbol** and **HGNC suggested symbol** if
discrepancy was found. Also, an additional column indicating
significance status will be added (1 stands for significant and 0 for
insignificant).

With great cohesion, objects resulted from last step (e.g. **GeneInfo**)
could be fed into the subsequent loop processes.

With six databases been build-in beforehand, users can simply indicate
which database to use by providing a parameter - **Database = “xxx”** in
the function arguments.

``` r
# Enrichment of each database might take a few mintues to finish.
HyperGEnrich(GeneSet = GeneInfo,
             Database = 'kegg', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 0.1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = F)
```

As results, the function does not return any object in the R
environment, however, all the results would be packed into an *.RData*
object and saved in the current working directory.

There are two elements in the resulting *.RData* object:

  - **results** contains significant results based on the prometers
    provided.
  - **results\_raw** records every single term/pathway tested

<!-- end list -->

``` r
# Here is a demo of results format
data(SampleResults)
class(results) # it's a list
#> [1] "list"
length(results) # number of elements equals to numbers of (significant) gene list provided
#> [1] 3
dim(results[[1]]) # e.g. here are 144 significant pathways/terms and each has 9 attributes/statistics
#> [1] 144   9
```

``` r
names(results[[1]]) # Specific attributes/statistics
#> [1] "Term"               "totalG"             "sigG"              
#> [4] "pvalue"             "ExternalLoss_total" "ExternalLoss_sig"  
#> [7] "findG"              "hitsPerc"           "adj.pvalue"
```

Totally nine columns will be documented in the outputs:

  - **Term**: term ID and annotation/explanations.
  - **totalG**: ***m*** (number of total genes in the pathway)
  - **sigG**: ***k*** (number of significant genes in the pathway)
  - **pvalue**: pvalue of fisher’s exact test
  - **ExternalLoss\_total**: \# of total genes **NOT** annotated in the
    database
  - **ExternalLoss\_sig**: \# of significant genes **NOT** annotated in
    the database
  - **findG**: enumerating **k**, significant genes found in the pathway
  - **hitsPerc**: **k/m**, percentage of sigG
  - **adj.pvalue**: adjusted pvalues for multiple testing correction

## About Databases

Although databases will be updated on regular bases (tentatively
half-yearly), users are free to request an update or update/download
databases by themselves using the built-in database updating functions.

There are totally six dataset being subject to update.

``` r
# Note that these functions are potentially time-comsuming.
# New databases (in .RData format) will be stored in current working directory
GO_DB_Update()
KEGG_DB_Update()
Interpro_DB_Update()
MeSH_DB_Update()
Msig_DB_Update()
Reactome_DB_Update()
```

Please note, when you use new databases, please make sure:

  - Put the newly gathered database in your current working directory,
    they could not be loaded otherwise.  
  - Make sure to set *NewDB* parameter to **T**.

<!-- end list -->

``` r
HyperGEnrich(GeneSet = GeneInfo,
             Database = 'kegg', #c("go","kegg","interpro","mesh","msig","reactome")
             minOverlap = 4, # minimum overlap of pathway genes and total genes
             pvalue_thres = 0.05, # pvalue of fisher's exact test
             adj_pvalue_thres = 0.1, # adjusted pvalues based on multiple testing correction
             padj_method = "fdr", #c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
             NewDB = T) ### Set to T
```
