---
title: Functional Enrichment Analysis 
author: "First/last name (first.last@ucr.edu)"
date: "Last update: 05 May, 2021" 
output:
  html_document:
    toc: true
    toc_float:
        collapsed: true
        smooth_scroll: true
    toc_depth: 3
    fig_caption: yes
    code_folding: show
    number_sections: true

fontsize: 14pt
bibliography: bibtex.bib
weight: 11
type: docs
---

<!--
- Compile from command-line
Rscript -e "rmarkdown::render('rfea.Rmd', c('html_document'), clean=F); knitr::knit('rfea.Rmd', tangle=TRUE)"
-->
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  document.querySelector("h1").className = "title";
});
</script>
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  var links = document.links;  
  for (var i = 0, linksLength = links.length; i < linksLength; i++)
    if (links[i].hostname != window.location.hostname)
      links[i].target = '_blank';
});
</script>

<div style="text-align: right">

Source code downloads:    
\[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rfea/rfea.Rmd) \]    
\[ [.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rfea/rfea.R) \]

</div>

## Introduction

The following introduces several widely used gene and protein annotation
systems that are commonly used for functional enrichment analysis (FEA). These
include among many other annotation systems: Gene Ontology (GO), Disease
Ontology (DO) and pathway annotations, such as KEGG and Reactome. Examples of
widely used statistical enrichment methods are introduced as well. These
statistical FEA methods assess whether functional annotation terms are
over-represented in a query gene set. In case of so called over-represention
analysis (ORA) methods the query is a list of unranked gene identifiers
(Falcon and Gentleman 2007). In contrast to this, Gene Set Enrichment Analysis (GSEA)
algorithms use as query score ranked lists, *e.g.* all genes profiled by an assay
(Subramanian et al. 2005; Sergushichev 2016; Duan et al. 2020). The actual sets can
be composed of genes, proteins and even compounds. For simplicity. the term
gene sets is used throughtout this text.

## Functional Annotations Systems

This section introduces a small selection of functional annotation systems provided
by Bioconductor. This includes code to inspect how the annotation data is organized
and how to access it.

## Gene Ontology DB

``` r
## Load GOstats library
library(GOstats); library(GO.db)
## Print complete GO term information for "GO:0003700"
GOTERM$"GO:0003700"
## Print parent and children terms for a GO ID
GOMFPARENTS$"GO:0003700"; GOMFCHILDREN$"GO:0003700"
## Print complete lineages of parents and children for a GO ID
GOMFANCESTOR$"GO:0003700"; GOMFOFFSPRING$"GO:0003700"
## Print number of GO terms in each of the 3 ontologies
zz <- eapply(GOTERM, function(x) x@Ontology); table(unlist(zz))
## Gene to GO mappings for an organism (here Arabidopsis)
library(org.At.tair.db) # For human use org.Hs.eg.db
xx <- as.list(org.At.tairGO2ALLTAIRS)
```

## Pathways

### KEGG

#### `KEGG.db`

### Reactome

#### `reactome.db`

#### `ReactomeContentService4R`

## Functional Enrichment Analysis Methods

### Over-representation analysis (ORA)

#### `GOstats` Package

``` r
## Load required packages
library(GOstats); library(GO.db); library(org.At.tair.db)
## Define universe and test sample set
geneUniverse <- keys(org.At.tairGENENAME)
geneSample <- c("AT2G46210", "AT2G19880", "AT2G38910", "AT5G25140", "AT2G44525")
## Generate params object
params <- new("GOHyperGParams", geneIds = geneSample,
                universeGeneIds = geneUniverse,
                annotation="org.At.tair", ontology = "MF", pvalueCutoff = 0.5,
                conditional = FALSE, testDirection = "over")
## Run enrichment test
hgOver <- hyperGTest(params)
## Viewing of results
summary(hgOver)[1:4,]
htmlReport(hgOver, file = "MyhyperGresult.html")
```

#### `GOHyperGAll` and `GOCluster_Report`

``` r
library(systemPipeR)
## Obtain annotations from BioMart
listMarts() # To choose BioMart database
m <- useMart("ENSEMBL_MART_PLANT"); listDatasets(m) 
m <- useMart("ENSEMBL_MART_PLANT", dataset="athaliana_eg_gene")
listAttributes(m) # Choose data types you want to download
go <- getBM(attributes=c("go_accession", "tair_locus", "go_namespace_1003"), mart=m)
go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
write.table(go, "GOannotationsBiomart_mod.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

## Create catDB instance (takes a while but needs to be done only once)
catdb <- makeCATdb(myfile="GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)
catdb

## Create catDB from Bioconductor annotation package
# catdb <- makeCATdb(myfile=NULL, lib="ath1121501.db", org="", colno=c(1,2,3), idconv=NULL)

## AffyID-to-GeneID mappings when working with AffyIDs 
# affy2locusDF <- systemPipeR:::.AffyID2GeneID(map = "ftp://ftp.arabidopsis.org/home/tair/Microarrays/Affymetrix/affy_ATH1_array_elements-2010-12-20.txt", download=TRUE)
# catdb_conv <- makeCATdb(myfile="GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=list(affy=affy2locusDF))
# systemPipeR:::.AffyID2GeneID(catdb=catdb_conv, affyIDs=c("244901_at", "244902_at"))

## Next time catDB can be loaded from file
save(catdb, file="catdb.RData") 
load("catdb.RData")

## Perform enrichment test on single gene set
test_sample <- unique(as.character(catmap(catdb)$D_MF[1:100,"GeneID"]))
GOHyperGAll(catdb=catdb, gocat="MF", sample=test_sample, Nannot=2)[1:20,]

## GO Slim analysis by subsetting results accordingly
GOHyperGAll_result <- GOHyperGAll(catdb=catdb, gocat="MF", sample=test_sample, Nannot=2)
GOHyperGAll_Subset(catdb, GOHyperGAll_result, sample=test_sample, type="goSlim") 

## Reduce GO term redundancy in 'GOHyperGAll_results'
simplifyDF <- GOHyperGAll_Simplify(GOHyperGAll_result, gocat="MF", cutoff=0.001, correct=T)
# Returns the redundancy reduced data set. 
data.frame(GOHyperGAll_result[GOHyperGAll_result[,1] 

## Batch Analysis of Gene Clusters
testlist <- list(Set1=test_sample)
GOBatchResult <- GOCluster_Report(catdb=catdb, setlist=testlist, method="all", id_type="gene", CLSZ=10, cutoff=0.001, gocats=c("MF", "BP", "CC"), recordSpecGO=c("GO:0003674", "GO:0008150", "GO:0005575"))

## Plot 'GOBatchResult' as bar plot
goBarplot(GOBatchResult, gocat="MF")
```

### Set enrichment analysis (SEA)

#### `fgsea` Package

## Version Information

``` r
sessionInfo()
```

    ## R version 4.0.5 (2021-03-31)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Debian GNU/Linux 10 (buster)
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
    ##  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    ## [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] fgsea_1.16.0     ggplot2_3.3.2    BiocStyle_2.18.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.5          bslib_0.2.4         compiler_4.0.5      pillar_1.4.7       
    ##  [5] BiocManager_1.30.10 jquerylib_0.1.3     tools_4.0.5         digest_0.6.27      
    ##  [9] lattice_0.20-41     jsonlite_1.7.1      evaluate_0.14       lifecycle_0.2.0    
    ## [13] tibble_3.0.4        gtable_0.3.0        pkgconfig_2.0.3     rlang_0.4.8        
    ## [17] Matrix_1.3-2        fastmatch_1.1-0     parallel_4.0.5      yaml_2.2.1         
    ## [21] blogdown_1.2        xfun_0.22           gridExtra_2.3       withr_2.3.0        
    ## [25] stringr_1.4.0       dplyr_1.0.2         knitr_1.30          generics_0.1.0     
    ## [29] sass_0.3.1          vctrs_0.3.5         grid_4.0.5          tidyselect_1.1.0   
    ## [33] data.table_1.13.2   glue_1.4.2          R6_2.5.0            BiocParallel_1.24.1
    ## [37] rmarkdown_2.7       bookdown_0.21       purrr_0.3.4         magrittr_2.0.1     
    ## [41] codetools_0.2-18    scales_1.1.1        htmltools_0.5.1.1   ellipsis_0.3.1     
    ## [45] colorspace_2.0-0    stringi_1.5.3       munsell_0.5.0       crayon_1.3.4

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Duan2020-wz" class="csl-entry">

Duan, Yuzhu, Daniel S Evans, Richard A Miller, Nicholas J Schork, Steven R Cummings, and Thomas Girke. 2020. “<span class="nocase">signatureSearch: environment for gene expression signature searching and functional interpretation</span>.” *Nucleic Acids Res.*, October. <https://doi.org/10.1093/nar/gkaa878>.

</div>

<div id="ref-Falcon2007-eb" class="csl-entry">

Falcon, S, and R Gentleman. 2007. “<span class="nocase">Using GOstats to test gene lists for GO term association</span>.” *Bioinformatics* 23 (2): 257–58. <https://doi.org/10.1093/bioinformatics/btl567>.

</div>

<div id="ref-Sergushichev2016-ms" class="csl-entry">

Sergushichev, Alexey. 2016. “<span class="nocase">An algorithm for fast preranked gene set enrichment analysis using cumulative statistic calculation</span>.” *bioRxiv*. <https://doi.org/10.1101/060012>.

</div>

<div id="ref-Subramanian2005-kx" class="csl-entry">

Subramanian, A, P Tamayo, V K Mootha, S Mukherjee, B L Ebert, M A Gillette, A Paulovich, et al. 2005. “<span class="nocase">Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles</span>.” *Proc. Natl. Acad. Sci. U. S. A.* 102 (43): 15545–50. <https://doi.org/10.1073/pnas.0506580102>.

</div>

</div>
