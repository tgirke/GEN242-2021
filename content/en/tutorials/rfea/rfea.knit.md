---
title: Functional Enrichment Analysis 
author: "First/last name (first.last@ucr.edu)"
date: "Last update: 02 May, 2021" 
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
weight: 9
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
Source code downloads: &nbsp; &nbsp;
[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rfea/rfea.Rmd) ] &nbsp; &nbsp; 
[ [.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rfea/rfea.R) ]
</div>


## Introduction

TBA ...

## Version Information


```r
sessionInfo()
```

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
```

## References
