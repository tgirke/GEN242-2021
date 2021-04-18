---
title: "Genome Summary Graphics"
linkTitle: "Genome Graphics"
description: >
type: docs
weight: 409
---

<br></br>

## ChIP-Seq Workflow  

1. Read quality assessment, filtering and trimming
2. Align reads to reference genome
3. Compute read coverage across genome
4. Peak calling with different methods and consensus peak identification
5. Annotate peaks
6. Differential binding analysis
7. Gene set enrichment analysis
8. Motif prediction to identify putative TF binding sites

## Challenge Project: Programmable graphics for visualizing genomic features

+ Run workflow from start to finish (steps 1-8) on ChIP-Seq data set from Kaufman et al. (2010)
+ Challenge project tasks
    + This project focuses on the visualization of patterns in NGS experiments (_e.g._ consensus motifs in ChIP-Seq peaks) to discover novel features in genomes. The visualization backend should be based on one of the programmable and extendable R/Bioconductor environments such as [ggplot2](https://ggplot2.tidyverse.org/) ([ggplotly](https://plotly.com/ggplot2/)), [ggbio](https://www.bioconductor.org/packages/release/bioc/html/ggbio.html), [Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html), [RCircos](https://cran.r-project.org/web/packages/RCircos/index.html), etc. For instance, this could include:
        + The generation of motif logos (_e.g._ for ChIP-Seq peaks) for any number of sequence ranges of interest. 
        + Integration of the results with functional annotation information (_e.g._ protein families from Pfam, exonic regions coding for disordered structures), pathways and/or GO. 
        + Incorporation of quantitative information such as relative or differential abundance information obtained from the corresponding NGS profiling technology. 
        + If there is interest, a Shiny App could be included to run the developed R functions interactively from a web browser. 
    
## References

+ Hahne F, Ivanek R (2016). “Statistical Genomics: Methods and Protocols.” In Mathé E, Davis S (eds.), chapter Visualizing Genomic Data Using Gviz and Bioconductor, 335–351. Springer New York, New York, NY. ISBN 978-1-4939-3578-9, doi: 10.1007/978-1-4939-3578-9_16. [PubMed](https://link.springer.com/protocol/10.1007%2F978-1-4939-3578-9_16)
+ Kaufmann, K, F Wellmer, J M Muiño, T Ferrier, S E Wuest, V Kumar, A Serrano-Mislata, et al. 2010. “Orchestration of Floral Initiation by APETALA1.” Science 328 (5974): 85–89. [PubMed](https://pubmed.ncbi.nlm.nih.gov/20360106/)
+ Yin T, Cook D, Lawrence M (2012). “ggbio: an R package for extending the grammar of graphics for genomic data.” Genome Biology, 13(8), R77. [PubMed](https://pubmed.ncbi.nlm.nih.gov/22937822/) 
+ Zhang H, Meltzer P, Davis S (2013) RCircos: an R package for Circos 2D track plots. BMC Bioinformatics 14: 244–244. [PubMed](https://pubmed.ncbi.nlm.nih.gov/23937229/)


