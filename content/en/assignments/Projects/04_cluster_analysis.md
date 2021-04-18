---
title: "Cluster and network analysis methods"
linkTitle: "Cluster Analysis 1"
description: >
type: docs
weight: 405
---

<br></br>

## RNA-Seq Workflow  

1. Read quality assessment, filtering and trimming 
2. Map reads against reference genome 
3. Perform read counting for required ranges (_e.g._ exonic gene ranges)
4. Normalization of read counts
5. Identification of differentially expressed genes (DEGs)
6. Clustering of gene expression profiles 
7. Gene set enrichment analysis

## Challenge Project: Cluster and network analysis methods

+ Run workflow from start to finish (steps 1-7) on RNA-Seq data set from Howard et al. (2013)
+ Challenge project tasks
    + Compare at least 2-3 cluster analysis methods (e.g. hierarchical, k-means, Fuzzy C-Means, WGCNA, other) and assess the performance differences as follows:
        + Analyze the similarities and differences in the cluster groupings obtained from the two methods. 
        + Do the differences affect the results of the downstream gene set enrichment analysis?
        + Plot the performance of the clustering methods in form of ROC curves and/or record their AUC values. Functional annotations (e.g. GO, KEGG, Pfam) could be used as a benchmark for defining true results.

## References

+ Howard, B.E. et al., 2013. High-throughput RNA sequencing of pseudomonas-infected Arabidopsis reveals hidden transcriptome complexity and novel splice variants. PloS one, 8(10), p.e74183. [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24098335)
+ Langfelder P, Luo R, Oldham MC, Horvath S (2011) Is my network module preserved and reproducible? PLoS Comput Biol 7: e1001057. [PubMed](https://pubmed.ncbi.nlm.nih.gov/21283776/)
+ Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9: 559–559. [PubMed](https://pubmed.ncbi.nlm.nih.gov/19114008/)
+ Rodriguez MZ, Comin CH, Casanova D, Bruno OM, Amancio DR, Costa L da F, Rodrigues FA (2019) Clustering algorithms: A comparative approach. PLoS One 14: e0210236. [PubMed](https://pubmed.ncbi.nlm.nih.gov/30645617/)





