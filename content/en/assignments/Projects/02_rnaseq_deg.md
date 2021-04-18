---
title: "RNA-Seq - DEG Analysis Methods 1-2"
linkTitle: "RNA-Seq - DEG 1-2"
description: >
type: docs
weight: 402
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

## Challenge Projects

### 1. Comparison of DEG analysis methods

+ Run workflow from start to finish (steps 1-7) on RNA-Seq data set from Howard et al. (2013)
+ Challenge project tasks
    + Compare the DEG analysis method chosen for paper presentation with at least 1-2 additional methods (_e.g._ one student compares _edgeR_ _vs._ _baySeq_, and other student _DESeq2_ _vs._ _limma/voom_). Assess the results as follows:
        + Analyze the similarities and differences in the DEG lists obtained from the two methods using intersect matrices, venn diagrams and/or upset plots.
        + Assess the impact of the DEG method on the downstream gene set enrichment analysis?
        + Plot the performance of the DEG methods in form of ROC curves and/or record their AUC values. A consensus DEG set or the one from the Howard et al. (2013) paper could be used as a ‘pseudo’ ground truth result.

### 2. Comparison of DEG analysis methods

+ Similar as above but with different combination of DEG methods and/or performance testing approach.

## References

+ Howard, B.E. et al., 2013. High-throughput RNA sequencing of pseudomonas-infected Arabidopsis reveals hidden transcriptome complexity and novel splice variants. PloS one, 8(10), p.e74183. [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24098335)
+ Guo Y, Li C-I, Ye F, Shyr Y (2013) Evaluation of read count based RNAseq analysis methods. BMC Genomics 14 Suppl 8: S2 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24564449)
+ Hardcastle TJ, Kelly KA (2010) baySeq: empirical Bayesian methods for identifying differential expression in sequence count data. BMC Bioinformatics 11: 422 [PubMed](https://pubmed.ncbi.nlm.nih.gov/20698981/)
+ Liu R, Holik AZ, Su S, Jansz N, Chen K, Leong HS, Blewitt ME, Asselin-Labat M-L, Smyth GK, Ritchie ME (2015) Why weight? Modelling sample and observational level variability improves power in RNA-seq analyses. Nucleic Acids Res. doi: 10.1093/nar/gkv412. [PubMed](https://pubmed.ncbi.nlm.nih.gov/25925576/)
+ Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15: 550 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/25516281)
+ Zhou X, Lindsay H, Robinson MD (2014) Robustly detecting differential expression in RNA sequencing data using observation weights. Nucleic Acids Res 42: e91 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24753412)




