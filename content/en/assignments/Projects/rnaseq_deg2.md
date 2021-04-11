---
title: "RNA-Seq - DEG analysis methods 2"
linkTitle: "RNA-Seq-DEG2"
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

## Challenge Project: Comparison of DEG analysis methods

+ Run workflow from start to finish (steps 1-7) on RNA-Seq data set from Howard et al. (2013)
+ Challenge project tasks
    + Each student compares at least 2 RNA-Seq DEG analysis methods (e.g. edgeR vs baySeq and DESeq2 vs limma/voom) and assesses the results as follows:
        + Analyze the similarities and differences in the DEG lists obtained from the two methods
        + Does it affect the results from the downstream gene set enrichment analysis?
        + Plot the performance of the DEG methods in form of an ROC curve. The DEG set from the Howard et al., 2013 paper could be used as benchmark (true result). 

## References

+ Howard, B.E. et al., 2013. High-throughput RNA sequencing of pseudomonas-infected Arabidopsis reveals hidden transcriptome complexity and novel splice variants. PloS one, 8(10), p.e74183. [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24098335)
+ Guo Y, Li C-I, Ye F, Shyr Y (2013) Evaluation of read count based RNAseq analysis methods. BMC Genomics 14 Suppl 8: S2 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24564449)
+ Hardcastle TJ, Kelly KA (2010) baySeq: empirical Bayesian methods for identifying differential expression in sequence count data. BMC Bioinformatics 11: 422 [PubMed](https://pubmed.ncbi.nlm.nih.gov/20698981/)
+ Liu R, Holik AZ, Su S, Jansz N, Chen K, Leong HS, Blewitt ME, Asselin-Labat M-L, Smyth GK, Ritchie ME (2015) Why weight? Modelling sample and observational level variability improves power in RNA-seq analyses. Nucleic Acids Res. doi: 10.1093/nar/gkv412. [PubMed](https://pubmed.ncbi.nlm.nih.gov/25925576/)
+ Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15: 550 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/25516281)
+ Zhou X, Lindsay H, Robinson MD (2014) Robustly detecting differential expression in RNA sequencing data using observation weights. Nucleic Acids Res 42: e91 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24753412)



