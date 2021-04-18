---
title: "Embedding Methods for scRNA-Seq 1-3"
linkTitle: "Embedding Methods 1-3"
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

## Challenge Projects

### 1. Embedding Methods for scRNA-Seq 

+ Run workflow from start to finish (steps 1-7) on RNA-Seq data set from Howard et al. (2013)
+ Challenge project tasks
    + Compare the partition performance of at least 3 embedding methods for high-dimensional gene expression data using single cell RNA-Seq data. The dimensionality reduction methods can include PCA, MDS, [SC3](http://bioconductor.org/packages/release/bioc/html/SC3.html), [isomap](https://bioconductor.org/packages/release/bioc/html/RDRToolbox.html), [t-SNE](https://cran.r-project.org/web/packages/Rtsne/), [FIt-SNE](https://github.com/KlugerLab/FIt-SNE), [UMAP](https://cran.r-project.org/web/packages/umap/index.html), etc. 
    + To obtain meaningful test results, choose an scRNA-Seq data set (here pre-processed count data) where the correct cell clustering is known (ground truth). For simplicity the data could be obtained from the scRNAseq package (Risso and Cole, 2020) or loaded from GEO (e.g. Shulse et al., 2019). For learning purposes, organize the data in a SingleCellExperiment object.
    + Optional: plot the (partitioning) performance in form of ROC curves and/or record their AUC values.
    + Compare your test results with published performance test results, e.g. Sun et al. (2019) or Duò et al. (2018).

### 2. Embedding Methods for scRNA-Seq 

+ Similar as above but with different combination of embedding methods and/or performance testing approach.

### 3. Embedding Methods for scRNA-Seq 

+ Similar as above but with different combination of embedding methods and/or performance testing approach.

## References

+ Duò A, Robinson MD, Soneson C (2018) A systematic performance evaluation of clustering methods for single-cell RNA-seq data. F1000Res 7: 1141. [PubMed](https://pubmed.ncbi.nlm.nih.gov/30271584/)
+ Howard, B.E. et al., 2013. High-throughput RNA sequencing of pseudomonas-infected Arabidopsis reveals hidden transcriptome complexity and novel splice variants. PloS one, 8(10), p.e74183. [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24098335)
+ Kiselev VY, Kirschner K, Schaub MT, Andrews T, Yiu A, Chandra T, Natarajan KN, Reik W, Barahona M, Green AR, et al (2017) SC3: consensus clustering of single-cell RNA-seq data. Nat Methods 14: 483–486. [PubMed](https://pubmed.ncbi.nlm.nih.gov/28346451/)
+ L.J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research 9 (Nov) : 2579-2605, 2008. 
+ Linderman GC, Rachh M, Hoskins JG, Steinerberger S, Kluger Y (2019) Fast interpolation-based t-SNE for improved visualization of single-cell RNA-seq data. Nat Methods 16: 243–245 [PubMed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6402590/) (Note: this could be used as a more recent pub on t-SNE; the speed improved version is also available for R with a C)
+ McInnes L, Healy J, Melville J (2018) UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. [arXiv](https://arxiv.org/abs/1802.03426) 
+ Risso D, Cole M (2020). scRNAseq: Collection of Public Single-Cell RNA-Seq Datasets. R package version 2.4.0. -> Choose one scRNA-Seq data set from this Bioc data package for testing embedding methods. [URL](https://bioconductor.org/packages/release/data/experiment/html/scRNAseq.html)
+ Senabouth A, Lukowski SW, Hernandez JA, Andersen SB, Mei X, Nguyen QH, Powell JE (2019) ascend: R package for analysis of single-cell RNA-seq data. Gigascience. doi: 10.1093/gigascience/giz087. [PubMed](https://pubmed.ncbi.nlm.nih.gov/31505654/)
+ Shulse CN, Cole BJ, Ciobanu D, Lin J, Yoshinaga Y, Gouran M, Turco GM, Zhu Y, O’Malley RC, Brady SM, et al (2019) High-Throughput Single-Cell Transcriptome Profiling of Plant Cell Types. Cell Rep 27: 2241–2247.e4 [PubMed](https://pubmed.ncbi.nlm.nih.gov/31091459/)
+ Sun S, Zhu J, Ma Y, Zhou X (2019) Accuracy, robustness and scalability of dimensionality reduction methods for single-cell RNA-seq analysis. Genome Biol 20: 269. [PubMed](https://pubmed.ncbi.nlm.nih.gov/31823809/)
+ Sun S, Zhu J, Zhou X (2020) Statistical analysis of spatial expression patterns for spatially resolved transcriptomic studies. Nat Methods. doi: 10.1038/s41592-019-0701-7. [PubMed](https://pubmed.ncbi.nlm.nih.gov/31988518/)








