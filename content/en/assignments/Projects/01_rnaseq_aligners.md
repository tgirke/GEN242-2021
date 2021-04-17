---
title: RNA-Seq - NGS Aligners 
linkTitle: "RNA-Seq Aligners"
description: >
type: docs
weight: 401
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

## Challenge Project: Comparison of RNA-Seq Aligners 

+ Run workflow from start to finish (steps 1-7) on RNA-Seq data set from Howard et al. (2013)
+ Challenge project tasks
    + Compare the RNA-Seq aligner HISAT2 with at least 1-2 other aligners, such as  Rsubread, Star or Kallisto. Evaluate the impact of the aligner on the downstream analysis results including:
        + Read counts
        + Differentially expressed genes (DEGs)
        + Generate plots to compare the results efficiently

## References

+ Bray NL, Pimentel H, Melsted P, Pachter L (2016) Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol. doi: 10.1038/nbt.3519 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/27043002)
+ Howard, B.E. et al., 2013. High-throughput RNA sequencing of pseudomonas-infected Arabidopsis reveals hidden transcriptome complexity and novel splice variants. PloS one, 8(10), p.e74183. [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24098335)
+ Kim D, Pertea G, Trapnell C, Pimentel H, Kelley R, Salzberg SL (2013) TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions. Genome Biol. doi: 10.1186/gb-2013-14-4-r36 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/23618408)
+ Kim D, Langmead B, Salzberg SL (2015) HISAT: a fast spliced aligner with low memory requirements. Nat Methods 12: 357–360 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/25751142)
+ Liao Y, Smyth GK, Shi W (2013) The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Res 41: e108 [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/23558742)







