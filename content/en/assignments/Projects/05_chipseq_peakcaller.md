---
title: "ChIP-Seq Peak Callers"
linkTitle: "ChIP-Seq Peak Callers"
description: >
type: docs
weight: 406
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

## Challenge Projects

### 1. Comparison of peak calling methods

+ Run workflow from start to finish (steps 1-8) on ChIP-Seq data set from Kaufman et al. (2010)
+ Challenge project tasks
    + Call peaks with at least 2-3 software tools, such as MACS2, PeakSeq, F-Seq, Homer, ChIPseqR, or CSAR.
	+ Compare the results from with peaks identified by Kaufmann et al (2010)
	+ Report unique and common peaks among three methods and plot the results as venn diagrams
	+ Plot the performance of the peak callers in form of ROC plots. As true result set one can use the intersect of the peaks identified by all methods.

### 2. Comparison of peak calling methods

+ Similar as above but with different combination of peak calling methods and/or performance testing approach.

## References

+ Feng J, Liu T, Qin B, Zhang Y, Liu XS (2012) Identifying ChIP-seq enrichment using MACS. Nat Protoc 7: 1728–1740. [PubMed](https://pubmed.ncbi.nlm.nih.gov/22936215/)
+ Kaufmann, K, F Wellmer, J M Muiño, T Ferrier, S E Wuest, V Kumar, A Serrano-Mislata, et al. 2010. “Orchestration of Floral Initiation by APETALA1.” Science 328 (5974): 85–89. [PubMed](https://pubmed.ncbi.nlm.nih.gov/20360106/)
+ Landt SG, Marinov GK, Kundaje A, Kheradpour P, Pauli F, Batzoglou S, Bernstein BE, Bickel P, Brown JB, Cayting P, et al (2012) ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. Genome Res 22: 1813–1831. [PubMed](https://pubmed.ncbi.nlm.nih.gov/22955991/)
+ Lun ATL, Smyth GK (2014) De novo detection of differentially bound regions for ChIP-seq data using peaks and windows: controlling error rates correctly. Nucleic Acids Res 42: e95. [PubMed](https://pubmed.ncbi.nlm.nih.gov/24852250/)
+ Muiño JM, Kaufmann K, van Ham RC, Angenent GC, Krajewski P (2011) ChIP-seq Analysis in R (CSAR): An R package for the statistical detection of protein-bound genomic regions. Plant Methods 7: 11. [PubMed](https://pubmed.ncbi.nlm.nih.gov/21554688/)
+ Wilbanks EG, Facciotti MT (2010) Evaluation of algorithm performance in ChIP-seq peak detection. PLoS One. doi: 10.1371/journal.pone.0011471. [PubMed](https://pubmed.ncbi.nlm.nih.gov/20628599/)
+ Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nussbaum C, Myers RM, Brown M, Li W, et al (2008) Model-based analysis of ChIP-Seq (MACS). Genome Biol. doi: 10.1186/gb-2008-9-9-r137. [PubMed](https://pubmed.ncbi.nlm.nih.gov/18798982/)





