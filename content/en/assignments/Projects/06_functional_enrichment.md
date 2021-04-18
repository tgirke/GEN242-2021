---
title: "Functional enrichment analysis (FEA)"
linkTitle: "Enrichment Analysis"
description: >
type: docs
weight: 407
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

## Challenge Project: Functional enrichment analysis (FEA)

+ Run workflow from start to finish (steps 1-8) on ChIP-Seq data set from Kaufman et al. (2010)
+ Challenge project tasks
    + Perform functional enrichment analysis on the genes overlapping or downstream of the peak ranges discovered by the ChIP-Seq workflow.   
    + Compare at least 2 functional enrichment methods (_e.g._ [GOCluster_Report](http://bioconductor.org/packages/devel/bioc/html/systemPipeR.html), [fgsea](https://bioconductor.org/packages/3.12/bioc/html/fgsea.html), chipenrich, goseq, GOstats) using KEGG/Reactome or Gene Ontology as functional annotation systems. Among the FEA methods include one based on  the hypergeometric distribution (ORA) and one on the Gene Set Enrichment Analysis (GSEA) algorithm. Assess the results as follows:
        + Quantify the rank-based similarities of the functional categories among the chosen enrichment methods.
        + Determine whether the enrichment results match the biological expectations of the experiment (e.g. are certain biological processes enriched)?
        + Optional: visualize the results with one of the pathway or GO graph viewing tools. 

## References

+ Kaufmann, K, F Wellmer, J M Muiño, T Ferrier, S E Wuest, V Kumar, A Serrano-Mislata, et al. 2010. “Orchestration of Floral Initiation by APETALA1.” Science 328 (5974): 85–89. [PubMed](https://pubmed.ncbi.nlm.nih.gov/20360106/)
+ Sergushichev A (2016) An algorithm for fast preranked gene set enrichment analysis using cumulative statistic calculation. [bioRxiv 060012](https://www.biorxiv.org/content/10.1101/060012v3)
+ Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, et al (2005) Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A 102: 15545–15550. [PubMed](https://pubmed.ncbi.nlm.nih.gov/16199517/)
+ Welch RP, Lee C, Imbriano PM, Patil S, Weymouth TE, Smith RA, Scott LJ, Sartor MA (2014) ChIP-Enrich: gene set enrichment testing for ChIP-seq data. Nucleic Acids Res 42: e105. [PubMed](https://pubmed.ncbi.nlm.nih.gov/24878920/)
+ Young MD, Wakefield MJ, Smyth GK, Oshlack A (2010) Gene ontology analysis for RNA-seq: accounting for selection bias. Genome Biol 11: R14. [PubMed](https://pubmed.ncbi.nlm.nih.gov/20132535/)



