## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----create_se_sce1a, eval=TRUE, message=FALSE, warning=FALSE-----------------
library(SummarizedExperiment); library(SingleCellExperiment)                                                                                                                        
targetspath <- "results/targetsPE.txt"                                                                                                                                                      
countpath <- "results/countDFeByg.xls"                                                                                                                                              
targets <- read.delim(targetspath, comment.char = "#")                                                                                                                              
rownames(targets) <- targets$SampleName                                                                                                                                             
countDF <- read.delim(countpath, row.names=1, check.names=FALSE)                                                                                                                    
(se <- SummarizedExperiment(assays=list(counts=countDF), colData=targets))                                                                                                          
(sce <- as(se, "SingleCellExperiment"))


## ----create_se_sce1b, eval=TRUE-----------------------------------------------
sce2 <- SingleCellExperiment(assays=list(counts=countDF), colData=targets)


## ----preprocess1, eval=TRUE, message=FALSE, warning=FALSE---------------------
library(scran); library(scater)
sce <- logNormCounts(sce)
colLabels(sce) <- factor(colData(sce)$Factor) # This uses replicate info from above targets file as pseudo-clusters


## ----run_tsne1, eval=TRUE-----------------------------------------------------
sce <- runTSNE(sce)
reducedDimNames(sce)
plotTSNE(sce, colour_by="label", text_by="label")


## ----run_mds1, eval=TRUE------------------------------------------------------
sce <- runMDS(sce)
reducedDimNames(sce)
plotMDS(sce, colour_by="label", text_by="label")


## ----run_umap1, eval=TRUE-----------------------------------------------------
sce <- runUMAP(sce) 
reducedDimNames(sce)
plotUMAP(sce, colour_by="label", text_by="label")


## ----run_pca1, eval=TRUE, message=FALSE, warning=FALSE------------------------
sce <- runPCA(sce) # gives a warning due to small size of data set but it still works 
reducedDimNames(sce)
plotPCA(sce, colour_by="label", text_by="label")


## ----create_sce2, eval=FALSE, message=FALSE, warning=FALSE--------------------
## library(scRNAseq)
## sce <- AztekinTailData()


## ----preprocess2, eval=FALSE--------------------------------------------------
## library(scran); library(scater)
## sce <- logNormCounts(sce)
## clusters <- quickCluster(sce)
## # sce <- computeSumFactors(sce, clusters=clusters)
## colLabels(sce) <- factor(clusters)
## table(colLabels(sce))


## ----filter2, eval=FALSE------------------------------------------------------
## filter <- colSums(assays(sce)$counts) >= 10^4
## sce <- sce[, filter]


## ----collor_by_celltype2, eval=TRUE-------------------------------------------
# colLabels(sce) <- colData(sce)$cluster


## ----run_tsne2, eval=FALSE----------------------------------------------------
## sce <- runTSNE(sce)
## reducedDimNames(sce)
## plotTSNE(sce, colour_by="label", text_by="label")


## ----run_mds2, eval=FALSE-----------------------------------------------------
## sce <- runMDS(sce)
## reducedDimNames(sce)
## plotMDS(sce, colour_by="label", text_by="label")


## ----run_umap2, eval=FALSE----------------------------------------------------
## sce <- runUMAP(sce) # Note, the UMAP embedding is already stored in downloaded SingleCellExperiment object by authers. So one can just use this one or recompute it.
## reducedDimNames(sce)
## plotUMAP(sce, colour_by="label", text_by="label")


## ----run_pca2, eval=FALSE-----------------------------------------------------
## sce <- runPCA(sce)
## reducedDimNames(sce)
## plotPCA(sce, colour_by="label", text_by="label")


## ----sessionInfo--------------------------------------------------------------
sessionInfo()

