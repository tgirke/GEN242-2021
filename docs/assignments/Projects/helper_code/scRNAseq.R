###############################################################
## Embedding Methods for RNA-Seq and scRNA-Seq Samples/Cells ##
###############################################################

## (1) Generate SummarizedExperiment and SingleCellExperiment objects for
## RNA-Seq data from Howard et al (2013)

###############################################################
## Embedding Methods for RNA-Seq and scRNA-Seq Samples/Cells ##
###############################################################

## (A) Bulk RNA-Seq Data

## (A1) Generate SummarizedExperiment and SingleCellExperiment objects for
## RNA-Seq data from Howard et al (2013)

## (A1a) Creaat first Summarized Experiment object and then coerce to SingleCellExperiment
library(SummarizedExperiment); library(SingleCellExperiment)                                                                                                                        
targetspath <- "targetsPE.txt"                                                                                                                                                      
countpath <- "results/countDFeByg.xls"                                                                                                                                              
targets <- read.delim(targetspath, comment.char = "#")                                                                                                                              
rownames(targets) <- targets$SampleName                                                                                                                                             
countDF <- read.delim(countpath, row.names=1, check.names=FALSE)                                                                                                                    
(se <- SummarizedExperiment(assays=list(counts=countDF), colData=targets))                                                                                                          
(sce <- as(se, "SingleCellExperiment"))
## (1a) The SingleCellExperiment can also be created directly from the two tabular input files like this
(sce2 <- SingleCellExperiment(assays=list(counts=countDF), colData=targets))

## (A2) Prepare data for plotting with run* embedding functions from scran
library(scran); library(scater)
sce <- logNormCounts(sce)
colLabels(sce) <- factor(colData(sce)$Factor) # This uses replicate info from above targets file as pseudo-clusters

## (A3) Run and plot with different embedding methods Note: the embedding
## results are sequentially appended to the SingleCellExperiment object,
## meaning one can use the plot function whenever necessary.

## (A3a) Run tSNE on Howard et al data
sce <- runTSNE(sce)
reducedDimNames(sce)
plotTSNE(sce, colour_by="label", text_by="label")

## (A3b) Run MDS on Howard et al data
sce <- runMDS(sce)
reducedDimNames(sce)
plotMDS(sce, colour_by="label", text_by="label")

## (A3c) Run UMAP on Howard et al data
sce <- runUMAP(sce) 
reducedDimNames(sce)
plotUMAP(sce, colour_by="label", text_by="label")

## (A3d) Run PCA on Howard et al data
sce <- runPCA(sce) # gives a warning due to small size of data set but it still works 
reducedDimNames(sce)
plotPCA(sce, colour_by="label", text_by="label")


## (B) Single Cell RNA-Seq Data

## (B1) Obtain scRNA-Seq Data from Xenopus tail (Aztekin et al, 2019)
library(scRNAseq)
sce <- AztekinTailData()

## (B2) Prepare data for plotting with run* embedding functions from scran
library(scran); library(scater)
sce <- logNormCounts(sce)
clusters <- quickCluster(sce)
# sce <- computeSumFactors(sce, clusters=clusters)
colLabels(sce) <- factor(clusters)
table(colLabels(sce))

## (B3) Run and plot with different embedding methods Note: the embedding
## results are sequentially appended to the SingleCellExperiment object,
## meaning one can use the plot function whenever necessary.

## (B3a) Run tSNE on Xenopus scRNA-Seq data
sce <- runTSNE(sce)
reducedDimNames(sce)
plotTSNE(sce, colour_by="label", text_by="label")

## (B3b) Run MDS on Xenopus scRNA-Seq data
sce <- runMDS(sce)
reducedDimNames(sce)
plotMDS(sce, colour_by="label", text_by="label")

## (B3c) Run UMAP on Xenopus scRNA-Seq data
# sce <- runUMAP(sce) # Note, the UMAP embedding is already stored in downloaded SingleCellExperiment object by authers. So one can just use this one or recompute it. 
reducedDimNames(sce)
plotUMAP(sce, colour_by="label", text_by="label")

## (B3d) Run PCA on Xenopus scRNA-Seq data
sce <- runPCA(sce) 
reducedDimNames(sce)
plotPCA(sce, colour_by="label", text_by="label")



