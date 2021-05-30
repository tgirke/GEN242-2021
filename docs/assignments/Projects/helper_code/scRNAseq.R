###############################################################
## Embedding Methods for RNA-Seq and scRNA-Seq Samples/Cells ##
###############################################################

## (1) Generate SummarizedExperiment and SingleCellExperiment objects for
## RNA-Seq data from Howard et al (2013)

## (1a) Creaat first Summarized Experiment object and then coerce to SingleCellExperiment
library(SummarizedExperiment); library(SingleCellExperiment)                                                                                                                        
targetspath <- "targetsPE.txt" # Assumes targetsPE file from course project in your current directory                                                                                                                                                     
countpath <- "results/countDFeByg.xls" # Assumes full count table from course project in results sub-directory                                                                                                                                            
targets <- read.delim(targetspath, comment.char = "#")                                                                                                                              
rownames(targets) <- targets$SampleName                                                                                                                                             
countDF <- read.delim(countpath, row.names=1, check.names=FALSE)                                                                                                                    
(se <- SummarizedExperiment(assays=list(counts=countDF), colData=targets))                                                                                                          
(sce <- as(se, "SingleCellExperiment"))
## (1a) The SingleCellExperiment can also be created directly from the two tabular input files like this
(sce2 <- SingleCellExperiment(assays=list(counts=countDF), colData=targets))

## (2) Prepare data for plotting with run* embedding functions from scran/scater
library(scran); library(scater)
sce <- logNormCounts(sce)
colLabels(sce) <- factor(colData(sce)$Factor) # This uses replicate info from above targets file as pseudo-clusters

## (3) Run and plot with different embedding methods. Note: the embedding
## results are sequentially appended to the SingleCellExperiment object,
## meaning one can use the plot function whenever necessary.

## (3a) Run tSNE on Howard et al data
sce <- runTSNE(sce)
reducedDimNames(sce)
plotTSNE(sce, colour_by="label", text_by="label")

## (3b) Run MDS on Howard et al data
sce <- runMDS(sce)
reducedDimNames(sce)
plotMDS(sce, colour_by="label", text_by="label")

## (3c) Run UMAP on Howard et al data
sce <- runUMAP(sce) # gives a warning but it still works 
reducedDimNames(sce)
plotUMAP(sce, colour_by="label", text_by="label")

## (3d) Run PCA on Howard et al data
sce <- runPCA(sce) # gives a warning but it still works 
reducedDimNames(sce)
plotPCA(sce, colour_by="label", text_by="label")

## Importantly: now compare the results in the four embedding plots with the sample
## tree from the RNA-Seq workflow here: https://bit.ly/3c3z1j5
