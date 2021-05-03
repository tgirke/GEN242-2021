---
title: "ChIP-Seq Workflow Template" 
author: "Author: Daniela Cassol, Le Zhang and Thomas Girke"
date: "Last update: 02 May, 2021" 
output:
  BiocStyle::html_document:
    toc_float: true
    code_folding: show
package: systemPipeR
vignette: |
  %\VignetteEncoding{UTF-8}
  %\VignetteIndexEntry{WF: ChIP-Seq Workflow Template}
  %\VignetteEngine{knitr::rmarkdown}
fontsize: 14pt
bibliography: bibtex.bib
weight: 9
type: docs
---

<style type="text/css">
pre code {
white-space: pre !important;
overflow-x: scroll !important;
word-break: keep-all !important;
word-wrap: initial !important;
}
</style>
<!--
- Compile from command-line
Rscript -e "rmarkdown::render('spchipseq.Rmd', c('BiocStyle::html_document'), clean=F); knitr::knit('spchipseq.Rmd', tangle=TRUE)"; Rscript -e "rmarkdown::render('spchipseq.Rmd', c('BiocStyle::pdf_document'))"
-->
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  document.querySelector("h1").className = "title";
});
</script>
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  var links = document.links;  
  for (var i = 0, linksLength = links.length; i < linksLength; i++)
    if (links[i].hostname != window.location.hostname)
      links[i].target = '_blank';
});
</script>

<div style="text-align: right">

Source code downloads:    
\[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/spchipseq/spchipseq.Rmd) \]    
\[ [.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/spchipseq/spchipseq.R) \]

</div>

## Introduction

Users want to provide here background information about the design of their ChIP-Seq project.

### Background and objectives

This report describes the analysis of several ChIP-Seq experiments
studying the DNA binding patterns of the transcriptions factors … from *organism* ….

### Experimental design

Typically, users want to specify here all information relevant for the
analysis of their NGS study. This includes detailed descriptions of
FASTQ files, experimental design, reference genome, gene annotations,
etc.

### Workflow environment

<font color="red">NOTE: this section</font> describes how to set up the proper environment (directory structure) for running
`systemPipeR` workflows. After mastering this task the workflow run instructions <font color="red">can be deleted</font> since they are not expected
to be included in a final HTML/PDF report of a workflow.

1.  If a remote system or cluster is used, then users need to log in to the
    remote system first. The following applies to an HPC cluster (*e.g.* HPCC
    cluster).

    A terminal application needs to be used to log in to a user’s cluster account. Next, one
    can open an interactive session on a computer node with `srun`. More details about
    argument settings for `srun` are available in this [HPCC
    manual](http://hpcc.ucr.edu/manuals_linux-cluster_jobs.html#partitions) or
    the HPCC section of this website
    [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/linux/linux/#job-submission-with-sbatch).
    Next, load the R version required for running the workflow with `module load`. Sometimes it may be necessary to
    first unload an active software version before loading another version, *e.g.* `module unload R`.

``` sh
srun --x11 --partition=short --mem=2gb --cpus-per-task 4 --ntasks 1 --time 2:00:00 --pty bash -l
module load R/4.0.3
```

2.  Load a workflow template with the `genWorkenvir` function. This can be done from the command-line or from within R.
    However, only one of the two options needs to be used.

From command-line

``` sh
$ Rscript -e "systemPipeRdata::genWorkenvir(workflow='chipseq')"
$ cd chipseq
```

From R

``` r
library(systemPipeRdata)
genWorkenvir(workflow = "chipseq")
setwd("chipseq")
```

3.  Optional: if the user wishes to use another `Rmd` file than the template instance provided by the `genWorkenvir` function, then it can be copied or downloaded
    into the root directory of the workflow environment (*e.g.* with `cp` or `wget`).

4.  Now one can open from the root directory of the workflow the corresponding R Markdown script (*e.g.* systemPipeChIPseq.Rmd) using an R IDE, such as *nvim-r*, *ESS* or RStudio.
    Subsequently, the workflow can be run as outlined below. For learning purposes it is recommended to run workflows for the first time interactively. Once all workflow steps are
    understood and possibly modified to custom needs, one can run the workflow from start to finish with a single command using `rmarkdown::render()` or `runWF()`.

### Load packages

The `systemPipeR` package needs to be loaded to perform the analysis
steps shown in this report (H Backman and Girke 2016). The package allows users
to run the entire analysis workflow interactively or with a single command
while also generating the corresponding analysis report. For details
see `systemPipeR's` main [vignette](http://www.bioconductor.org/packages/devel/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html).

``` r
library(systemPipeR)
```

## Read preprocessing

### Experiment definition provided by `targets` file

The `targets` file defines all FASTQ files and sample comparisons of the analysis workflow.

``` r
targetspath <- system.file("extdata", "targetsPE_chip.txt", package = "systemPipeR")
targets <- read.delim(targetspath, comment.char = "#")
targets[1:4, -c(5, 6)]
```

    ##                     FileName1                   FileName2
    ## 1 ./data/SRR446027_1.fastq.gz ./data/SRR446027_2.fastq.gz
    ## 2 ./data/SRR446028_1.fastq.gz ./data/SRR446028_2.fastq.gz
    ## 3 ./data/SRR446029_1.fastq.gz ./data/SRR446029_2.fastq.gz
    ## 4 ./data/SRR446030_1.fastq.gz ./data/SRR446030_2.fastq.gz
    ##   SampleName Factor        Date SampleReference
    ## 1        M1A     M1 23-Mar-2012                
    ## 2        M1B     M1 23-Mar-2012                
    ## 3        A1A     A1 23-Mar-2012             M1A
    ## 4        A1B     A1 23-Mar-2012             M1B

### Read quality filtering and trimming

The following example shows how one can design a custom read
preprocessing function using utilities provided by the `ShortRead` package, and then
apply it with `preprocessReads` in batch mode to all FASTQ samples referenced in the
corresponding `SYSargs2` instance (`trim` object below). More detailed information on
read preprocessing is provided in `systemPipeR's` main vignette.

First, we construct *`SYSargs2`* object from *`cwl`* and *`yml`* param and *`targets`* files.

``` r
dir_path <- system.file("extdata/cwl/preprocessReads/trim-pe", 
    package = "systemPipeR")
trim <- loadWF(targets = targetspath, wf_file = "trim-pe.cwl", 
    input_file = "trim-pe.yml", dir_path = dir_path)
trim <- renderWF(trim, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
    FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
trim
output(trim)[1:2]
```

Next, we execute the code for trimming all the raw data.

``` r
filterFct <- function(fq, cutoff = 20, Nexceptions = 0) {
    qcount <- rowSums(as(quality(fq), "matrix") <= cutoff, na.rm = TRUE)
    fq[qcount <= Nexceptions]
    # Retains reads where Phred scores are >= cutoff with N
    # exceptions
}
preprocessReads(args = trim, Fct = "filterFct(fq, cutoff=20, Nexceptions=0)", 
    batchsize = 1e+05)
writeTargetsout(x = trim, file = "targets_chip_trimPE.txt", step = 1, 
    new_col = c("FileName1", "FileName2"), new_col_output_index = c(1, 
        2), overwrite = TRUE)
```

### FASTQ quality report

The following `seeFastq` and `seeFastqPlot` functions generate and plot a series
of useful quality statistics for a set of FASTQ files including per cycle quality box
plots, base proportions, base-level quality trends, relative k-mer
diversity, length and occurrence distribution of reads, number of reads
above quality cutoffs and mean quality distribution. The results are
written to a PDF file named `fastqReport.pdf`.
Parallelization of FASTQ quality report via scheduler (*e.g.* Slurm) across several compute nodes.

``` r
library(BiocParallel)
library(batchtools)
f <- function(x) {
    library(systemPipeR)
    targets <- system.file("extdata", "targetsPE_chip.txt", package = "systemPipeR")
    dir_path <- system.file("extdata/cwl/preprocessReads/trim-pe", 
        package = "systemPipeR")
    trim <- loadWorkflow(targets = targets, wf_file = "trim-pe.cwl", 
        input_file = "trim-pe.yml", dir_path = dir_path)
    trim <- renderWF(trim, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
        FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
    seeFastq(fastq = infile1(trim)[x], batchsize = 1e+05, klength = 8)
}

resources <- list(walltime = 120, ntasks = 1, ncpus = 4, memory = 1024)
param <- BatchtoolsParam(workers = 4, cluster = "slurm", template = "batchtools.slurm.tmpl", 
    resources = resources)
fqlist <- bplapply(seq(along = trim), f, BPPARAM = param)

pdf("./results/fastqReport.pdf", height = 18, width = 4 * length(fqlist))
seeFastqPlot(unlist(fqlist, recursive = FALSE))
dev.off()
```

![](../results/fastqReport.png)

<div align="center">

Figure 1: FASTQ quality report for 18 samples

</div>

</br>

## Alignments

### Read mapping with `Bowtie2`

The NGS reads of this project will be aligned with `Bowtie2` against the
reference genome sequence (Langmead and Salzberg 2012). The parameter settings of the
aligner are defined in the `bowtie2-index.cwl` and `bowtie2-index.yml` files.
In ChIP-Seq experiments it is usually more appropriate to eliminate reads mapping
to multiple locations. To achieve this, users want to remove the argument setting
`-k 50 non-deterministic` in the configuration files.

Building the index:

``` r
dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-idx", package = "systemPipeR")
idx <- loadWorkflow(targets = NULL, wf_file = "bowtie2-index.cwl", 
    input_file = "bowtie2-index.yml", dir_path = dir_path)
idx <- renderWF(idx)
idx
cmdlist(idx)

### Run in single machine
runCommandline(idx, make_bam = FALSE)
```

The following submits 18 alignment jobs via a scheduler to a computer cluster.

``` r
targets <- system.file("extdata", "targetsPE_chip.txt", package = "systemPipeR")
dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-pe", package = "systemPipeR")
args <- loadWF(targets = targets, wf_file = "bowtie2-mapping-pe.cwl", 
    input_file = "bowtie2-mapping-pe.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
    FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
args
cmdlist(args)[1:2]
output(args)[1:2]
```

``` r
moduleload(modules(args))  # Skip if a module system is not used
resources <- list(walltime = 120, ntasks = 1, ncpus = 4, memory = 1024)
reg <- clusterRun(args, FUN = runCommandline, more.args = list(args = args, 
    dir = FALSE), conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl", 
    Njobs = 18, runid = "01", resourceList = resources)
getStatus(reg = reg)
waitForJobs(reg = reg)
args <- output_update(args, dir = FALSE, replace = TRUE, extension = c(".sam", 
    ".bam"))  ## Updates the output(args) to the right location in the subfolders
output(args)
```

Alternatively, one can run the alignments sequentially on a single system.

``` r
args <- runCommandline(args, force = F)
```

Check whether all BAM files have been created and write out the new targets file.

``` r
writeTargetsout(x = args, file = "targets_bam.txt", step = 1, 
    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE, 
    remove = TRUE)
outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
file.exists(outpaths)
```

### Read and alignment stats

The following provides an overview of the number of reads in each sample
and how many of them aligned to the reference.

``` r
read_statsDF <- alignStats(args = args)
write.table(read_statsDF, "results/alignStats.xls", row.names = FALSE, 
    quote = FALSE, sep = "\t")
read.delim("results/alignStats.xls")
```

### Create symbolic links for viewing BAM files in IGV

The `symLink2bam` function creates symbolic links to view the BAM alignment files in a
genome browser such as IGV without moving these large files to a local
system. The corresponding URLs are written to a file with a path
specified under `urlfile`, here `IGVurl.txt`. Please replace the directory and the user name.

``` r
symLink2bam(sysargs = args, htmldir = c("~/.html/", "somedir/"), 
    urlbase = "http://cluster.hpcc.ucr.edu/~tgirke/", urlfile = "./results/IGVurl.txt")
```

## Utilities for coverage data

The following introduces several utilities useful for ChIP-Seq data. They are not part of the actual workflow.

### Rle object stores coverage information

``` r
library(rtracklayer)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
aligns <- readGAlignments(outpaths[1])
cov <- coverage(aligns)
cov
```

### Resizing aligned reads

``` r
trim(resize(as(aligns, "GRanges"), width = 200))
```

### Naive peak calling

``` r
islands <- slice(cov, lower = 15)
islands[[1]]
```

### Plot coverage for defined region

``` r
library(ggbio)
myloc <- c("Chr1", 1, 1e+05)
ga <- readGAlignments(outpaths[1], use.names = TRUE, param = ScanBamParam(which = GRanges(myloc[1], 
    IRanges(as.numeric(myloc[2]), as.numeric(myloc[3])))))
autoplot(ga, aes(color = strand, fill = strand), facets = strand ~ 
    seqnames, stat = "coverage")
```

## Peak calling with MACS2

### Merge BAM files of replicates prior to peak calling

Merging BAM files of technical and/or biological replicates can improve
the sensitivity of the peak calling by increasing the depth of read
coverage. The `mergeBamByFactor` function merges BAM files based on grouping information
specified by a `factor`, here the `Factor` column of the imported targets file. It
also returns an updated `SYSargs2` object containing the paths to the
merged BAM files as well as to any unmerged files without replicates.
This step can be skipped if merging of BAM files is not desired.

``` r
dir_path <- system.file("extdata/cwl/mergeBamByFactor", package = "systemPipeR")
args <- loadWF(targets = "targets_bam.txt", wf_file = "merge-bam.cwl", 
    input_file = "merge-bam.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName = "_BAM_PATH_", 
    SampleName = "_SampleName_"))

args_merge <- mergeBamByFactor(args = args, overwrite = TRUE)
writeTargetsout(x = args_merge, file = "targets_mergeBamByFactor.txt", 
    step = 1, new_col = "FileName", new_col_output_index = 1, 
    overwrite = TRUE, remove = TRUE)
```

### Peak calling without input/reference sample

MACS2 can perform peak calling on ChIP-Seq data with and without input
samples (Zhang et al. 2008). The following performs peak calling without
input on all samples specified in the corresponding `args` object. Note, due to
the small size of the sample data, MACS2 needs to be run here with the
`nomodel` setting. For real data sets, users want to remove this parameter
in the corresponding `*.param` file(s).

``` r
dir_path <- system.file("extdata/cwl/MACS2/MACS2-noinput/", package = "systemPipeR")
args <- loadWF(targets = "targets_mergeBamByFactor.txt", wf_file = "macs2.cwl", 
    input_file = "macs2.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))

runCommandline(args, make_bam = FALSE, force = T)
outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
file.exists(outpaths)
writeTargetsout(x = args, file = "targets_macs.txt", step = 1, 
    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
```

### Peak calling with input/reference sample

To perform peak calling with input samples, they can be most
conveniently specified in the `SampleReference` column of the initial
`targets` file. The `writeTargetsRef` function uses this information to create a `targets`
file intermediate for running MACS2 with the corresponding input samples.

``` r
writeTargetsRef(infile = "targets_mergeBamByFactor.txt", outfile = "targets_bam_ref.txt", 
    silent = FALSE, overwrite = TRUE)
dir_path <- system.file("extdata/cwl/MACS2/MACS2-input/", package = "systemPipeR")
args_input <- loadWF(targets = "targets_bam_ref.txt", wf_file = "macs2-input.cwl", 
    input_file = "macs2.yml", dir_path = dir_path)
args_input <- renderWF(args_input, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
    FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
cmdlist(args_input)[1]
### Run
args_input <- runCommandline(args_input, make_bam = FALSE, force = T)
outpaths_input <- subsetWF(args_input, slot = "output", subset = 1, 
    index = 1)
file.exists(outpaths_input)
writeTargetsout(x = args_input, file = "targets_macs_input.txt", 
    step = 1, new_col = "FileName", new_col_output_index = 1, 
    overwrite = TRUE)
```

The peak calling results from MACS2 are written for each sample to
separate files in the `results` directory. They are named after the corresponding
files with extensions used by MACS2.

### Identify consensus peaks

The following example shows how one can identify consensus preaks among two peak sets sharing either a minimum absolute overlap and/or minimum relative overlap using the `subsetByOverlaps` or `olRanges` functions, respectively. Note, the latter is a custom function imported below by sourcing it.

``` r
# source('http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/rangeoverlapper.R')
outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)  ## escolher um dos outputs index
peak_M1A <- outpaths["M1A"]
peak_M1A <- as(read.delim(peak_M1A, comment = "#")[, 1:3], "GRanges")
peak_A1A <- outpaths["A1A"]
peak_A1A <- as(read.delim(peak_A1A, comment = "#")[, 1:3], "GRanges")
(myol1 <- subsetByOverlaps(peak_M1A, peak_A1A, minoverlap = 1))
# Returns any overlap
myol2 <- olRanges(query = peak_M1A, subject = peak_A1A, output = "gr")
# Returns any overlap with OL length information
myol2[values(myol2)["OLpercQ"][, 1] >= 50]
# Returns only query peaks with a minimum overlap of 50%
```

## Annotate peaks with genomic context

### Annotation with `ChIPpeakAnno` package

The following annotates the identified peaks with genomic context information using the `ChIPpeakAnno` and `ChIPseeker` packages, respectively (Zhu et al. 2010; Yu, Wang, and He 2015).

``` r
library(ChIPpeakAnno)
library(GenomicFeatures)
dir_path <- system.file("extdata/cwl/annotate_peaks", package = "systemPipeR")
args <- loadWF(targets = "targets_macs.txt", wf_file = "annotate-peaks.cwl", 
    input_file = "annotate-peaks.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))

txdb <- makeTxDbFromGFF(file = "data/tair10.gff", format = "gff", 
    dataSource = "TAIR", organism = "Arabidopsis thaliana")
ge <- genes(txdb, columns = c("tx_name", "gene_id", "tx_type"))
for (i in seq(along = args)) {
    peaksGR <- as(read.delim(infile1(args)[i], comment = "#"), 
        "GRanges")
    annotatedPeak <- annotatePeakInBatch(peaksGR, AnnotationData = genes(txdb))
    df <- data.frame(as.data.frame(annotatedPeak), as.data.frame(values(ge[values(annotatedPeak)$feature, 
        ])))
    outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
    write.table(df, outpaths[i], quote = FALSE, row.names = FALSE, 
        sep = "\t")
}
writeTargetsout(x = args, file = "targets_peakanno.txt", step = 1, 
    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
```

The peak annotation results are written for each peak set to separate
files in the `results` directory. They are named after the corresponding peak
files with extensions specified in the `annotate_peaks.param` file,
here `*.peaks.annotated.xls`.

### Annotation with `ChIPseeker` package

Same as in previous step but using the `ChIPseeker` package for annotating the peaks.

``` r
library(ChIPseeker)
for (i in seq(along = args)) {
    peakAnno <- annotatePeak(infile1(args)[i], TxDb = txdb, verbose = FALSE)
    df <- as.data.frame(peakAnno)
    outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
    write.table(df, outpaths[i], quote = FALSE, row.names = FALSE, 
        sep = "\t")
}
writeTargetsout(x = args, file = "targets_peakanno.txt", step = 1, 
    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
```

Summary plots provided by the `ChIPseeker` package. Here applied only to one sample
for demonstration purposes.

``` r
peak <- readPeakFile(infile1(args)[1])
covplot(peak, weightCol = "X.log10.pvalue.")
outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
peakHeatmap(outpaths[1], TxDb = txdb, upstream = 1000, downstream = 1000, 
    color = "red")
plotAvgProf2(outpaths[1], TxDb = txdb, upstream = 1000, downstream = 1000, 
    xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency")
```

## Count reads overlapping peaks

The `countRangeset` function is a convenience wrapper to perform read counting
iteratively over serveral range sets, here peak range sets. Internally,
the read counting is performed with the `summarizeOverlaps` function from the
`GenomicAlignments` package. The resulting count tables are directly saved to
files, one for each peak set.

``` r
library(GenomicRanges)
dir_path <- system.file("extdata/cwl/count_rangesets", package = "systemPipeR")
args <- loadWF(targets = "targets_macs.txt", wf_file = "count_rangesets.cwl", 
    input_file = "count_rangesets.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))

### Bam Files
targets <- system.file("extdata", "targetsPE_chip.txt", package = "systemPipeR")
dir_path <- system.file("extdata/cwl/bowtie2/bowtie2-pe", package = "systemPipeR")
args_bam <- loadWF(targets = targets, wf_file = "bowtie2-mapping-pe.cwl", 
    input_file = "bowtie2-mapping-pe.yml", dir_path = dir_path)
args_bam <- renderWF(args_bam, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))
args_bam <- output_update(args_bam, dir = FALSE, replace = TRUE, 
    extension = c(".sam", ".bam"))
outpaths <- subsetWF(args_bam, slot = "output", subset = 1, index = 1)

bfl <- BamFileList(outpaths, yieldSize = 50000, index = character())
countDFnames <- countRangeset(bfl, args, mode = "Union", ignore.strand = TRUE)
writeTargetsout(x = args, file = "targets_countDF.txt", step = 1, 
    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
```

## Differential binding analysis

The `runDiff` function performs differential binding analysis in batch mode for
several count tables using `edgeR` or `DESeq2` (Robinson, McCarthy, and Smyth 2010; Love, Huber, and Anders 2014).
Internally, it calls the functions `run_edgeR` and `run_DESeq2`. It also returns
the filtering results and plots from the downstream `filterDEGs` function using
the fold change and FDR cutoffs provided under the `dbrfilter` argument.

``` r
dir_path <- system.file("extdata/cwl/rundiff", package = "systemPipeR")
args_diff <- loadWF(targets = "targets_countDF.txt", wf_file = "rundiff.cwl", 
    input_file = "rundiff.yml", dir_path = dir_path)
args_diff <- renderWF(args_diff, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))

cmp <- readComp(file = args_bam, format = "matrix")
dbrlist <- runDiff(args = args_diff, diffFct = run_edgeR, targets = targets.as.df(targets(args_bam)), 
    cmp = cmp[[1]], independent = TRUE, dbrfilter = c(Fold = 2, 
        FDR = 1))
writeTargetsout(x = args_diff, file = "targets_rundiff.txt", 
    step = 1, new_col = "FileName", new_col_output_index = 1, 
    overwrite = TRUE)
```

## GO term enrichment analysis

The following performs GO term enrichment analysis for each annotated peak set.

``` r
dir_path <- system.file("extdata/cwl/annotate_peaks", package = "systemPipeR")
args <- loadWF(targets = "targets_bam_ref.txt", wf_file = "annotate-peaks.cwl", 
    input_file = "annotate-peaks.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
    FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))

args_anno <- loadWF(targets = "targets_macs.txt", wf_file = "annotate-peaks.cwl", 
    input_file = "annotate-peaks.yml", dir_path = dir_path)
args_anno <- renderWF(args_anno, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))
annofiles <- subsetWF(args_anno, slot = "output", subset = 1, 
    index = 1)
gene_ids <- sapply(names(annofiles), function(x) unique(as.character(read.delim(annofiles[x])[, 
    "geneId"])), simplify = FALSE)
load("data/GO/catdb.RData")
BatchResult <- GOCluster_Report(catdb = catdb, setlist = gene_ids, 
    method = "all", id_type = "gene", CLSZ = 2, cutoff = 0.9, 
    gocats = c("MF", "BP", "CC"), recordSpecGO = NULL)
```

## Motif analysis

### Parse DNA sequences of peak regions from genome

Enrichment analysis of known DNA binding motifs or *de novo* discovery
of novel motifs requires the DNA sequences of the identified peak
regions. To parse the corresponding sequences from the reference genome,
the `getSeq` function from the `Biostrings` package can be used. The
following example parses the sequences for each peak set and saves the
results to separate FASTA files, one for each peak set. In addition, the
sequences in the FASTA files are ranked (sorted) by increasing p-values
as expected by some motif discovery tools, such as `BCRANK`.

``` r
library(Biostrings)
library(seqLogo)
library(BCRANK)
dir_path <- system.file("extdata/cwl/annotate_peaks", package = "systemPipeR")
args <- loadWF(targets = "targets_macs.txt", wf_file = "annotate-peaks.cwl", 
    input_file = "annotate-peaks.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))

rangefiles <- infile1(args)
for (i in seq(along = rangefiles)) {
    df <- read.delim(rangefiles[i], comment = "#")
    peaks <- as(df, "GRanges")
    names(peaks) <- paste0(as.character(seqnames(peaks)), "_", 
        start(peaks), "-", end(peaks))
    peaks <- peaks[order(values(peaks)$X.log10.pvalue., decreasing = TRUE)]
    pseq <- getSeq(FaFile("./data/tair10.fasta"), peaks)
    names(pseq) <- names(peaks)
    writeXStringSet(pseq, paste0(rangefiles[i], ".fasta"))
}
```

### Motif discovery with `BCRANK`

The Bioconductor package `BCRANK` is one of the many tools available for
*de novo* discovery of DNA binding motifs in peak regions of ChIP-Seq
experiments. The given example applies this method on the first peak
sample set and plots the sequence logo of the highest ranking motif.

``` r
set.seed(0)
BCRANKout <- bcrank(paste0(rangefiles[1], ".fasta"), restarts = 25, 
    use.P1 = TRUE, use.P2 = TRUE)
toptable(BCRANKout)
topMotif <- toptable(BCRANKout, 1)
weightMatrix <- pwm(topMotif, normalize = FALSE)
weightMatrixNormalized <- pwm(topMotif, normalize = TRUE)
pdf("results/seqlogo.pdf")
seqLogo(weightMatrixNormalized)
dev.off()
```

![](../results/seqlogo.png)

<div align="center">

Figure 2: One of the motifs identified by `BCRANK`

</div>

</br>

## Version Information

``` r
sessionInfo()
```

    ## R version 4.0.5 (2021-03-31)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Debian GNU/Linux 10 (buster)
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices
    ## [6] utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] systemPipeR_1.24.5          ShortRead_1.48.0           
    ##  [3] GenomicAlignments_1.26.0    SummarizedExperiment_1.20.0
    ##  [5] Biobase_2.50.0              MatrixGenerics_1.2.0       
    ##  [7] matrixStats_0.57.0          BiocParallel_1.24.1        
    ##  [9] Rsamtools_2.6.0             Biostrings_2.58.0          
    ## [11] XVector_0.30.0              GenomicRanges_1.42.0       
    ## [13] GenomeInfoDb_1.26.1         IRanges_2.24.0             
    ## [15] S4Vectors_0.28.0            BiocGenerics_0.36.0        
    ## [17] BiocStyle_2.18.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] colorspace_2.0-0         rjson_0.2.20            
    ##   [3] hwriter_1.3.2            ellipsis_0.3.1          
    ##   [5] bit64_4.0.5              AnnotationDbi_1.52.0    
    ##   [7] xml2_1.3.2               codetools_0.2-18        
    ##   [9] splines_4.0.5            knitr_1.30              
    ##  [11] jsonlite_1.7.1           annotate_1.68.0         
    ##  [13] GO.db_3.12.1             dbplyr_2.0.0            
    ##  [15] png_0.1-7                pheatmap_1.0.12         
    ##  [17] graph_1.68.0             BiocManager_1.30.10     
    ##  [19] compiler_4.0.5           httr_1.4.2              
    ##  [21] backports_1.2.0          GOstats_2.56.0          
    ##  [23] assertthat_0.2.1         Matrix_1.3-2            
    ##  [25] limma_3.46.0             formatR_1.7             
    ##  [27] htmltools_0.5.1.1        prettyunits_1.1.1       
    ##  [29] tools_4.0.5              gtable_0.3.0            
    ##  [31] glue_1.4.2               GenomeInfoDbData_1.2.4  
    ##  [33] Category_2.56.0          dplyr_1.0.2             
    ##  [35] rsvg_2.1                 batchtools_0.9.14       
    ##  [37] rappdirs_0.3.1           V8_3.4.0                
    ##  [39] Rcpp_1.0.5               jquerylib_0.1.3         
    ##  [41] vctrs_0.3.5              blogdown_1.2            
    ##  [43] rtracklayer_1.50.0       xfun_0.22               
    ##  [45] stringr_1.4.0            lifecycle_0.2.0         
    ##  [47] XML_3.99-0.5             edgeR_3.32.0            
    ##  [49] zlibbioc_1.36.0          scales_1.1.1            
    ##  [51] BSgenome_1.58.0          VariantAnnotation_1.36.0
    ##  [53] hms_0.5.3                RBGL_1.66.0             
    ##  [55] RColorBrewer_1.1-2       yaml_2.2.1              
    ##  [57] curl_4.3                 memoise_1.1.0           
    ##  [59] ggplot2_3.3.2            sass_0.3.1              
    ##  [61] biomaRt_2.46.0           latticeExtra_0.6-29     
    ##  [63] stringi_1.5.3            RSQLite_2.2.1           
    ##  [65] genefilter_1.72.0        checkmate_2.0.0         
    ##  [67] GenomicFeatures_1.42.1   DOT_0.1                 
    ##  [69] rlang_0.4.8              pkgconfig_2.0.3         
    ##  [71] bitops_1.0-6             evaluate_0.14           
    ##  [73] lattice_0.20-41          purrr_0.3.4             
    ##  [75] bit_4.0.4                tidyselect_1.1.0        
    ##  [77] GSEABase_1.52.0          AnnotationForge_1.32.0  
    ##  [79] magrittr_2.0.1           bookdown_0.21           
    ##  [81] R6_2.5.0                 generics_0.1.0          
    ##  [83] base64url_1.4            DelayedArray_0.16.0     
    ##  [85] DBI_1.1.0                withr_2.3.0             
    ##  [87] pillar_1.4.7             survival_3.2-10         
    ##  [89] RCurl_1.98-1.2           tibble_3.0.4            
    ##  [91] crayon_1.3.4             BiocFileCache_1.14.0    
    ##  [93] rmarkdown_2.7            jpeg_0.1-8.1            
    ##  [95] progress_1.2.2           locfit_1.5-9.4          
    ##  [97] grid_4.0.5               data.table_1.13.2       
    ##  [99] blob_1.2.1               Rgraphviz_2.34.0        
    ## [101] digest_0.6.27            xtable_1.8-4            
    ## [103] brew_1.0-6               openssl_1.4.3           
    ## [105] munsell_0.5.0            bslib_0.2.4             
    ## [107] askpass_1.1

## Funding

This project is funded by NSF award [ABI-1661152](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1661152).

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-H_Backman2016-bt" class="csl-entry">

H Backman, Tyler W, and Thomas Girke. 2016. “<span class="nocase">systemPipeR: NGS workflow and report generation environment</span>.” *BMC Bioinformatics* 17 (1): 388. <https://doi.org/10.1186/s12859-016-1241-0>.

</div>

<div id="ref-Langmead2012-bs" class="csl-entry">

Langmead, Ben, and Steven L Salzberg. 2012. “Fast Gapped-Read Alignment with Bowtie 2.” *Nat. Methods* 9 (4): 357–59. <https://doi.org/10.1038/nmeth.1923>.

</div>

<div id="ref-Love2014-sh" class="csl-entry">

Love, Michael, Wolfgang Huber, and Simon Anders. 2014. “Moderated Estimation of Fold Change and Dispersion for <span class="nocase">RNA-seq</span> Data with DESeq2.” *Genome Biol.* 15 (12): 550. <https://doi.org/10.1186/s13059-014-0550-8>.

</div>

<div id="ref-Robinson2010-uk" class="csl-entry">

Robinson, M D, D J McCarthy, and G K Smyth. 2010. “edgeR: A Bioconductor Package for Differential Expression Analysis of Digital Gene Expression Data.” *Bioinformatics* 26 (1): 139–40. <https://doi.org/10.1093/bioinformatics/btp616>.

</div>

<div id="ref-Yu2015-xu" class="csl-entry">

Yu, Guangchuang, Li-Gen Wang, and Qing-Yu He. 2015. “ChIPseeker: An R/Bioconductor Package for ChIP Peak Annotation, Comparison and Visualization.” *Bioinformatics* 31 (14): 2382–83. <https://doi.org/10.1093/bioinformatics/btv145>.

</div>

<div id="ref-Zhang2008-pc" class="csl-entry">

Zhang, Y, T Liu, C A Meyer, J Eeckhoute, D S Johnson, B E Bernstein, C Nussbaum, et al. 2008. “Model-Based Analysis of ChIP-Seq (MACS).” *Genome Biol.* 9 (9). <https://doi.org/10.1186/gb-2008-9-9-r137>.

</div>

<div id="ref-Zhu2010-zo" class="csl-entry">

Zhu, Lihua J, Claude Gazin, Nathan D Lawson, Hervé Pagès, Simon M Lin, David S Lapointe, and Michael R Green. 2010. “ChIPpeakAnno: A Bioconductor Package to Annotate <span class="nocase">ChIP-seq</span> and <span class="nocase">ChIP-chip</span> Data.” *BMC Bioinformatics* 11: 237. <https://doi.org/10.1186/1471-2105-11-237>.

</div>

</div>
