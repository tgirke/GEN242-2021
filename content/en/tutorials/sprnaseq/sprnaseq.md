---
title: "RNA-Seq Workflow Template" 
author: "Author: First Last"
date: "Last update: 29 May, 2021" 
output:
  BiocStyle::html_document:
    toc_float: true
    code_folding: show
package: systemPipeR
vignette: |
  %\VignetteIndexEntry{WF: RNA-Seq Workflow Template}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
fontsize: 14pt
bibliography: bibtex.bib
weight: 8
type: docs
---

<script src="/rmarkdown-libs/htmlwidgets/htmlwidgets.js"></script>
<script src="/rmarkdown-libs/jquery/jquery.min.js"></script>
<link href="/rmarkdown-libs/datatables-css/datatables-crosstalk.css" rel="stylesheet" />
<script src="/rmarkdown-libs/datatables-binding/datatables.js"></script>
<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.min.css" rel="stylesheet" />
<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.extra.css" rel="stylesheet" />
<script src="/rmarkdown-libs/dt-core/js/jquery.dataTables.min.js"></script>
<link href="/rmarkdown-libs/crosstalk/css/crosstalk.css" rel="stylesheet" />
<script src="/rmarkdown-libs/crosstalk/js/crosstalk.min.js"></script>
<script src="/rmarkdown-libs/htmlwidgets/htmlwidgets.js"></script>
<script src="/rmarkdown-libs/jquery/jquery.min.js"></script>
<link href="/rmarkdown-libs/datatables-css/datatables-crosstalk.css" rel="stylesheet" />
<script src="/rmarkdown-libs/datatables-binding/datatables.js"></script>
<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.min.css" rel="stylesheet" />
<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.extra.css" rel="stylesheet" />
<script src="/rmarkdown-libs/dt-core/js/jquery.dataTables.min.js"></script>
<link href="/rmarkdown-libs/crosstalk/css/crosstalk.css" rel="stylesheet" />
<script src="/rmarkdown-libs/crosstalk/js/crosstalk.min.js"></script>
<script src="/rmarkdown-libs/htmlwidgets/htmlwidgets.js"></script>
<script src="/rmarkdown-libs/jquery/jquery.min.js"></script>
<link href="/rmarkdown-libs/datatables-css/datatables-crosstalk.css" rel="stylesheet" />
<script src="/rmarkdown-libs/datatables-binding/datatables.js"></script>
<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.min.css" rel="stylesheet" />
<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.extra.css" rel="stylesheet" />
<script src="/rmarkdown-libs/dt-core/js/jquery.dataTables.min.js"></script>
<link href="/rmarkdown-libs/crosstalk/css/crosstalk.css" rel="stylesheet" />
<script src="/rmarkdown-libs/crosstalk/js/crosstalk.min.js"></script>
<style type="text/css">
pre code {
white-space: pre !important;
overflow-x: scroll !important;
word-break: keep-all !important;
word-wrap: initial !important;
}
</style>
<!--
Rscript -e "rmarkdown::render('sprnaseq.Rmd', c('BiocStyle::html_document'), clean=F); knitr::knit('sprnaseq.Rmd', tangle=TRUE)"
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
\[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/sprnaseq/sprnaseq.Rmd) \]    
\[ [.html](http://cluster.hpcc.ucr.edu/~ttest/rnaseq/systemPipeRNAseq.html) \]    
\[ [old version .Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/sprnaseq/back_sprnaseq_rmd) \]

</div>

## Introduction

This report describes the analysis of the RNA-Seq data set from
Howard et al (2013). The corresponding FASTQ files were downloaded from
GEO (Accession: [SRP010938](http://www.ncbi.nlm.nih.gov/sra/?term=SRP010938)).
This data set contains 18 paired-end (PE) read sets from *Arabidposis thaliana*.
The details about all download steps are provided [here](https://girke.bioinformatics.ucr.edu/GEN242/assignments/projects/project_data/).

Users want to provide here additional background information about the design of their
RNA-Seq project.

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
srun --x11 --partition=short --mem=8gb --cpus-per-task 4 --ntasks 1 --time 2:00:00 --pty bash -l
module unload R; module load R/4.0.3_gcc-8.3.0
```

2.  Load a workflow template with the `genWorkenvir` function. This can be done from the command-line or from within R.
    However, only one of the two options needs to be used.

From command-line

``` sh
$ Rscript -e "systemPipeRdata::genWorkenvir(workflow='rnaseq')"
$ cd rnaseq
```

From R

``` r
library(systemPipeRdata)
genWorkenvir(workflow = "rnaseq")
setwd("rnaseq")
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

To apply workflows to custom data, the user needs to modify the *`targets`* file and if
necessary update the corresponding parameter (*`.cwl`* and *`.yml`*) files.
A collection of pre-generated *`.cwl`* and *`.yml`* files are provided in the *`param/cwl`* subdirectory
of each workflow template. They are also viewable in the GitHub repository of *`systemPipeRdata`* ([see
here](https://github.com/tgirke/systemPipeRdata/tree/master/inst/extdata/param/cwl)).
For more information of the structure of the *targets* file, please consult the documentation
[here](http://www.bioconductor.org/packages/release/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html#25_structure_of_targets_file). More details about the new parameter files from systemPipeR can be found [here](http://www.bioconductor.org/packages/release/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html#26_structure_of_the_new_param_files_and_construct_sysargs2_container).

### Import custom functions

Custem functions for the challenge projects can be imported with the source
command from a local R script (here `custom_Fct.R`). Skip this step if such a
script is not available. Alternatively, these functions can be loaded from a
custom R package.

``` r
source("challengeProject_Fct.R")
```

### Experiment definition provided by `targets` file

The `targets` file defines all FASTQ files and sample comparisons of the analysis workflow.
If needed the tab separated (TSV) version of this file can be downloaded from [here](https://github.com/tgirke/GEN242/tree/main/content/en/assignments/Projects/targets_files)
and the corresponding Google Sheet is [here](https://docs.google.com/spreadsheets/d/1DTgTGlZZscSPjlHOGdJC8QK4vvimN1BORjXKXzd_cfA/edit#gid=472150521).

``` r
targetspath <- "targetsPE.txt"
targets <- read.delim(targetspath, comment.char = "#")
DT::datatable(targets, options = list(scrollX = TRUE, autoWidth = TRUE))
```

<div id="htmlwidget-1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-1">{"x":{"filter":"none","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"],["./data/SRR446027_1.fastq.gz","./data/SRR446028_1.fastq.gz","./data/SRR446029_1.fastq.gz","./data/SRR446030_1.fastq.gz","./data/SRR446031_1.fastq.gz","./data/SRR446032_1.fastq.gz","./data/SRR446033_1.fastq.gz","./data/SRR446034_1.fastq.gz","./data/SRR446035_1.fastq.gz","./data/SRR446036_1.fastq.gz","./data/SRR446037_1.fastq.gz","./data/SRR446038_1.fastq.gz","./data/SRR446039_1.fastq.gz","./data/SRR446040_1.fastq.gz","./data/SRR446041_1.fastq.gz","./data/SRR446042_1.fastq.gz","./data/SRR446043_1.fastq.gz","./data/SRR446044_1.fastq.gz"],["./data/SRR446027_2.fastq.gz","./data/SRR446028_2.fastq.gz","./data/SRR446029_2.fastq.gz","./data/SRR446030_2.fastq.gz","./data/SRR446031_2.fastq.gz","./data/SRR446032_2.fastq.gz","./data/SRR446033_2.fastq.gz","./data/SRR446034_2.fastq.gz","./data/SRR446035_2.fastq.gz","./data/SRR446036_2.fastq.gz","./data/SRR446037_2.fastq.gz","./data/SRR446038_2.fastq.gz","./data/SRR446039_2.fastq.gz","./data/SRR446040_2.fastq.gz","./data/SRR446041_2.fastq.gz","./data/SRR446042_2.fastq.gz","./data/SRR446043_2.fastq.gz","./data/SRR446044_2.fastq.gz"],["M1A","M1B","A1A","A1B","V1A","V1B","M6A","M6B","A6A","A6B","V6A","V6B","M12A","M12B","A12A","A12B","V12A","V12B"],["M1","M1","A1","A1","V1","V1","M6","M6","A6","A6","V6","V6","M12","M12","A12","A12","V12","V12"],["Mock.1h.A","Mock.1h.B","Avr.1h.A","Avr.1h.B","Vir.1h.A","Vir.1h.B","Mock.6h.A","Mock.6h.B","Avr.6h.A","Avr.6h.B","Vir.6h.A","Vir.6h.B","Mock.12h.A","Mock.12h.B","Avr.12h.A","Avr.12h.B","Vir.12h.A","Vir.12h.B"],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],["23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>FileName1<\/th>\n      <th>FileName2<\/th>\n      <th>SampleName<\/th>\n      <th>Factor<\/th>\n      <th>SampleLong<\/th>\n      <th>Experiment<\/th>\n      <th>Date<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"autoWidth":true,"columnDefs":[{"className":"dt-right","targets":6},{"orderable":false,"targets":0}],"order":[],"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

## Read preprocessing

### Read quality filtering and trimming

The function `preprocessReads` allows to apply predefined or custom
read preprocessing functions to all FASTQ files referenced in a
`SYSargs2` container, such as quality filtering or adapter trimming
routines. The paths to the resulting output FASTQ files are stored in the
`output` slot of the `SYSargs2` object. The following example performs adapter trimming with
the `trimLRPatterns` function from the `Biostrings` package.
After the trimming step a new targets file is generated (here
`targets_trim.txt`) containing the paths to the trimmed FASTQ files.
The new targets file can be used for the next workflow step with an updated
`SYSargs2` instance, *e.g.* running the NGS alignments using the
trimmed FASTQ files.

Construct *`SYSargs2`* object from *`cwl`* and *`yml`* param and *`targets`* files.

``` r
dir_path <- "param/cwl/preprocessReads/trim-pe"
trim <- loadWorkflow(targets = targetspath, wf_file = "trim-pe.cwl", 
    input_file = "trim-pe.yml", dir_path = dir_path)
trim <- renderWF(trim, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
    FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
trim
output(trim)[1:2]
```

``` r
preprocessReads(args = trim, Fct = "trimLRPatterns(Rpattern='GCCCGGGTAA', 
                subject=fq)", 
    batchsize = 1e+05, overwrite = TRUE, compress = TRUE)
writeTargetsout(x = trim, file = "targets_trim.txt", step = 1, 
    new_col = c("FileName1", "FileName2"), new_col_output_index = c(1, 
        2), overwrite = TRUE)
```

### FASTQ quality report

The following `seeFastq` and `seeFastqPlot` functions generate and plot a series of useful
quality statistics for a set of FASTQ files including per cycle quality box
plots, base proportions, base-level quality trends, relative k-mer
diversity, length and occurrence distribution of reads, number of reads
above quality cutoffs and mean quality distribution. The results are
written to a PDF file named `fastqReport.pdf`.

``` r
fqlist <- seeFastq(fastq = infile1(trim), batchsize = 10000, 
    klength = 8)
png("./results/fastqReport.png", height = 18 * 96, width = 4 * 
    96 * length(fqlist))
seeFastqPlot(fqlist)
dev.off()
```

![](../results/fastqReport.png)

<div align="center">

Figure 1: FASTQ quality report for 18 samples

</div>

</br>

## Alignments

### Read mapping with `HISAT2`

The following steps will demonstrate how to use the short read aligner `Hisat2`
(Kim, Langmead, and Salzberg 2015) in both interactive job submissions and batch submissions to
queuing systems of clusters using the *`systemPipeR's`* new CWL command-line interface.

Build `Hisat2` index.

``` r
dir_path <- "param/cwl/hisat2/hisat2-idx"
idx <- loadWorkflow(targets = NULL, wf_file = "hisat2-index.cwl", 
    input_file = "hisat2-index.yml", dir_path = dir_path)
idx <- renderWF(idx)
idx
cmdlist(idx)

## Run
runCommandline(idx, make_bam = FALSE)
```

The parameter settings of the aligner are defined in the `hisat2-mapping-se.cwl`
and `hisat2-mapping-se.yml` files. The following shows how to construct the
corresponding *SYSargs2* object, here *args*.

``` r
dir_path <- "param/cwl/hisat2/hisat2-pe"
args <- loadWorkflow(targets = targetspath, wf_file = "hisat2-mapping-pe.cwl", 
    input_file = "hisat2-mapping-pe.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
    FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
args
```

    ## Instance of 'SYSargs2':
    ##    Slot names/accessors: 
    ##       targets: 18 (M1A...V12B), targetsheader: 4 (lines)
    ##       modules: 1
    ##       wf: 0, clt: 1, yamlinput: 8 (components)
    ##       input: 18, output: 18
    ##       cmdlist: 18
    ##    WF Steps:
    ##       1. hisat2-mapping-pe (rendered: TRUE)

``` r
cmdlist(args)[1:2]
```

    ## $M1A
    ## $M1A$`hisat2-mapping-pe`
    ## [1] "hisat2 -S ./results/M1A.sam  -x ./data/tair10.fasta  -k 1  --min-intronlen 30  --max-intronlen 3000  -1 ./data/SRR446027_1.fastq.gz -2 ./data/SRR446027_2.fastq.gz --threads 4"
    ## 
    ## 
    ## $M1B
    ## $M1B$`hisat2-mapping-pe`
    ## [1] "hisat2 -S ./results/M1B.sam  -x ./data/tair10.fasta  -k 1  --min-intronlen 30  --max-intronlen 3000  -1 ./data/SRR446028_1.fastq.gz -2 ./data/SRR446028_2.fastq.gz --threads 4"

``` r
output(args)[1:2]
```

    ## $M1A
    ## $M1A$`hisat2-mapping-pe`
    ## [1] "./results/M1A.sam"
    ## 
    ## 
    ## $M1B
    ## $M1B$`hisat2-mapping-pe`
    ## [1] "./results/M1B.sam"

#### Interactive job submissions on a single machine

To simplify the short read alignment execution for the user, the command-line
can be run with the *`runCommandline`* function.
The execution will be on a single machine without submitting to a queuing system
of a computer cluster. This way, the input FASTQ files will be processed sequentially.
By default *`runCommandline`* auto detects SAM file outputs and converts them
to sorted and indexed BAM files, using internally the `Rsamtools` package.
Besides, *`runCommandline`* allows the user to create a dedicated
results folder for each workflow and a sub-folder for each sample
defined in the *targets* file. This includes all the output and log files for each
step. When these options are used, the output location will be updated by default
and can be assigned to the same object.

``` r
## Run single Machine
args <- runCommandline(args)
```

#### Parallelization on clusters

Alternatively, the computation can be greatly accelerated by processing many files
in parallel using several compute nodes of a cluster, where a scheduling/queuing
system is used for load balancing. For this the *`clusterRun`* function submits
the computing requests to the scheduler using the run specifications
defined by *`runCommandline`*.

To avoid over-subscription of CPU cores on the compute nodes, the value from
*`yamlinput(args)['thread']`* is passed on to the submission command, here *`ncpus`*
in the *`resources`* list object. The number of independent parallel cluster
processes is defined under the *`Njobs`* argument. The following example will run
18 processes in parallel using for each 4 CPU cores. If the resources available
on a cluster allow running all 18 processes at the same time then the shown sample
submission will utilize in total 72 CPU cores. Note, *`clusterRun`* can be used
with most queueing systems as it is based on utilities from the *`batchtools`*
package which supports the use of template files (*`*.tmpl`*) for defining the
run parameters of different schedulers. To run the following code, one needs to
have both a conf file (see *`.batchtools.conf.R`* samples [here](https://mllg.github.io/batchtools/))
and a template file (see *`*.tmpl`* samples [here](https://github.com/mllg/batchtools/tree/master/inst/templates))
for the queueing available on a system. The following example uses the sample
conf and template files for the Slurm scheduler provided by this package.

``` r
library(batchtools)
resources <- list(walltime = 120, ntasks = 1, ncpus = 4, memory = 1024)
reg <- clusterRun(args, FUN = runCommandline, more.args = list(args = args, 
    make_bam = TRUE, dir = FALSE), conffile = ".batchtools.conf.R", 
    template = "batchtools.slurm.tmpl", Njobs = 18, runid = "01", 
    resourceList = resources)
getStatus(reg = reg)
waitForJobs(reg = reg)
args <- output_update(args, dir = FALSE, replace = TRUE, extension = c(".sam", 
    ".bam"))  ## Updates the output(args) to the right location in the subfolders
output(args)
```

Check whether all BAM files have been created.

``` r
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
```

The following shows the alignment statistics for a sample file provided by the `systemPipeR` package.

``` r
read.table("results/alignStats.xls", header = TRUE)[1:4, ]
```

    ##   FileName Nreads2x   Nalign Perc_Aligned Nalign_Primary
    ## 1      M1A 33609678 32136300     95.61621       32136300
    ## 2      M1B 53002402 43620124     82.29839       43620124
    ## 3      A1A 50223496 48438407     96.44571       48438407
    ## 4      A1B 43650000 35549889     81.44304       35549889
    ##   Perc_Aligned_Primary
    ## 1             95.61621
    ## 2             82.29839
    ## 3             96.44571
    ## 4             81.44304

### Create symbolic links for viewing BAM files in IGV

The `symLink2bam` function creates symbolic links to view the BAM alignment files in a
genome browser such as IGV. The corresponding URLs are written to a file
with a path specified under `urlfile` in the `results` directory.

``` r
symLink2bam(sysargs = args, htmldir = c("~/.html/", "somedir/"), 
    urlbase = "http://cluster.hpcc.ucr.edu/~<username>/", urlfile = "./results/IGVurl.txt")
```

## Read quantification

### Read counting with `summarizeOverlaps` in parallel mode using multiple cores

Reads overlapping with annotation ranges of interest are counted for
each sample using the `summarizeOverlaps` function (Lawrence et al. 2013). The read counting is
preformed for exonic gene regions in a non-strand-specific manner while
ignoring overlaps among different genes. Subsequently, the expression
count values are normalized by *reads per kp per million mapped reads*
(RPKM). The raw read count table (`countDFeByg.xls`) and the corresponding
RPKM table (`rpkmDFeByg.xls`) are written to separate files in the directory of this project. Parallelization is achieved with the `BiocParallel` package, here using 8 CPU cores.

``` r
library("GenomicFeatures")
library(BiocParallel)
txdb <- makeTxDbFromGFF(file = "data/tair10.gff", format = "gff", 
    dataSource = "TAIR", organism = "Arabidopsis thaliana")
saveDb(txdb, file = "./data/tair10.sqlite")
txdb <- loadDb("./data/tair10.sqlite")
outpaths <- subsetWF(args, slot = "output", subset = 1, index = 1)
# (align <- readGAlignments(outpaths[1])) # Demonstrates how
# to read bam file into R
eByg <- exonsBy(txdb, by = c("gene"))
bfl <- BamFileList(outpaths, yieldSize = 50000, index = character())
multicoreParam <- MulticoreParam(workers = 4)
register(multicoreParam)
registered()
counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, 
    x, mode = "Union", ignore.strand = TRUE, inter.feature = FALSE, 
    singleEnd = FALSE))
countDFeByg <- sapply(seq(along = counteByg), function(x) assays(counteByg[[x]])$counts)
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]]))
colnames(countDFeByg) <- names(bfl)
rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts = x, 
    ranges = eByg))
write.table(countDFeByg, "results/countDFeByg.xls", col.names = NA, 
    quote = FALSE, sep = "\t")
write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names = NA, 
    quote = FALSE, sep = "\t")
```

Shows count table generated in previous step (`countDFeByg.xls`).
To avoid slowdowns of the load time of this page, ony 200 rows of the source
table are imported into the below `datatable` view .

``` r
countDF <- read.delim("results/countDFeByg.xls", row.names = 1, 
    check.names = FALSE)[1:200, ]
library(DT)
datatable(countDF, options = list(scrollX = TRUE, autoWidth = TRUE))
```

<div id="htmlwidget-2" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2">{"x":{"filter":"none","data":[["AT1G01010","AT1G01020","AT1G01030","AT1G01040","AT1G01046","AT1G01050","AT1G01060","AT1G01070","AT1G01073","AT1G01080","AT1G01090","AT1G01100","AT1G01110","AT1G01115","AT1G01120","AT1G01130","AT1G01140","AT1G01150","AT1G01160","AT1G01170","AT1G01180","AT1G01183","AT1G01190","AT1G01200","AT1G01210","AT1G01220","AT1G01225","AT1G01230","AT1G01240","AT1G01250","AT1G01260","AT1G01270","AT1G01280","AT1G01290","AT1G01300","AT1G01305","AT1G01310","AT1G01320","AT1G01340","AT1G01350","AT1G01355","AT1G01360","AT1G01370","AT1G01380","AT1G01390","AT1G01400","AT1G01410","AT1G01420","AT1G01430","AT1G01440","AT1G01448","AT1G01450","AT1G01453","AT1G01460","AT1G01470","AT1G01471","AT1G01480","AT1G01490","AT1G01500","AT1G01510","AT1G01520","AT1G01530","AT1G01540","AT1G01550","AT1G01560","AT1G01570","AT1G01580","AT1G01590","AT1G01600","AT1G01610","AT1G01620","AT1G01630","AT1G01640","AT1G01650","AT1G01660","AT1G01670","AT1G01680","AT1G01690","AT1G01695","AT1G01700","AT1G01710","AT1G01720","AT1G01725","AT1G01730","AT1G01740","AT1G01750","AT1G01760","AT1G01770","AT1G01780","AT1G01790","AT1G01800","AT1G01810","AT1G01820","AT1G01830","AT1G01840","AT1G01860","AT1G01870","AT1G01880","AT1G01890","AT1G01900","AT1G01910","AT1G01920","AT1G01930","AT1G01940","AT1G01950","AT1G01960","AT1G01970","AT1G01980","AT1G01990","AT1G02000","AT1G02010","AT1G02020","AT1G02030","AT1G02040","AT1G02050","AT1G02060","AT1G02065","AT1G02070","AT1G02074","AT1G02080","AT1G02090","AT1G02100","AT1G02110","AT1G02120","AT1G02130","AT1G02136","AT1G02140","AT1G02145","AT1G02150","AT1G02160","AT1G02170","AT1G02180","AT1G02190","AT1G02205","AT1G02210","AT1G02220","AT1G02230","AT1G02250","AT1G02260","AT1G02270","AT1G02280","AT1G02290","AT1G02300","AT1G02305","AT1G02310","AT1G02320","AT1G02330","AT1G02335","AT1G02340","AT1G02350","AT1G02360","AT1G02370","AT1G02380","AT1G02390","AT1G02391","AT1G02400","AT1G02405","AT1G02410","AT1G02420","AT1G02430","AT1G02440","AT1G02450","AT1G02460","AT1G02470","AT1G02475","AT1G02480","AT1G02490","AT1G02500","AT1G02510","AT1G02520","AT1G02530","AT1G02540","AT1G02550","AT1G02560","AT1G02570","AT1G02575","AT1G02580","AT1G02590","AT1G02600","AT1G02610","AT1G02620","AT1G02630","AT1G02640","AT1G02650","AT1G02660","AT1G02670","AT1G02680","AT1G02681","AT1G02690","AT1G02700","AT1G02710","AT1G02720","AT1G02730","AT1G02740","AT1G02750","AT1G02760","AT1G02770","AT1G02780","AT1G02790","AT1G02800"],[405,155,185,1342,30,283,139,5,0,546,1732,1263,16,0,4059,16,2494,0,393,417,100,0,62,11,77,236,16,562,515,19,396,0,0,104,845,0,3,1209,209,305,5,165,9,1,23,0,0,153,306,212,56,30,1,0,15440,1,102,990,451,243,1,0,929,3116,105,21,28,7,22,1287,1195,411,15,1510,20,29,14,8,2,27,426,1304,189,378,42,7,19,73,34,2483,604,6,567,48,39,102,0,17,0,57,631,130,258,101,154,819,243,0,268,330,101,187,1,1,58,43,0,0,0,1227,473,368,173,1280,390,0,1199,116,1030,314,995,81,2,387,0,97,12,0,95,587,229,25,123,442,2,0,349,43,30,15,288,36,13,87,0,2419,35,225,51,30,10,5,205,1,337,0,0,3825,0,22,3,0,0,2644,0,0,1,0,0,151,6,0,131,12,1247,26,98,1,24,1,7,76,43,70,80,0,6,4107,1,0],[577,327,202,1732,29,367,698,2,0,1025,2934,1246,20,0,6526,25,5619,0,514,319,74,0,59,25,131,312,16,782,326,18,894,0,0,203,2940,0,4,3405,334,585,3,401,6,2,106,0,0,392,682,215,95,72,4,7,17055,0,240,2455,668,570,68,0,2339,3329,176,24,100,1,145,4938,1642,679,20,3121,29,20,35,25,4,34,519,3186,306,384,134,25,35,149,72,3641,697,26,657,163,26,84,0,47,1,232,1443,292,365,113,328,1638,355,2,838,608,237,478,11,9,103,103,0,0,0,2506,745,592,272,1980,418,0,1184,220,1139,255,948,36,5,441,0,166,11,2,148,782,459,55,121,461,31,0,616,68,54,16,322,86,18,123,0,4715,43,343,177,92,10,14,351,7,549,0,0,5851,0,40,11,0,0,2686,0,0,1,0,0,94,1,0,130,23,2081,34,107,0,100,1,32,160,68,107,118,0,1,5334,0,4],[625,235,281,1737,30,435,145,7,0,630,2284,1609,34,0,5461,35,4252,0,614,507,182,0,20,12,141,342,20,1095,688,34,636,0,0,164,2141,0,2,1144,667,644,6,310,12,3,36,0,0,185,424,301,106,62,1,2,28423,0,172,1588,854,302,1,1,1649,6380,319,34,64,2,48,2413,1387,747,28,2892,18,28,44,14,2,38,725,2014,293,947,105,11,29,94,41,3901,958,10,1214,117,40,100,0,15,0,99,1107,265,422,224,243,1311,247,1,421,602,188,297,2,2,71,48,0,2,0,1890,901,840,179,2567,673,0,1920,217,1218,429,2195,104,8,418,0,134,17,4,169,1152,237,35,257,861,34,0,678,55,45,20,370,54,11,195,0,5360,82,362,106,56,16,9,488,5,452,0,0,11353,2,29,6,0,0,4789,1,0,1,0,0,238,6,0,285,24,2475,41,207,0,44,4,14,102,35,120,136,0,7,5843,0,2],[329,243,100,1200,17,306,148,1,0,882,2389,904,19,0,3851,39,2804,0,333,189,68,0,40,36,95,228,24,611,394,39,707,0,0,197,2585,0,4,2326,326,460,0,239,11,3,115,0,0,226,455,182,68,38,0,0,11281,0,272,2214,550,398,26,0,1519,2903,210,6,44,1,47,2825,1241,627,18,2133,25,31,20,16,27,8,432,945,65,293,76,16,20,119,73,2949,639,7,714,226,10,39,0,11,0,287,1160,184,291,139,256,1254,210,0,713,437,108,297,3,0,54,69,0,3,0,1835,718,552,219,2130,380,0,1019,120,745,224,856,63,21,232,0,163,25,0,224,585,299,49,89,372,35,0,525,80,40,1,537,68,14,88,0,2662,28,346,118,78,7,90,204,25,429,0,0,6069,0,27,3,0,0,2535,0,0,0,0,0,80,0,0,72,13,1136,27,79,0,62,7,19,127,67,93,71,0,2,4596,1,6],[903,282,225,1605,28,356,52,16,0,759,1850,1039,12,0,4089,18,3356,0,651,625,154,1,13,2,123,343,28,1005,771,27,545,0,0,146,1238,0,1,1351,565,472,11,359,14,0,19,0,0,127,432,259,153,75,1,2,23646,2,347,1689,628,234,2,0,1478,5078,204,21,56,2,33,1890,2392,662,49,2389,24,48,28,15,2,27,618,2113,355,794,119,15,22,141,33,3396,842,7,1367,117,27,107,0,26,0,66,1245,175,373,146,183,1590,271,3,447,505,166,212,1,0,64,51,1,2,0,2247,619,757,123,2198,493,0,1650,193,1012,338,1409,185,3,557,0,110,19,1,197,891,291,31,289,770,40,0,582,50,50,14,429,56,43,171,0,4698,30,243,119,82,8,12,485,11,528,0,1,6339,1,29,3,0,0,4141,0,1,1,0,0,328,25,1,242,19,2030,36,189,0,42,2,2,89,35,102,121,0,11,4173,0,3],[517,269,329,1580,41,663,49,6,0,1128,4010,1980,54,0,6599,37,4687,0,563,709,210,0,54,18,86,370,31,1156,434,19,1125,0,3,234,3599,0,1,1849,668,653,8,205,23,1,32,0,0,198,640,337,212,135,2,4,27160,0,233,1243,881,478,3,0,2338,6861,617,43,70,0,28,1475,1904,726,52,3834,19,40,42,29,2,59,1011,1986,330,647,92,19,45,160,78,5021,1093,4,1509,156,23,83,0,16,0,183,1231,184,391,209,434,1852,367,0,526,568,208,468,7,2,97,126,0,1,0,3198,746,866,315,2033,669,0,1145,261,1661,347,2403,85,10,458,0,112,4,5,417,1534,531,38,186,1030,43,0,503,81,22,22,1383,36,19,297,0,10035,69,506,116,69,12,8,208,8,612,0,0,10358,0,37,1,0,0,6639,0,1,1,2,2,54,7,2,193,15,2832,127,132,0,67,6,4,200,90,156,137,0,3,6188,0,3],[597,349,169,2032,49,1210,12,26,0,1091,4073,4821,13,0,2965,14,3156,0,1492,1190,394,0,59,13,300,710,25,703,642,45,250,0,0,233,1014,0,2,12417,434,508,11,481,12,8,82,0,0,337,1165,462,318,131,0,0,5538,0,44,2214,1265,1000,0,1,2217,2118,119,54,8,8,71,2827,5632,1222,49,2401,27,97,175,30,5,41,844,1123,284,852,100,23,63,581,231,2885,1759,15,2107,337,156,381,0,12,0,754,1674,455,366,268,515,2506,447,0,802,463,203,765,4,2,303,135,0,0,2,2517,899,846,732,3793,976,0,2915,366,3192,1023,1951,127,4,1450,0,107,53,0,648,376,856,101,637,3444,6,0,488,34,11,15,211,188,38,48,0,483,13,307,98,49,1,47,399,17,1116,0,0,13817,1,104,6,0,0,5769,1,0,0,0,0,304,83,0,195,16,365,12,258,0,72,8,25,132,80,169,289,0,15,16024,1,0],[339,294,69,1399,26,1398,10,31,0,1716,4449,4861,28,0,1697,25,4793,0,1451,1201,418,0,148,27,250,866,42,769,542,87,321,0,0,274,1347,0,5,12836,402,467,1,479,23,4,102,0,0,259,1524,555,264,142,0,0,9000,0,123,2109,1009,1113,3,3,2099,2544,414,31,11,30,37,1979,5189,1275,65,2469,38,104,313,40,2,34,959,1000,230,608,120,27,72,727,213,4425,1853,15,2160,365,144,384,0,21,0,472,2244,394,483,349,483,3369,544,1,695,507,187,655,2,6,153,132,0,1,0,2698,1024,865,645,3255,964,0,3402,342,3014,978,1706,157,2,1065,0,178,111,1,440,470,870,135,809,2760,24,0,489,64,9,13,267,172,43,43,0,456,13,394,121,65,8,331,218,17,1428,0,1,13827,1,104,6,0,0,6038,0,0,1,0,0,865,41,1,286,19,1477,21,275,0,82,8,25,179,191,168,365,0,6,15871,2,2],[1185,294,63,2038,63,909,10,64,0,602,2469,5652,38,0,915,19,1601,1,1370,850,465,0,21,15,368,547,23,767,461,94,438,0,0,291,959,0,4,5363,2667,517,10,579,42,17,22,0,0,263,400,844,328,141,6,1,4742,0,215,1293,761,807,0,11,1411,4843,2017,42,13,14,50,1071,1628,1075,124,3004,82,146,13221,226,13,59,767,1838,338,1608,99,7,39,277,117,2245,1657,12,2147,384,109,254,0,24,0,344,1437,467,315,484,526,2546,232,0,806,343,130,686,2,6,340,168,2,2,0,2875,1190,1239,678,5816,919,0,3294,331,2339,782,3263,134,1,382,0,746,157,1,789,3303,582,143,370,2581,42,1,590,68,8,14,537,270,32,781,0,831,26,439,118,114,32,1244,474,78,1007,0,0,35482,2,684,4,1,19,5124,0,0,7,0,0,51,5,0,61,5,443,19,308,0,116,4,24,93,133,475,691,0,21,17511,0,1],[861,297,190,1764,40,1150,4,24,0,924,3411,3829,31,0,901,8,3028,0,1168,914,375,1,73,19,291,610,12,679,546,20,453,0,0,248,1002,0,12,8487,1295,544,2,335,38,5,22,0,0,173,658,464,294,159,1,0,5118,0,316,2069,697,897,1,5,1817,2780,767,35,2,5,24,1386,4496,738,83,2235,55,81,2635,91,6,41,642,1519,310,894,88,25,78,589,179,2391,1552,11,1813,280,90,247,0,14,0,539,1571,367,367,323,521,2998,309,0,803,333,173,473,2,7,233,178,0,2,0,2358,837,999,575,4360,836,0,2779,273,2272,685,2064,131,2,701,0,685,177,3,537,961,792,80,427,2259,39,0,495,79,6,12,324,118,21,246,0,612,17,337,112,64,19,530,407,38,1123,0,1,27669,1,218,7,0,4,4146,0,0,1,0,0,153,17,1,224,13,510,34,236,0,87,10,30,156,110,243,387,0,3,12786,8,0],[2141,529,290,4308,71,1674,3,71,0,971,3969,7384,47,0,2748,30,5419,1,2103,1145,511,0,47,30,516,953,26,1028,941,104,683,0,0,354,2190,0,4,16385,1892,768,13,785,29,15,30,0,0,485,1151,846,493,281,3,1,11323,0,288,2592,1011,1379,0,6,2204,5832,1049,56,17,45,85,2982,4825,1513,220,5593,80,162,4248,201,13,97,1182,2218,470,1967,223,35,107,920,369,3022,2078,28,3230,587,142,509,0,22,0,1006,1916,738,552,550,870,4648,471,0,1460,571,231,880,2,8,565,346,0,1,0,4836,1351,1924,1158,8165,1794,0,4991,560,4036,962,3048,170,12,1835,0,441,142,3,803,1368,1236,204,892,5357,29,0,731,71,20,12,636,358,45,194,0,954,27,503,156,110,16,263,980,76,1445,0,0,44059,1,719,5,0,2,6821,1,0,12,0,0,149,9,0,143,9,612,21,387,1,151,3,20,166,175,445,679,0,13,26651,0,2],[1370,575,519,3800,93,2401,14,47,0,1841,5859,6912,84,0,2380,25,6432,1,2110,1603,415,0,59,28,585,931,32,954,769,37,876,0,1,443,2415,0,15,20227,530,693,17,600,20,14,23,0,1,264,1460,470,596,341,1,3,8181,0,152,3200,1027,1773,0,1,3000,3506,236,36,12,26,85,3948,9477,1164,154,4766,52,111,549,57,5,69,1206,2792,611,1335,219,35,148,1514,477,3720,2197,18,2968,437,178,504,0,23,0,985,2425,710,594,410,679,5112,649,0,1633,611,285,851,5,12,520,392,0,4,0,4002,1203,1706,1381,4740,1366,0,4650,588,4974,1118,1735,228,1,1719,0,263,50,0,643,529,2115,174,1023,5301,33,1,553,111,23,14,129,292,51,50,0,448,26,324,114,52,2,95,716,32,1732,0,0,32119,0,260,8,0,0,6977,1,2,1,0,0,282,39,0,403,10,604,21,409,0,163,21,33,249,128,218,563,0,11,18362,0,1],[310,104,79,983,37,1013,321,14,0,639,2029,2185,21,0,950,1,1433,0,771,500,243,0,18,10,86,433,10,333,292,59,172,0,0,152,1004,0,0,9085,344,145,8,279,5,0,47,0,0,230,751,260,205,122,1,0,1426,0,63,975,264,690,0,15,1208,1743,195,25,8,32,17,1697,3967,626,71,1235,26,45,686,21,1,38,500,623,115,477,53,4,31,280,151,2431,808,3,732,440,75,200,0,17,0,243,546,210,190,82,301,1508,268,0,579,226,123,447,12,0,249,154,0,3,0,2250,337,665,588,1464,383,0,831,214,1676,229,900,61,2,470,0,116,26,0,381,250,521,70,390,2247,6,0,193,45,11,10,93,57,5,61,0,77,3,195,45,36,1,54,131,58,630,0,0,12780,0,91,0,0,0,3297,0,0,10,0,0,44,4,0,83,9,106,7,105,0,55,7,59,93,69,97,143,0,1,6524,0,6],[562,359,136,2724,43,3321,770,125,0,2603,6336,5860,109,0,2097,40,4947,3,2656,2554,607,0,194,35,417,931,29,830,680,119,538,0,0,454,959,0,14,21880,327,635,19,987,19,11,141,0,0,424,2018,732,908,543,0,0,3016,0,56,3829,664,2012,2,10,4135,3390,163,47,15,33,28,2729,15812,1285,185,2019,79,140,193,36,6,78,1306,1125,704,998,174,42,143,1115,380,6592,2814,13,2639,701,228,641,0,39,0,284,1356,538,620,307,775,3934,788,0,1440,546,475,974,4,26,874,516,0,4,1,3694,891,2043,2586,6195,1168,0,4203,586,5586,1177,1880,242,1,352,0,109,41,2,1151,398,2854,180,1208,6661,32,1,488,94,37,23,124,162,29,58,0,157,18,471,111,69,8,48,168,45,1581,0,0,21376,0,93,3,0,0,7496,1,1,4,0,0,283,25,2,397,23,243,18,335,2,96,44,98,255,172,184,380,0,11,16012,5,3],[1269,326,91,3359,86,1836,1207,161,0,887,2774,6327,32,0,1179,31,3638,0,2312,1475,512,0,20,22,314,711,18,827,933,318,513,0,0,400,956,0,12,16928,3228,648,24,1018,79,15,74,0,0,483,606,975,414,210,14,0,4146,0,459,2787,544,996,25,40,1892,5389,1238,77,22,26,80,1354,2461,1182,203,4382,132,197,7587,99,19,95,1103,2102,623,1829,149,16,61,750,323,3257,2369,23,2229,1039,146,321,0,18,0,393,1533,477,463,439,622,3349,368,0,1129,342,237,937,8,17,624,323,1,1,0,4097,1371,2093,1103,7531,1486,0,4323,545,2702,990,3080,183,1,871,0,1091,106,0,1208,1009,699,199,1040,5809,26,0,606,32,15,14,316,156,39,326,0,208,77,621,118,130,16,558,269,644,1655,0,0,44136,2,488,5,0,9,6928,0,0,19,0,0,95,2,1,61,3,462,26,485,0,73,12,22,127,71,414,800,0,23,20836,0,2],[1424,372,164,2245,40,1732,651,141,0,984,3166,2761,52,0,665,19,3642,3,1743,1350,605,0,48,5,317,521,13,584,992,31,626,0,1,291,615,0,9,13377,2505,609,6,604,45,3,46,0,0,265,619,780,440,343,2,1,2716,0,849,3391,338,863,2,29,1879,4161,1130,11,6,11,20,758,5843,667,242,2619,85,119,6642,69,5,46,860,1629,621,930,118,34,95,1039,340,2671,1879,4,2061,443,80,214,0,21,0,197,1485,409,405,274,417,2845,401,0,922,279,201,569,2,13,398,363,0,1,0,3108,788,1993,857,6542,940,0,2842,371,2495,541,1857,151,0,108,0,1125,225,2,856,879,1566,119,701,4391,49,0,460,82,23,11,278,103,22,145,0,155,38,453,156,111,18,418,324,328,1613,0,0,25994,0,316,4,0,8,4379,3,1,31,0,0,202,12,0,176,1,862,25,311,0,71,20,43,100,91,287,731,0,4,9236,0,0],[1017,143,280,1076,37,997,254,35,0,406,1864,2513,21,0,1283,5,1318,0,694,483,339,0,13,11,95,320,12,354,214,87,219,0,0,116,877,0,4,6539,910,213,4,207,11,3,38,0,0,200,423,400,195,104,2,0,3821,0,262,1140,271,560,5,5,991,2392,764,36,12,14,44,1964,3898,482,74,2012,32,53,2174,35,3,35,426,1745,165,510,60,6,44,312,144,1927,812,4,785,302,39,127,0,3,0,572,504,192,174,112,197,1446,201,0,506,270,121,392,2,0,192,172,0,0,0,1991,435,757,582,1801,419,0,1132,220,1437,164,1305,69,6,1180,0,208,36,1,361,496,675,64,265,1816,2,0,197,31,15,3,80,49,12,132,0,182,62,211,40,35,1,189,344,84,607,0,0,15287,0,120,3,0,4,3287,0,0,2,0,0,31,4,1,39,1,188,3,137,0,38,8,25,79,31,91,271,0,3,5998,0,3],[1571,366,1006,2326,54,2296,386,46,0,940,3444,4968,96,0,1683,17,2811,1,1556,1456,295,0,30,26,477,420,10,487,440,75,590,0,0,234,1549,1,6,10184,739,574,4,483,10,12,30,0,0,221,478,425,299,182,0,2,11223,1,341,3426,412,986,18,2,1519,2954,467,14,8,23,41,2446,14560,470,130,3008,75,92,1317,38,1,70,692,6477,709,1003,110,20,93,778,259,2095,1575,10,1649,377,101,220,0,16,0,434,958,418,344,227,316,1867,400,0,1030,277,212,471,2,14,420,272,0,2,0,2323,793,1495,1417,3456,867,0,3361,405,3729,619,1603,217,5,1165,0,344,229,1,393,378,2219,82,473,3445,8,0,382,78,26,7,46,88,29,33,0,271,140,319,88,45,10,160,840,57,823,0,0,21832,1,98,1,0,1,4802,0,0,2,0,0,105,5,0,90,3,258,13,390,1,66,13,34,100,76,126,568,0,5,10817,0,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>M1A<\/th>\n      <th>M1B<\/th>\n      <th>A1A<\/th>\n      <th>A1B<\/th>\n      <th>V1A<\/th>\n      <th>V1B<\/th>\n      <th>M6A<\/th>\n      <th>M6B<\/th>\n      <th>A6A<\/th>\n      <th>A6B<\/th>\n      <th>V6A<\/th>\n      <th>V6B<\/th>\n      <th>M12A<\/th>\n      <th>M12B<\/th>\n      <th>A12A<\/th>\n      <th>A12B<\/th>\n      <th>V12A<\/th>\n      <th>V12B<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"autoWidth":true,"columnDefs":[{"className":"dt-right","targets":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]},{"orderable":false,"targets":0}],"order":[],"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

A data slice of RPKM table (`rpkmDFeByg.xls`) is shown here.

``` r
read.delim("results/rpkmDFeByg.xls", row.names = 1, check.names = FALSE)[1:4, 
    1:4]
```

    ##                 M1A       M1B       A1A       A1B
    ## AT1G01010 15.552350 15.855557 15.515099 11.482534
    ## AT1G01020  5.663586  8.550121  5.550872  8.069877
    ## AT1G01030  6.294920  4.918521  6.180994  3.092568
    ## AT1G01040 13.909390 12.846007 11.638283 11.304143

Note, for most statistical differential expression or abundance analysis
methods, such as `edgeR` or `DESeq2`, the raw count values should be used as input. The
usage of RPKM values should be restricted to specialty applications
required by some users, *e.g.* manually comparing the expression levels
among different genes or features.

### Sample-wise correlation analysis

The following computes the sample-wise Spearman correlation coefficients from
the `rlog` transformed expression values generated with the `DESeq2` package. After
transformation to a distance matrix, hierarchical clustering is performed with
the `hclust` function and the result is plotted as a dendrogram
(also see file `sample_tree.pdf`).

``` r
library(DESeq2, quietly = TRUE)
library(ape, warn.conflicts = FALSE)
countDF <- as.matrix(read.table("./results/countDFeByg.xls"))
colData <- data.frame(row.names = targets.as.df(targets(args))$SampleName, 
    condition = targets.as.df(targets(args))$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, 
    design = ~condition)
d <- cor(assay(rlog(dds)), method = "spearman")
hc <- hclust(dist(1 - d))
png("results/sample_tree.png")
plot.phylo(as.phylo(hc), type = "p", edge.col = "blue", edge.width = 2, 
    show.node.label = TRUE, no.margin = TRUE)
dev.off()
```

![](../results/sample_tree.png)

<div align="center">

Figure 2: Correlation dendrogram of samples

</div>

</br>

## Analysis of DEGs

The analysis of differentially expressed genes (DEGs) is performed with
the glm method of the `edgeR` package (Robinson, McCarthy, and Smyth 2010). The sample
comparisons used by this analysis are defined in the header lines of the
`targets.txt` file starting with `<CMP>`.

### Run `edgeR`

``` r
library(edgeR)
countDF <- read.delim("results/countDFeByg.xls", row.names = 1, 
    check.names = FALSE)
targets <- read.delim("targetsPE.txt", comment = "#")
cmp <- readComp(file = "targetsPE.txt", format = "matrix", delim = "-")
edgeDF <- run_edgeR(countDF = countDF, targets = targets, cmp = cmp[[1]], 
    independent = FALSE, mdsplot = "")
```

Add gene descriptions

``` r
library("biomaRt")
m <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "plants.ensembl.org")
desc <- getBM(attributes = c("tair_locus", "description"), mart = m)
desc <- desc[!duplicated(desc[, 1]), ]
descv <- as.character(desc[, 2])
names(descv) <- as.character(desc[, 1])
edgeDF <- data.frame(edgeDF, Desc = descv[rownames(edgeDF)], 
    check.names = FALSE)
write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote = FALSE, 
    sep = "\t", col.names = NA)
```

### Plot DEG results

Filter and plot DEG results for up and down regulated genes. The
definition of *up* and *down* is given in the corresponding help
file. To open it, type `?filterDEGs` in the R console.

``` r
edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names = 1, 
    check.names = FALSE)
png("results/DEGcounts.png")
DEG_list <- filterDEGs(degDF = edgeDF, filter = c(Fold = 2, FDR = 20))
dev.off()
write.table(DEG_list$Summary, "./results/DEGcounts.xls", quote = FALSE, 
    sep = "\t", row.names = FALSE)
```

![](../results/DEGcounts.png)

<div align="center">

Figure 3: Up and down regulated DEGs with FDR of 1%

</div>

</br>

### Venn diagrams of DEG sets

The `overLapper` function can compute Venn intersects for large numbers of sample
sets (up to 20 or more) and plots 2-5 way Venn diagrams. A useful
feature is the possibility to combine the counts from several Venn
comparisons with the same number of sample sets in a single Venn diagram
(here for 4 up and down DEG sets).

``` r
vennsetup <- overLapper(DEG_list$Up[6:9], type = "vennsets")
vennsetdown <- overLapper(DEG_list$Down[6:9], type = "vennsets")
png("results/vennplot.png")
vennPlot(list(vennsetup, vennsetdown), mymain = "", mysub = "", 
    colmode = 2, ccol = c("blue", "red"))
dev.off()
```

![](../results/vennplot.png)

<div align="center">

Figure 4: Venn Diagram for 4 Up and Down DEG Sets

</div>

</br>

## GO term enrichment analysis

### Obtain gene-to-GO mappings

The following shows how to obtain gene-to-GO mappings from `biomaRt` (here for *A.
thaliana*) and how to organize them for the downstream GO term
enrichment analysis. Alternatively, the gene-to-GO mappings can be
obtained for many organisms from Bioconductor’s `*.db` genome annotation
packages or GO annotation files provided by various genome databases.
For each annotation this relatively slow preprocessing step needs to be
performed only once. Subsequently, the preprocessed data can be loaded
with the `load` function as shown in the next subsection.

``` r
library("biomaRt")
listMarts()  # To choose BioMart database
listMarts(host = "plants.ensembl.org")
m <- useMart("plants_mart", host = "plants.ensembl.org")
listDatasets(m)
m <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "plants.ensembl.org")
listAttributes(m)  # Choose data types you want to download
go <- getBM(attributes = c("go_id", "tair_locus", "namespace_1003"), 
    mart = m)
go <- go[go[, 3] != "", ]
go[, 3] <- as.character(go[, 3])
go[go[, 3] == "molecular_function", 3] <- "F"
go[go[, 3] == "biological_process", 3] <- "P"
go[go[, 3] == "cellular_component", 3] <- "C"
go[1:4, ]
dir.create("./data/GO")
write.table(go, "data/GO/GOannotationsBiomart_mod.txt", quote = FALSE, 
    row.names = FALSE, col.names = FALSE, sep = "\t")
catdb <- makeCATdb(myfile = "data/GO/GOannotationsBiomart_mod.txt", 
    lib = NULL, org = "", colno = c(1, 2, 3), idconv = NULL)
save(catdb, file = "data/GO/catdb.RData")
```

### Batch GO term enrichment analysis

Apply the enrichment analysis to the DEG sets obtained the above differential
expression analysis. Note, in the following example the `FDR` filter is set
here to an unreasonably high value, simply because of the small size of the toy
data set used in this vignette. Batch enrichment analysis of many gene sets is
performed with the function. When `method=all`, it returns all GO terms passing
the p-value cutoff specified under the `cutoff` arguments. When `method=slim`,
it returns only the GO terms specified under the `myslimv` argument. The given
example shows how a GO slim vector for a specific organism can be obtained from
BioMart.

``` r
library("biomaRt")
load("data/GO/catdb.RData")
DEG_list <- filterDEGs(degDF = edgeDF, filter = c(Fold = 2, FDR = 50), 
    plot = FALSE)
up_down <- DEG_list$UporDown
names(up_down) <- paste(names(up_down), "_up_down", sep = "")
up <- DEG_list$Up
names(up) <- paste(names(up), "_up", sep = "")
down <- DEG_list$Down
names(down) <- paste(names(down), "_down", sep = "")
DEGlist <- c(up_down, up, down)
DEGlist <- DEGlist[sapply(DEGlist, length) > 0]
BatchResult <- GOCluster_Report(catdb = catdb, setlist = DEGlist, 
    method = "all", id_type = "gene", CLSZ = 2, cutoff = 0.9, 
    gocats = c("MF", "BP", "CC"), recordSpecGO = NULL)
write.table(BatchResult, "results/GOBatchAll.xls", row.names = FALSE, 
    quote = FALSE, sep = "\t")
library("biomaRt")
m <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "plants.ensembl.org")
goslimvec <- as.character(getBM(attributes = c("goslim_goa_accession"), 
    mart = m)[, 1])
BatchResultslim <- GOCluster_Report(catdb = catdb, setlist = DEGlist, 
    method = "slim", id_type = "gene", myslimv = goslimvec, CLSZ = 10, 
    cutoff = 0.01, gocats = c("MF", "BP", "CC"), recordSpecGO = NULL)
write.table(BatchResultslim, "results/GOBatchSlim.xls", row.names = FALSE, 
    quote = FALSE, sep = "\t")
```

Shows GO term enrichment results from previous step. The last gene identifier column (10)
of this table has been excluded in this viewing instance to minimze the complexity of the
result.
To avoid slowdowns of the load time of this page, ony 200 rows of the source
table are imported into the below `datatable` view .

``` r
BatchResult <- read.delim("results/GOBatchAll.xls")[1:199, ]
library(DT)
datatable(BatchResult[, -10], options = list(scrollX = TRUE, 
    autoWidth = TRUE))
```

<div id="htmlwidget-2" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2">{"x":{"filter":"none","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199"],["M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-A1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down","M1-V1_up_down"],[26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164,1164],["GO:0050291","GO:0004345","GO:0050664","GO:0052593","GO:0052594","GO:0052595","GO:0052596","GO:0052793","GO:0008131","GO:0016018","GO:0004144","GO:0008144","GO:0016411","GO:0009881","GO:0016641","GO:0048038","GO:0016638","GO:0008171","GO:0016747","GO:0004674","GO:0016301","GO:0005515","GO:0003755","GO:0016859","GO:0030247","GO:0016746","GO:0005488","GO:0016410","GO:0050661","GO:0042277","no_annot_MF","GO:0007623","GO:0048511","GO:0042752","GO:0009637","GO:0048580","GO:0009409","GO:0019438","GO:0032091","GO:0009646","GO:0010114","GO:0051100","GO:0006355","GO:1903506","GO:2001141","GO:2000028","GO:2000026","GO:0009662","GO:0043393","GO:2001007","GO:0051252","GO:0009737","GO:0097305","GO:0019219","GO:0048582","GO:0006351","GO:0097659","GO:0009266","GO:0051239","GO:0032774","GO:2000112","GO:0010556","GO:0044249","GO:1901362","GO:0031326","GO:0009889","GO:0009806","GO:0009807","GO:0051240","GO:2000241","GO:0009058","GO:0042759","GO:0048573","GO:0080167","GO:0051094","GO:2000113","GO:0010558","GO:0010677","GO:0010468","GO:0034654","GO:0050793","GO:0031327","GO:0009642","GO:0009648","GO:0033993","GO:0009890","GO:0042221","GO:0051171","GO:0042793","GO:0045912","GO:0018130","GO:0080090","GO:0009051","GO:0051098","GO:2001006","GO:0043207","GO:0051707","GO:0044271","GO:0032950","GO:0032951","GO:0045995","GO:0010033","GO:0009607","GO:0046513","GO:0050896","GO:1901576","GO:0044419","GO:0031323","GO:0002238","GO:0060255","GO:0010030","GO:0031324","GO:0009628","GO:0010228","GO:0009639","GO:0001676","GO:0010962","GO:0044237","GO:0006098","GO:0019432","GO:1901700","GO:0006071","GO:0006672","GO:0019400","GO:0006740","GO:0032885","GO:0046460","GO:0046463","GO:0019222","GO:0009682","GO:0006725","GO:0006641","GO:0051156","GO:0006638","GO:0006639","GO:0006739","GO:0006950","GO:0016070","GO:0071704","GO:0048608","GO:0061458","GO:0009911","GO:0030148","GO:0032881","GO:0009416","GO:0009605","GO:0010605","GO:0009314","GO:0045892","GO:1902679","GO:1903507","GO:0034645","no_annot_BP","GO:0000427","GO:0009505","GO:0031969","GO:0042170","GO:0042644","GO:0042646","GO:0005811","GO:0009295","GO:0009707","GO:0009527","no_annot_CC","GO:0016491","GO:0016829","GO:0003729","GO:0016830","GO:0016843","GO:0016987","GO:0035251","GO:0003913","GO:0009973","GO:0033741","GO:0046527","GO:0071949","GO:0046914","GO:0005506","GO:0004497","GO:0070330","GO:0080043","GO:0008171","GO:0080044","GO:0016840","GO:0004096","GO:0050373","GO:0016831","GO:0043295","GO:1900750","GO:0016712","GO:1901681","GO:0046906","GO:0004375","GO:0016594","GO:0016642","GO:0052638","GO:0097243","GO:2001147","GO:2001227"],[4,6,11,11,11,11,11,12,15,16,17,18,25,27,30,37,39,44,376,880,1498,8016,63,68,70,432,14627,77,80,93,null,140,149,54,88,397,416,3312,1,46,63,2,2447,2447,2447,71,574,3,3,3,2545,603,610,2605,87,2614,2615,632,632,2642,2710,2719,5519,3434,2804,2830,6,6,119,383,5767,7,127,127,130,401,402,8,3008,3012,811,422,140,140,813,428,3149,3186,12,12,3239,3263,14,14,15,1460,1460,4039,16,16,16,2061,1485,17,6430,5628,1511,3405,19,3482,21,591,2212,229,230,24,25,11348,26,26,1690,27,27,27,28,28,29,29,3733,30,5313,33,35,37,37,37,3899,3919,11761,1263,1266,40,40,40,745,1915,767,769,328,329,329,4065,null,9,237,293,301,51,55,56,62,66,68,null,1672,469,1028,113,20,6,142,3,3,3,180,80,1281,354,325,4,64,44,55,35,5,5,71,19,19,6,40,432,2,2,2,2,2,2,2],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,12,1,1,1,2,18,1,1,1,5,4,4,3,3,4,4,10,1,2,2,1,8,8,8,2,4,1,1,1,8,4,4,8,2,8,8,4,4,8,8,8,12,9,8,8,1,1,2,3,12,1,2,2,2,3,3,1,8,8,4,3,2,2,4,3,8,8,1,1,8,8,1,1,1,5,5,9,1,1,1,6,5,1,12,11,5,8,1,8,1,3,6,2,2,1,1,17,1,1,5,1,1,1,1,1,1,1,8,1,10,1,1,1,1,1,8,8,17,4,4,1,1,1,3,5,3,3,2,2,2,8,2,1,2,2,2,1,1,1,1,1,1,5,115,43,76,17,7,4,18,3,3,3,20,12,84,31,29,3,10,8,9,7,3,3,10,5,5,3,7,33,2,2,2,2,2,2,2],[0.00396207240277524,0.00593749775518193,0.0108597499578255,0.0108597499578255,0.0108597499578255,0.0108597499578255,0.0108597499578255,0.011841410886506,0.0147808296665232,0.0157587840329912,0.0167358138384895,0.0177119199134345,0.0245188675894235,0.0264554433826766,0.0293534431388996,0.0360835258339547,0.0379982289149229,0.0427691521706767,0.0529073433456739,0.0545185204287504,0.056840052129226,0.0570632552980758,0.0606939608584034,0.0653576141384002,0.0672168937533096,0.0676011553685207,0.0726842568703283,0.0736966714483606,0.0764605755519373,0.0883469924359658,null,1.25484487597364e-05,1.60417470534208e-05,2.36776231188795e-05,0.000102350563099042,0.000698025863889206,0.000831087636599189,0.00103954361600856,0.00104004160166371,0.00104337239345001,0.00194768376635433,0.00207904654214364,0.00244035763538266,0.00244035763538266,0.00244035763538266,0.00246570013326211,0.00270973011552939,0.00311701580984691,0.00311701580984691,0.00311701580984691,0.00313156468046933,0.0032359677894733,0.00337278628031118,0.00362623680986932,0.00367451032305311,0.00370549526082861,0.00371438580088941,0.00382859056100743,0.00382859056100743,0.00396088979088355,0.00463928739860478,0.00473553190134841,0.00489653263459281,0.00554221581530132,0.00572411427706303,0.00605667219829232,0.00622471944690062,0.00622471944690062,0.00675797982076891,0.00708837654842155,0.0071506357854646,0.00725855588960855,0.00766241908184978,0.00766241908184978,0.00801495919796392,0.00803945399463063,0.00809440137601334,0.00829136158810706,0.00875720523106411,0.00882704277956502,0.00918470782652796,0.00924058875147424,0.00924203000058363,0.00924203000058363,0.00926291078104471,0.00960213637762476,0.0114891163851419,0.0123032764641655,0.0124122965870369,0.0124122965870369,0.0135449711992255,0.0141374484402185,0.0144666051373641,0.0144666051373641,0.0154922231009709,0.0155928707931346,0.0155928707931346,0.0160170086551973,0.0165168181595368,0.0165168181595368,0.0165168181595368,0.0166276534542345,0.0166815208694746,0.0175403912889455,0.0175904889886221,0.0177411639040473,0.0178673934569459,0.0180517765376037,0.0195844756593918,0.0204854779553819,0.0216244840016301,0.0226802870978921,0.0228785891697048,0.0234442884774654,0.0236350756420669,0.0246768709663364,0.0256923036878706,0.0265339486436073,0.0267067232233131,0.0267067232233131,0.0276016226747585,0.0277201305396971,0.0277201305396971,0.0277201305396971,0.0287325266031724,0.0287325266031724,0.0297439123790111,0.0297439123790111,0.030110369468071,0.0307542888316037,0.0325161732740539,0.0337793718806442,0.035791064938108,0.0377987437733969,0.0377987437733969,0.0377987437733969,0.0380432483076367,0.0390897338292566,0.0392726582011206,0.0395493526036962,0.0398436805323981,0.040802752117775,0.040802752117775,0.040802752117775,0.040982971677698,0.0440528511268703,0.0440680930821638,0.0443542965956777,0.0452612501494667,0.0455096065774526,0.0455096065774526,0.0473492100647661,null,0.00907818411455563,0.0236739455602806,0.0350080505847825,0.0367713298992154,0.0504151732037068,0.054265208102662,0.0552253913665304,0.0609670227854618,0.0647762983300598,0.0666754128376102,null,1.81515788700584e-06,6.750632216831e-06,1.03087384433365e-05,1.05924957878821e-05,1.59423079482841e-05,5.48248118467457e-05,6.28173292241127e-05,8.83836670914496e-05,8.83836670914496e-05,8.83836670914496e-05,0.000162070431673474,0.000211772857725349,0.000242676015481734,0.000284786439439802,0.000323837226419927,0.000341747674214289,0.000502732091758004,0.000643718606381178,0.000672503027486255,0.000770557719623788,0.000825948591497468,0.000825948591497468,0.00115768974409642,0.00120190615916265,0.00120190615916265,0.00159707237114819,0.0017604887930816,0.00194893355872084,0.00198619304528545,0.00198619304528545,0.00198619304528545,0.00198619304528545,0.00198619304528545,0.00198619304528545,0.00198619304528545],[0.0396207240277524,0.0593749775518193,0.108597499578255,0.108597499578255,0.108597499578255,0.108597499578255,0.108597499578255,0.11841410886506,0.147808296665232,0.157587840329912,0.167358138384895,0.177119199134345,0.245188675894235,0.264554433826766,0.293534431388996,0.360835258339547,0.379982289149229,0.427691521706767,0.529073433456739,0.545185204287504,0.56840052129226,0.570632552980759,0.606939608584034,0.653576141384002,0.672168937533096,0.676011553685207,0.726842568703283,0.736966714483606,0.764605755519373,0.883469924359658,null,0.000238420526434991,0.000304793194014996,0.00044987483925871,0.00194466069888179,0.0132624914138949,0.0157906650953846,0.0197513287041626,0.0197607904316105,0.0198240754755501,0.0370059915607322,0.0395018843007291,0.0463667950722706,0.0463667950722706,0.0463667950722706,0.0468483025319801,0.0514848721950584,0.0592233003870913,0.0592233003870913,0.0592233003870913,0.0594997289289172,0.0614833879999926,0.0640829393259124,0.0688984993875171,0.0698156961380091,0.0704044099557437,0.0705733302168989,0.0727432206591412,0.0727432206591412,0.0752569060267874,0.0881464605734908,0.0899751061256199,0.0930341200572635,0.105302100490725,0.108758171264198,0.115076771767554,0.118269669491112,0.118269669491112,0.128401616594609,0.134679154420009,0.135862079923827,0.137912561902563,0.145585962555146,0.145585962555146,0.152284224761315,0.152749625897982,0.153793626144254,0.157535870174034,0.166386899390218,0.167713812811735,0.174509448704031,0.175571186278011,0.175598570011089,0.175598570011089,0.17599530483985,0.18244059117487,0.218293211317696,0.233762252819144,0.235833635153701,0.235833635153701,0.257354452785285,0.268611520364152,0.274865497609918,0.274865497609918,0.294352238918447,0.296264545069557,0.296264545069557,0.304323164448749,0.313819545031198,0.313819545031198,0.313819545031198,0.315925415630455,0.316948896520017,0.333267434489964,0.33421929078382,0.337082114176898,0.339480475681972,0.342983754214471,0.372105037528444,0.389224081152257,0.410865196030972,0.430925454859951,0.434693194224391,0.445441481071842,0.449066437199272,0.468860548360392,0.488153770069541,0.504145024228539,0.50742774124295,0.50742774124295,0.524430830820411,0.526682480254245,0.526682480254245,0.526682480254245,0.545918005460275,0.545918005460275,0.565134335201211,0.565134335201211,0.572097019893349,0.58433148780047,0.617807292207025,0.641808065732239,0.680030233824052,0.718176131694541,0.718176131694541,0.718176131694541,0.722821717845098,0.742704942755876,0.746180505821291,0.751437699470228,0.757029930115564,0.775252290237724,0.775252290237724,0.775252290237724,0.778676461876262,0.837004171410537,0.837293768561113,0.842731635317876,0.859963752839867,0.8646825249716,0.8646825249716,0.899634991230556,null,0.118016393489223,0.307761292283648,0.455104657602173,0.478027288689801,0.655397251648189,0.705447705334606,0.717930087764895,0.792571296211004,0.842091878290778,0.866780366888933,null,0.000490092629491576,0.00182267069854437,0.00278335937970086,0.00285997386272817,0.0043044231460367,0.0148026991986213,0.0169606788905104,0.0238635901146914,0.0238635901146914,0.0238635901146914,0.0437590165518379,0.0571786715858443,0.0655225241800682,0.0768923386487466,0.0874360511333803,0.0922718720378579,0.135737664774661,0.173804023722918,0.181575817421289,0.208050584298423,0.223006119704316,0.223006119704316,0.312576230906033,0.324514662973917,0.324514662973917,0.431209540210012,0.475331974132031,0.526212060854626,0.536272122227072,0.536272122227072,0.536272122227072,0.536272122227072,0.536272122227072,0.536272122227072,0.536272122227072],["sphingosine N-acyltransferase activity","glucose-6-phosphate dehydrogenase activity","oxidoreductase activity, acting on NAD(P)H, oxygen as acceptor","tryptamine:oxygen oxidoreductase (deaminating) activity","aminoacetone:oxygen oxidoreductase(deaminating) activity","aliphatic-amine oxidase activity","phenethylamine:oxygen oxidoreductase (deaminating) activity","pectin acetylesterase activity","primary amine oxidase activity","cyclosporin A binding","diacylglycerol O-acyltransferase activity","drug binding","acylglycerol O-acyltransferase activity","photoreceptor activity","oxidoreductase activity, acting on the CH-NH2 group of donors, oxygen as acceptor","quinone binding","oxidoreductase activity, acting on the CH-NH2 group of donors","O-methyltransferase activity","transferase activity, transferring acyl groups other than amino-acyl groups","protein serine/threonine kinase activity","kinase activity","protein binding","peptidyl-prolyl cis-trans isomerase activity","cis-trans isomerase activity","polysaccharide binding","transferase activity, transferring acyl groups","binding","N-acyltransferase activity","NADP binding","peptide binding",null,"circadian rhythm","rhythmic process","regulation of circadian rhythm","response to blue light","regulation of post-embryonic development","response to cold","aromatic compound biosynthetic process","negative regulation of protein binding","response to absence of light","response to red light","negative regulation of binding","regulation of transcription, DNA-templated","regulation of nucleic acid-templated transcription","regulation of RNA biosynthetic process","regulation of photoperiodism, flowering","regulation of multicellular organismal development","etioplast organization","regulation of protein binding","negative regulation of cellulose biosynthetic process","regulation of RNA metabolic process","response to abscisic acid","response to alcohol","regulation of nucleobase-containing compound metabolic process","positive regulation of post-embryonic development","transcription, DNA-templated","nucleic acid-templated transcription","response to temperature stimulus","regulation of multicellular organismal process","RNA biosynthetic process","regulation of cellular macromolecule biosynthetic process","regulation of macromolecule biosynthetic process","cellular biosynthetic process","organic cyclic compound biosynthetic process","regulation of cellular biosynthetic process","regulation of biosynthetic process","lignan metabolic process","lignan biosynthetic process","positive regulation of multicellular organismal process","regulation of reproductive process","biosynthetic process","long-chain fatty acid biosynthetic process","photoperiodism, flowering","response to karrikin","positive regulation of developmental process","negative regulation of cellular macromolecule biosynthetic process","negative regulation of macromolecule biosynthetic process","negative regulation of cellular carbohydrate metabolic process","regulation of gene expression","nucleobase-containing compound biosynthetic process","regulation of developmental process","negative regulation of cellular biosynthetic process","response to light intensity","photoperiodism","response to lipid","negative regulation of biosynthetic process","response to chemical","regulation of nitrogen compound metabolic process","plastid transcription","negative regulation of carbohydrate metabolic process","heterocycle biosynthetic process","regulation of primary metabolic process","pentose-phosphate shunt, oxidative branch","regulation of binding","regulation of cellulose biosynthetic process","response to external biotic stimulus","response to other organism","cellular nitrogen compound biosynthetic process","regulation of beta-glucan metabolic process","regulation of beta-glucan biosynthetic process","regulation of embryonic development","response to organic substance","response to biotic stimulus","ceramide biosynthetic process","response to stimulus","organic substance biosynthetic process","interspecies interaction between organisms","regulation of cellular metabolic process","response to molecule of fungal origin","regulation of macromolecule metabolic process","positive regulation of seed germination","negative regulation of cellular metabolic process","response to abiotic stimulus","vegetative to reproductive phase transition of meristem","response to red or far red light","long-chain fatty acid metabolic process","regulation of glucan biosynthetic process","cellular metabolic process","pentose-phosphate shunt","triglyceride biosynthetic process","response to oxygen-containing compound","glycerol metabolic process","ceramide metabolic process","alditol metabolic process","NADPH regeneration","regulation of polysaccharide biosynthetic process","neutral lipid biosynthetic process","acylglycerol biosynthetic process","regulation of metabolic process","induced systemic resistance","cellular aromatic compound metabolic process","triglyceride metabolic process","glucose 6-phosphate metabolic process","neutral lipid metabolic process","acylglycerol metabolic process","NADP metabolic process","response to stress","RNA metabolic process","organic substance metabolic process","reproductive structure development","reproductive system development","positive regulation of flower development","sphingolipid biosynthetic process","regulation of polysaccharide metabolic process","response to light stimulus","response to external stimulus","negative regulation of macromolecule metabolic process","response to radiation","negative regulation of transcription, DNA-templated","negative regulation of RNA biosynthetic process","negative regulation of nucleic acid-templated transcription","cellular macromolecule biosynthetic process",null,"plastid-encoded plastid RNA polymerase complex","plant-type cell wall","chloroplast membrane","plastid membrane","chloroplast nucleoid","plastid nucleoid","lipid droplet","nucleoid","chloroplast outer membrane","plastid outer membrane",null,"oxidoreductase activity","lyase activity","mRNA binding","carbon-carbon lyase activity","amine-lyase activity","sigma factor activity","UDP-glucosyltransferase activity","DNA photolyase activity","adenylyl-sulfate reductase activity","adenylyl-sulfate reductase (glutathione) activity","glucosyltransferase activity","FAD binding","transition metal ion binding","iron ion binding","monooxygenase activity","aromatase activity","quercetin 3-O-glucosyltransferase activity","O-methyltransferase activity","quercetin 7-O-glucosyltransferase activity","carbon-nitrogen lyase activity","catalase activity","UDP-arabinose 4-epimerase activity","carboxy-lyase activity","glutathione binding","oligopeptide binding","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen","sulfur compound binding","tetrapyrrole binding","glycine dehydrogenase (decarboxylating) activity","glycine binding","oxidoreductase activity, acting on the CH-NH2 group of donors, disulfide as acceptor","indole-3-butyrate beta-glucosyltransferase activity","flavonoid binding","camalexin binding","quercitrin binding"],["MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","CC","CC","CC","CC","CC","CC","CC","CC","CC","CC","CC","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>CLID<\/th>\n      <th>CLSZ<\/th>\n      <th>GOID<\/th>\n      <th>NodeSize<\/th>\n      <th>SampleMatch<\/th>\n      <th>Phyper<\/th>\n      <th>Padj<\/th>\n      <th>Term<\/th>\n      <th>Ont<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"autoWidth":true,"columnDefs":[{"className":"dt-right","targets":[2,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

### Plot batch GO term results

The `data.frame` generated by `GOCluster` can be plotted with the `goBarplot` function. Because of the
variable size of the sample sets, it may not always be desirable to show
the results from different DEG sets in the same bar plot. Plotting
single sample sets is achieved by subsetting the input data frame as
shown in the first line of the following example.

``` r
gos <- BatchResultslim
png("results/GOslimbarplotMF.png")
goBarplot(gos, gocat = "MF")
dev.off()
png("results/GOslimbarplotBP.png")
goBarplot(gos, gocat = "BP")
dev.off()
png("results/GOslimbarplotCC.png")
goBarplot(gos, gocat = "CC")
dev.off()
```

![](../results/GOslimbarplotMF.png)

<div align="center">

Figure 5: GO Slim Barplot for MF Ontology

</div>

</br>

## Clustering and heat maps

The following example performs hierarchical clustering on the `rlog`
transformed expression matrix subsetted by the DEGs identified in the above
differential expression analysis. It uses a Pearson correlation-based distance
measure and complete linkage for cluster joining.

``` r
library(pheatmap)
geneids <- unique(as.character(unlist(DEG_list[[1]])))
y <- assay(rlog(dds))[geneids, ]
png("results/heatmap1.png")
pheatmap(y, scale = "row", clustering_distance_rows = "correlation", 
    clustering_distance_cols = "correlation")
dev.off()
```

![](../results/heatmap1.png)

<div align="center">

Figure 6: Heat Map with Hierarchical Clustering Dendrograms of DEGs

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
    ##  [1] batchtools_0.9.14           ape_5.4-1                  
    ##  [3] ggplot2_3.3.2               systemPipeR_1.24.5         
    ##  [5] ShortRead_1.48.0            GenomicAlignments_1.26.0   
    ##  [7] SummarizedExperiment_1.20.0 Biobase_2.50.0             
    ##  [9] MatrixGenerics_1.2.0        matrixStats_0.57.0         
    ## [11] BiocParallel_1.24.1         Rsamtools_2.6.0            
    ## [13] Biostrings_2.58.0           XVector_0.30.0             
    ## [15] GenomicRanges_1.42.0        GenomeInfoDb_1.26.1        
    ## [17] IRanges_2.24.0              S4Vectors_0.28.0           
    ## [19] BiocGenerics_0.36.0         BiocStyle_2.18.0           
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
    ##  [35] rsvg_2.1                 rappdirs_0.3.1          
    ##  [37] V8_3.4.0                 Rcpp_1.0.5              
    ##  [39] jquerylib_0.1.3          vctrs_0.3.5             
    ##  [41] nlme_3.1-149             blogdown_1.2            
    ##  [43] rtracklayer_1.50.0       xfun_0.22               
    ##  [45] stringr_1.4.0            lifecycle_0.2.0         
    ##  [47] XML_3.99-0.5             edgeR_3.32.0            
    ##  [49] zlibbioc_1.36.0          scales_1.1.1            
    ##  [51] BSgenome_1.58.0          VariantAnnotation_1.36.0
    ##  [53] hms_0.5.3                RBGL_1.66.0             
    ##  [55] RColorBrewer_1.1-2       yaml_2.2.1              
    ##  [57] curl_4.3                 memoise_1.1.0           
    ##  [59] sass_0.3.1               biomaRt_2.46.0          
    ##  [61] latticeExtra_0.6-29      stringi_1.5.3           
    ##  [63] RSQLite_2.2.1            genefilter_1.72.0       
    ##  [65] checkmate_2.0.0          GenomicFeatures_1.42.1  
    ##  [67] DOT_0.1                  rlang_0.4.8             
    ##  [69] pkgconfig_2.0.3          bitops_1.0-6            
    ##  [71] evaluate_0.14            lattice_0.20-41         
    ##  [73] purrr_0.3.4              bit_4.0.4               
    ##  [75] tidyselect_1.1.0         GSEABase_1.52.0         
    ##  [77] AnnotationForge_1.32.0   magrittr_2.0.1          
    ##  [79] bookdown_0.21            R6_2.5.0                
    ##  [81] generics_0.1.0           base64url_1.4           
    ##  [83] DelayedArray_0.16.0      DBI_1.1.0               
    ##  [85] withr_2.3.0              pillar_1.4.7            
    ##  [87] survival_3.2-10          RCurl_1.98-1.2          
    ##  [89] tibble_3.0.4             crayon_1.3.4            
    ##  [91] BiocFileCache_1.14.0     rmarkdown_2.7           
    ##  [93] jpeg_0.1-8.1             progress_1.2.2          
    ##  [95] locfit_1.5-9.4           grid_4.0.5              
    ##  [97] data.table_1.13.2        blob_1.2.1              
    ##  [99] Rgraphviz_2.34.0         digest_0.6.27           
    ## [101] xtable_1.8-4             brew_1.0-6              
    ## [103] openssl_1.4.3            munsell_0.5.0           
    ## [105] bslib_0.2.4              askpass_1.1

## Funding

This project was supported by funds from the National Institutes of Health (NIH).

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-H_Backman2016-bt" class="csl-entry">

H Backman, Tyler W, and Thomas Girke. 2016. “<span class="nocase">systemPipeR: NGS workflow and report generation environment</span>.” *BMC Bioinformatics* 17 (1): 388. <https://doi.org/10.1186/s12859-016-1241-0>.

</div>

<div id="ref-Howard2013-fq" class="csl-entry">

Howard, Brian E, Qiwen Hu, Ahmet Can Babaoglu, Manan Chandra, Monica Borghi, Xiaoping Tan, Luyan He, et al. 2013. “High-Throughput RNA Sequencing of Pseudomonas-Infected Arabidopsis Reveals Hidden Transcriptome Complexity and Novel Splice Variants.” *PLoS One* 8 (10): e74183. <https://doi.org/10.1371/journal.pone.0074183>.

</div>

<div id="ref-Kim2015-ve" class="csl-entry">

Kim, Daehwan, Ben Langmead, and Steven L Salzberg. 2015. “HISAT: A Fast Spliced Aligner with Low Memory Requirements.” *Nat. Methods* 12 (4): 357–60.

</div>

<div id="ref-Lawrence2013-kt" class="csl-entry">

Lawrence, Michael, Wolfgang Huber, Hervé Pagès, Patrick Aboyoun, Marc Carlson, Robert Gentleman, Martin T Morgan, and Vincent J Carey. 2013. “Software for Computing and Annotating Genomic Ranges.” *PLoS Comput. Biol.* 9 (8): e1003118. <https://doi.org/10.1371/journal.pcbi.1003118>.

</div>

<div id="ref-Robinson2010-uk" class="csl-entry">

Robinson, M D, D J McCarthy, and G K Smyth. 2010. “edgeR: A Bioconductor Package for Differential Expression Analysis of Digital Gene Expression Data.” *Bioinformatics* 26 (1): 139–40. <https://doi.org/10.1093/bioinformatics/btp616>.

</div>

</div>
