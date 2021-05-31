---
title: "RNA-Seq Workflow Template" 
author: "Author: First Last"
date: "Last update: 31 May, 2021" 
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
Source code downloads: &nbsp; &nbsp;
[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/spWFtemplates/sprnaseq.Rmd) ] &nbsp; &nbsp; 
[ [.html](https://girke.bioinformatics.ucr.edu/GEN242/custom/spWFtemplates/sprnaseq.html) ] &nbsp; &nbsp; 
[ [old version .Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/sprnaseq/back_sprnaseq_rmd) ]
</div>


## Introduction

This report describes the analysis of the RNA-Seq data set from 
Howard et al [-@Howard2013-fq]. The corresponding FASTQ files were downloaded from
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

1. If a remote system or cluster is used, then users need to log in to the
   remote system first. The following applies to an HPC cluster (_e.g._ HPCC
   cluster). 

   A terminal application needs to be used to log in to a user's cluster account. Next, one
   can open an interactive session on a computer node with `srun`. More details about
   argument settings for `srun` are available in this [HPCC
   manual](http://hpcc.ucr.edu/manuals_linux-cluster_jobs.html#partitions) or
   the HPCC section of this website
   [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/linux/linux/#job-submission-with-sbatch).
   Next, load the R version required for running the workflow with `module load`. Sometimes it may be necessary to
   first unload an active software version before loading another version, _e.g._ `module unload R`.

```sh
srun --x11 --partition=gen242 --mem=20gb --cpus-per-task 8 --ntasks 1 --time 2:00:00 --pty bash -l
module unload R; module load R/4.0.3_gcc-8.3.0
```

2. Load a workflow template with the `genWorkenvir` function. This can be done from the command-line or from within R. 
   However, only one of the two options needs to be used.

From command-line

```sh
$ Rscript -e "systemPipeRdata::genWorkenvir(workflow='rnaseq')"
$ cd rnaseq
```

From R


```r
library(systemPipeRdata)
genWorkenvir(workflow = "rnaseq")
setwd("rnaseq")
```

3. Optional: if the user wishes to use another `Rmd` file than the template instance provided by the `genWorkenvir` function, then it can be copied or downloaded 
   into the root directory of the workflow environment (_e.g._ with `cp` or `wget`).

4. Now one can open from the root directory of the workflow the corresponding R Markdown script (_e.g._ systemPipeChIPseq.Rmd) using an R IDE, such as _nvim-r_, _ESS_ or RStudio. 
   Subsequently, the workflow can be run as outlined below. For learning purposes it is recommended to run workflows for the first time interactively. Once all workflow steps are 
   understood and possibly modified to custom needs, one can run the workflow from start to finish with a single command using `rmarkdown::render()` or `runWF()`.


### Load packages

The `systemPipeR` package needs to be loaded to perform the analysis 
steps shown in this report [@H_Backman2016-bt]. The package allows users
to run the entire analysis workflow interactively or with a single command 
while also generating the corresponding analysis report. For details
see `systemPipeR's` main [vignette](http://www.bioconductor.org/packages/devel/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html).


```r
library(systemPipeR)
```

To apply workflows to custom data, the user needs to modify the _`targets`_ file and if
necessary update the corresponding parameter (_`.cwl`_ and _`.yml`_) files. 
A collection of pre-generated _`.cwl`_ and _`.yml`_ files are provided in the _`param/cwl`_ subdirectory 
of each workflow template. They are also viewable in the GitHub repository of _`systemPipeRdata`_ ([see
here](https://github.com/tgirke/systemPipeRdata/tree/master/inst/extdata/param/cwl)).
For more information of the structure of the *targets* file, please consult the documentation 
[here](http://www.bioconductor.org/packages/release/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html#25_structure_of_targets_file). More details about the new parameter files from systemPipeR can be found [here](http://www.bioconductor.org/packages/release/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html#26_structure_of_the_new_param_files_and_construct_sysargs2_container). 

### Import custom functions

Custem functions for the challenge projects can be imported with the source
command from a local R script (here [challengeProject_Fct.R](https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/tutorials/spchipseq/challengeProject_Fct.R)).  Skip this step if such a
script is not available.  Alternatively, these functions can be loaded from a
custom R package. 


```r
source("challengeProject_Fct.R")
```

### Experiment definition provided by `targets` file

The `targets` file defines all FASTQ files and sample comparisons of the analysis workflow.
If needed the tab separated (TSV) version of this file can be downloaded from [here](https://github.com/tgirke/GEN242/tree/main/content/en/assignments/Projects/targets_files)
and the corresponding Google Sheet is [here](https://docs.google.com/spreadsheets/d/1DTgTGlZZscSPjlHOGdJC8QK4vvimN1BORjXKXzd_cfA/edit#gid=472150521).


```r
targetspath <- "targetsPE.txt"
targets <- read.delim(targetspath, comment.char = "#")
DT::datatable(targets, options = list(scrollX = TRUE, autoWidth = TRUE))
```

```{=html}
<div id="htmlwidget-34065808aa3b695107be" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-34065808aa3b695107be">{"x":{"filter":"none","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"],["./data/SRR446027_1.fastq.gz","./data/SRR446028_1.fastq.gz","./data/SRR446029_1.fastq.gz","./data/SRR446030_1.fastq.gz","./data/SRR446031_1.fastq.gz","./data/SRR446032_1.fastq.gz","./data/SRR446033_1.fastq.gz","./data/SRR446034_1.fastq.gz","./data/SRR446035_1.fastq.gz","./data/SRR446036_1.fastq.gz","./data/SRR446037_1.fastq.gz","./data/SRR446038_1.fastq.gz","./data/SRR446039_1.fastq.gz","./data/SRR446040_1.fastq.gz","./data/SRR446041_1.fastq.gz","./data/SRR446042_1.fastq.gz","./data/SRR446043_1.fastq.gz","./data/SRR446044_1.fastq.gz"],["./data/SRR446027_2.fastq.gz","./data/SRR446028_2.fastq.gz","./data/SRR446029_2.fastq.gz","./data/SRR446030_2.fastq.gz","./data/SRR446031_2.fastq.gz","./data/SRR446032_2.fastq.gz","./data/SRR446033_2.fastq.gz","./data/SRR446034_2.fastq.gz","./data/SRR446035_2.fastq.gz","./data/SRR446036_2.fastq.gz","./data/SRR446037_2.fastq.gz","./data/SRR446038_2.fastq.gz","./data/SRR446039_2.fastq.gz","./data/SRR446040_2.fastq.gz","./data/SRR446041_2.fastq.gz","./data/SRR446042_2.fastq.gz","./data/SRR446043_2.fastq.gz","./data/SRR446044_2.fastq.gz"],["M1A","M1B","A1A","A1B","V1A","V1B","M6A","M6B","A6A","A6B","V6A","V6B","M12A","M12B","A12A","A12B","V12A","V12B"],["M1","M1","A1","A1","V1","V1","M6","M6","A6","A6","V6","V6","M12","M12","A12","A12","V12","V12"],["Mock.1h.A","Mock.1h.B","Avr.1h.A","Avr.1h.B","Vir.1h.A","Vir.1h.B","Mock.6h.A","Mock.6h.B","Avr.6h.A","Avr.6h.B","Vir.6h.A","Vir.6h.B","Mock.12h.A","Mock.12h.B","Avr.12h.A","Avr.12h.B","Vir.12h.A","Vir.12h.B"],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],["23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012","23-Mar-2012"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>FileName1<\/th>\n      <th>FileName2<\/th>\n      <th>SampleName<\/th>\n      <th>Factor<\/th>\n      <th>SampleLong<\/th>\n      <th>Experiment<\/th>\n      <th>Date<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"autoWidth":true,"columnDefs":[{"className":"dt-right","targets":6},{"orderable":false,"targets":0}],"order":[],"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

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
`SYSargs2` instance, _e.g._ running the NGS alignments using the
trimmed FASTQ files.

Construct _`SYSargs2`_ object from _`cwl`_ and _`yml`_ param and _`targets`_ files.


```r
dir_path <- "param/cwl/preprocessReads/trim-pe"
trim <- loadWorkflow(targets = targetspath, wf_file = "trim-pe.cwl", 
    input_file = "trim-pe.yml", dir_path = dir_path)
trim <- renderWF(trim, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
    FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
trim
output(trim)[1:2]
```


```r
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


```r
fqlist <- seeFastq(fastq = infile1(trim), batchsize = 10000, 
    klength = 8)
png("./results/fastqReport.png", height = 18 * 96, width = 4 * 
    96 * length(fqlist))
seeFastqPlot(fqlist)
dev.off()
```

![](../results/fastqReport.png)
<div align="center">Figure 1: FASTQ quality report for 18 samples</div></br>

## Alignments

### Read mapping with `HISAT2`

The following steps will demonstrate how to use the short read aligner `Hisat2`
[@Kim2015-ve] in both interactive job submissions and batch submissions to
queuing systems of clusters using the _`systemPipeR's`_ new CWL command-line interface.

Build `Hisat2` index.


```r
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











































