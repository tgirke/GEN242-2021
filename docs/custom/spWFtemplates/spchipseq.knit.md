---
title: "ChIP-Seq Workflow Template" 
author: "Author: First Last"
date: "Last update: 30 May, 2021" 
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
Rscript -e "rmarkdown::render('spchipseq.Rmd', c('BiocStyle::html_document'), clean=F); knitr::knit('spchipseq.Rmd', tangle=TRUE)"
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




# Introduction

The following analyzes the ChIP-Seq data from Kaufman et al. [-@Kaufmann2010-me] using 
for peak calling MACS2 where the uninduced sample serves as input (reference). 
The details about all download steps are provided [here](https://girke.bioinformatics.ucr.edu/GEN242/assignments/projects/project_data/).

Users want to extend this section to provide all background information relevant for this 
ChIP-Seq project.


## Experimental design

Typically, users want to specify here all information relevant for the
analysis of their NGS study. This includes detailed descriptions of
FASTQ files, experimental design, reference genome, gene annotations,
etc.

## Workflow environment

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
$ Rscript -e "systemPipeRdata::genWorkenvir(workflow='chipseq')"
$ cd chipseq
```

From R


```r
library(systemPipeRdata)
genWorkenvir(workflow = "chipseq")
setwd("chipseq")
```

3. Optional: if the user wishes to use another `Rmd` file than the template instance provided by the `genWorkenvir` function, then it can be copied or downloaded 
   into the root directory of the workflow environment (_e.g._ with `cp` or `wget`).

4. Now one can open from the root directory of the workflow the corresponding R Markdown script (_e.g._ systemPipeChIPseq.Rmd) using an R IDE, such as _nvim-r_, _ESS_ or RStudio. 
   Subsequently, the workflow can be run as outlined below. For learning purposes it is recommended to run workflows for the first time interactively. Once all workflow steps are 
   understood and possibly modified to custom needs, one can run the workflow from start to finish with a single command using `rmarkdown::render()` or `runWF()`.


## Load packages

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

## Import custom functions

Custem functions for the challenge projects can be imported with the source
command from a local R script (here [challengeProject_Fct.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/spchipseq/challengeProject_Fct.R)).  Skip this step if such a
script is not available.  Alternatively, these functions can be loaded from a
custom R package. 


```r
source("challengeProject_Fct.R")
```

## Experiment definition provided by `targets` file

The `targets` file defines all FASTQ files and sample comparisons of the analysis workflow. 
If needed the tab separated (TSV) version of this file can be downloaded from [here](https://github.com/tgirke/GEN242/tree/main/content/en/assignments/Projects/targets_files)
and the corresponding Google Sheet is [here](https://docs.google.com/spreadsheets/d/1w9V3JDOsXR8qW_qNqoXO0E0UR_Oev4k7-o2y6NMDTt8/edit#gid=472150521).













































