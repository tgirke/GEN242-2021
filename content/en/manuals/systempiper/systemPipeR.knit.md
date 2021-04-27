---
title: "systemPipeR: Workflow design and reporting generation environment" 
author: "Author: Daniela Cassol, Le Zhang and Thomas Girke"
date: "Last update: 26 April, 2021" 
output:
  BiocStyle::html_document:
    toc_float: true
    code_folding: show
  BiocStyle::pdf_document: default
package: systemPipeR
vignette: |
  %\VignetteEncoding{UTF-8}
  %\VignetteIndexEntry{systemPipeR: Workflow design and reporting generation environment}
  %\VignetteEngine{knitr::rmarkdown}
fontsize: 14pt
bibliography: bibtex.bib
editor_options: 
  chunk_output_type: console
weight: 6
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
Rscript -e "rmarkdown::render('systemPipeR.Rmd', c('BiocStyle::html_document'), clean=F); knitr::knit('systemPipeR.Rmd', tangle=TRUE)"; Rscript ../md2jekyll.R systemPipeR.knit.md 2; Rscript -e "rmarkdown::render('systemPipeR.Rmd', c('BiocStyle::pdf_document'))"
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

[_`systemPipeR`_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html) provides flexible utilities for building and running automated end-to-end analysis workflows for a wide range of research applications, including next-generation sequencing (NGS) experiments, such as RNA-Seq, ChIP-Seq, VAR-Seq and Ribo-Seq [@H_Backman2016-bt]. Important features include a uniform workflow interface across different data analysis applications, automated report generation, and support for running both R and command-line software, such as NGS aligners or peak/variant callers, on local computers or compute clusters (Figure 1). The latter supports interactive job submissions and batch submissions to queuing systems of clusters. For instance, _`systemPipeR`_ can be used with most command-line aligners such as `BWA` [@Li2013-oy; @Li2009-oc], `HISAT2` [@Kim2015-ve], `TopHat2` [@Kim2013-vg] and `Bowtie2` [@Langmead2012-bs], as well as the R-based NGS aligners [_`Rsubread`_](http://www.bioconductor.org/packages/devel/bioc/html/Rsubread.html) [@Liao2013-bn] and [_`gsnap (gmapR)`_](http://www.bioconductor.org/packages/devel/bioc/html/gmapR.html) [@Wu2010-iq]. Efficient handling of complex sample sets (_e.g._ FASTQ/BAM files) and experimental designs are facilitated by a well-defined sample annotation infrastructure which improves reproducibility and user-friendliness of many typical analysis workflows in the NGS area [@Lawrence2013-kt]. 

The main motivation and advantages of using _`systemPipeR`_ for complex data analysis tasks are:

1. Facilitates the design of complex NGS workflows involving multiple R/Bioconductor packages
2. Common workflow interface for different NGS applications
3. Makes NGS analysis with Bioconductor utilities more accessible to new users
4. Simplifies usage of command-line software from within R
5. Reduces the complexity of using compute clusters for R and command-line software
6. Accelerates runtime of workflows via parallelization on computer systems with multiple CPU cores and/or multiple compute nodes
6. Improves reproducibility by automating analyses and generation of analysis reports 

<center><img src="../utilities.png"></center>

**Figure 1:** Relevant features in _`systemPipeR`_.
Workflow design concepts are illustrated under (A & B). Examples of
*systemPipeR’s* visualization functionalities are given under (C). </br>

A central concept for designing workflows within the _`systemPipeR`_ environment 
is the use of workflow management containers. In previous versions, _`systemPipeR`_ 
used a custom command-line interface called _`SYSargs`_ (see Figure 3) and for 
this purpose will continue to be supported for some time. With the latest [Bioconductor Release 3.9](http://www.bioconductor.org/packages/release/bioc/html/systemPipeR.html), 
we are adopting for this functionality the widely used community standard 
[Common Workflow Language](https://www.commonwl.org/) (CWL) for describing 
analysis workflows in a generic and reproducible manner, introducing _`SYSargs2`_
workflow control class (see Figure 2). Using this community standard in _`systemPipeR`_
has many advantages. For instance, the integration of CWL allows running _`systemPipeR`_
workflows from a single specification instance either entirely from within R, from various command-line
wrappers (e.g., *cwl-runner*) or from other languages (*, e.g.,* Bash or Python).
_`systemPipeR`_ includes support for both command-line and R/Bioconductor software 
as well as resources for containerization, parallel evaluations on computer clusters 
along with the automated generation of interactive analysis reports.

An important feature of _`systemPipeR's`_ CWL interface is that it provides two
options to run command-line tools and workflows based on CWL. First, one can
run CWL in its native way via an R-based wrapper utility for *cwl-runner* or
*cwl-tools* (CWL-based approach). Second, one can run workflows using CWL's
command-line and workflow instructions from within R (R-based approach). In the
latter case the same CWL workflow definition files (*e.g.* `*.cwl` and `*.yml`)
are used but rendered and executed entirely with R functions defined by
_`systemPipeR`_, and thus use CWL mainly as a command-line and workflow
definition format rather than software to run workflows. In this regard
_`systemPipeR`_ also provides several convenience functions that are useful for
designing and debugging workflows, such as a command-line rendering function to
retrieve the exact command-line strings for each data set and processing step
prior to running a command-line.

This overview introduces the design of a new CWL S4 class in _`systemPipeR`_, 
as well as the custom command-line interface, combined with the overview of all
the common analysis steps of NGS experiments.

## Workflow design structure using _`SYSargs2`_ 

The flexibility of _`systemPipeR's`_ new interface workflow control class is the driving factor behind 
the use of as many steps necessary for the analysis, as well as the connection 
between command-line- or R-based software. The connectivity among all
workflow steps is achieved by the _`SYSargs2`_ workflow control class (see Figure 3).
This S4 class is a list-like container where each instance stores all the
input/output paths and parameter components required for a particular data
analysis step. _`SYSargs2`_ instances are generated by two constructor
functions, *loadWorkflow* and *renderWF*, using as data input *targets* or
*yaml* files as well as two *cwl* parameter files (for details see below). When
running preconfigured workflows, the only input the user needs to provide is
the initial *targets* file containing the paths to the input files (*e.g.*
FASTQ) along with unique sample labels. Subsequent targets instances are
created automatically. The parameters required for running command-line
software is provided by the parameter (*.cwl*) files described below. 

We also introduce the *`SYSargs2Pipe`* class that organizes one or many
SYSargs2 containers in a single compound object capturing all information
required to run, control and monitor complex workflows from start to finish. This
design enhances the *`systemPipeR`* workflow framework with a generalized,
flexible, and robust design.

<center><img src="../SYS_WF.png"></center>

**Figure 2:** Workflow steps with input/output file operations are controlled by 
_`SYSargs2`_ objects. Each _`SYSargs2`_ instance is constructed from one *targets* 
and two *param* files. The only input provided by the user is the initial *targets* 
file. Subsequent *targets* instances are created automatically, from the previous 
output files. Any number of predefined or custom workflow steps are supported. One
or many _`SYSargs2`_ objects are organized in an *`SYSargs2Pipe`* container.

## Workflow Management using _`SYSargsList`_ 

**systemPipeR** allows creation (multi-step analyses) and execution of workflow entirely for R, with control, flexibility, and scalability of the all process. The execution of the workflow can be sent to a HPC, can be parallelizes, accelerating results acquisition. A workflow management system provides an infrastructure for the set-up, performance and monitoring of a defined sequence of tasks, arranged as a workflow application.

<center><img src="../sysargslist.png"></center>
**Figure 3:** Workflow Management using _`SYSargsList`_.

## Workflow design structure using _`SYSargs`_: Previous version

Instances of this S4 object class are constructed by the _`systemArgs`_ function 
from two simple tabular files: a _`targets`_ file and a _`param`_ file. The latter
is optional for workflow steps lacking command-line software. Typically, a 
_`SYSargs`_ instance stores all sample-level inputs as well as the paths to the 
corresponding outputs generated by command-line- or R-based software generating 
sample-level output files, such as read preprocessors (trimmed/filtered FASTQ 
files), aligners (SAM/BAM files), variant callers (VCF/BCF files) or peak callers 
(BED/WIG files). Each sample level input/output operation uses its own _`SYSargs`_ 
instance. The outpaths of _`SYSargs`_ usually define the sample inputs for the 
next _`SYSargs`_ instance. This connectivity is established by writing the 
outpaths with the _`writeTargetsout`_ function to a new _`targets`_ file that 
serves as input to the next _`systemArgs`_ call. Typically, the user has to 
provide only the initial _`targets`_ file. All downstream _`targets`_ files are 
generated automatically. By chaining several _`SYSargs`_ steps together one can 
construct complex workflows involving many sample-level input/output file 
operations with any combination of command-line or R-based software. 

<center><img src="../SystemPipeR_Workflow.png"></center>

**Figure 4:** Workflow design structure of _`systemPipeR`_ using _`SYSargs`_. 

# Getting Started

## Installation

The R software for running [_`systemPipeR`_](http://www.bioconductor.org/packages/devel/bioc/html/systemPipeR.html) can be downloaded from [_CRAN_](http://cran.at.r-project.org/). The _`systemPipeR`_ environment can be installed from the R console using the [_`BiocManager::install`_](https://cran.r-project.org/web/packages/BiocManager/index.html) command. The associated data package [_`systemPipeRdata`_](http://www.bioconductor.org/packages/devel/data/experiment/html/systemPipeRdata.html) can be installed the same way. The latter is a helper package for generating _`systemPipeR`_ workflow environments with a single command containing all parameter files and sample data required to quickly test and run workflows. 


```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("systemPipeR")
BiocManager::install("systemPipeRdata")
```

Please note that if you desire to use a third-party command line tool, the particular tool and dependencies need to be installed and exported in your PATH. See [details](#tools). 

## Loading package and documentation


```r
library("systemPipeR")  # Loads the package
library(help = "systemPipeR")  # Lists package info
vignette("systemPipeR")  # Opens vignette
```

## Load sample data and workflow templates

The mini sample FASTQ files used by this overview vignette as well as the 
associated workflow reporting vignettes can be loaded via the _`systemPipeRdata`_ 
package as shown below. The chosen data set [`SRP010938`](http://www.ncbi.nlm.nih.gov/sra/?term=SRP010938) 
obtains 18 paired-end (PE) read sets from _Arabidposis thaliana_ [@Howard2013-fq]. 
To minimize processing time during testing, each FASTQ file has been subsetted to 
90,000-100,000 randomly sampled PE reads that map to the first 100,000 nucleotides 
of each chromosome of the _A. thalina_ genome. The corresponding reference genome 
sequence (FASTA) and its GFF annotation files (provided in the same download) have 
been truncated accordingly. This way the entire test sample data set requires 
less than 200MB disk storage space. A PE read set has been chosen for this test 
data set for flexibility, because it can be used for testing both types of analysis 
routines requiring either SE (single-end) reads or PE reads. 

The following generates a fully populated _`systemPipeR`_ workflow environment 
(here for RNA-Seq) in the current working directory of an R session. At this time 
the package includes workflow templates for RNA-Seq, ChIP-Seq, VAR-Seq, and Ribo-Seq. 
Templates for additional NGS applications will be provided in the future.


```r
library(systemPipeRdata)
genWorkenvir(workflow = "rnaseq")
setwd("rnaseq")
```

If you desire run this tutorial with your data set, please follow the instruction here:


```r
library(systemPipeRdata)
genWorkenvir(workflow = "new", mydirname = "FEB_project")
```

### Workflow template from an individual's package

The package provides pre-configured workflows and reporting templates for a wide range of NGS applications that are listed [here](https://github.com/tgirke/systemPipeR/tree/devel#workflow). Additional workflow templates will be provided in the future. 
If you desire to use an individual package and version, follow the instruction below:


```r
library(systemPipeRdata)
genWorkenvir(workflow = NULL, package_repo = "systemPipeR/systemPipeRIBOseq", ref = "master", 
    subdir = NULL)
```


```r
library(systemPipeRdata)
genWorkenvir(workflow = NULL, package_repo = "systemPipeR/systemPipeRNAseq", ref = "singleMachine", 
    subdir = NULL)
```

## Directory Structure

The working environment of the sample data loaded in the previous step contains
the following pre-configured directory structure (Figure 4). Directory names are indicated
in  <span style="color:grey">***green***</span>. Users can change this
structure as needed, but need to adjust the code in their workflows
accordingly. 

* <span style="color:green">_**workflow/**_</span> (*e.g.* *rnaseq/*) 
    + This is the root directory of the R session running the workflow.
    + Run script ( *\*.Rmd*) and sample annotation (*targets.txt*) files are located here.
    + Note, this directory can have any name (*e.g.* <span style="color:green">_**rnaseq**_</span>, <span style="color:green">_**varseq**_</span>). Changing its name does not require any modifications in the run script(s).
  + **Important subdirectories**: 
    + <span style="color:green">_**param/**_</span> 
        + Stores non-CWL parameter files such as: *\*.param*, *\*.tmpl* and *\*.run.sh*. These files are only required for backwards compatibility to run old workflows using the previous custom command-line interface.
        + <span style="color:green">_**param/cwl/**_</span>: This subdirectory stores all the CWL parameter files. To organize workflows, each can have its own subdirectory, where all `CWL param` and `input.yml` files need to be in the same subdirectory. 
    + <span style="color:green">_**data/**_ </span>
        + FASTQ files
        + FASTA file of reference (*e.g.* reference genome)
        + Annotation files
        + etc.
    + <span style="color:green">_**results/**_</span>
        + Analysis results are usually written to this directory, including: alignment, variant and peak files (BAM, VCF, BED); tabular result files; and image/plot files
        + Note, the user has the option to organize results files for a given sample and analysis step in a separate subdirectory.

<center><img src="../SYSdir.png"></center>

**Figure 5:** *systemPipeR's* preconfigured directory structure.

The following parameter files are included in each workflow template:

1. *`targets.txt`*: initial one provided by user; downstream *`targets_*.txt`* files are generated automatically
2. *`*.param/cwl`*: defines parameter for input/output file operations, *e.g.*:
    + *`hisat2-se/hisat2-mapping-se.cwl`* 
    + *`hisat2-se/hisat2-mapping-se.yml`*
3. *`*_run.sh`*: optional bash scripts 
4. Configuration files for computer cluster environments (skip on single machines):
    + *`.batchtools.conf.R`*: defines the type of scheduler for *`batchtools`* pointing to template file of cluster, and located in user's home directory
    + *`*.tmpl`*: specifies parameters of scheduler used by a system, *e.g.* Torque, SGE, Slurm, etc.

## Structure of _`targets`_ file

The _`targets`_ file defines all input files (_e.g._ FASTQ, BAM, BCF) and sample 
comparisons of an analysis workflow. The following shows the format of a sample 
_`targets`_ file included in the package. It also can be viewed and downloaded 
from _`systemPipeR`_'s GitHub repository [here](https://github.com/tgirke/systemPipeR/blob/master/inst/extdata/targets.txt). 
In a target file with a single type of input files, here FASTQ files of single-end (SE) reads, the first three columns are mandatory including their column
names, while it is four mandatory columns for FASTQ files of PE reads. All 
subsequent columns are optional and any number of additional columns can be added as needed.

Users should note here, the usage of targets files is optional when using
*systemPipeR's* new CWL interface. They can be replaced by a standard YAML
input file used by CWL. Since for organizing experimental variables targets
files are extremely useful and user-friendly. Thus, we encourage users to keep using 
them. 

### Structure of _`targets`_ file for single-end (SE) samples


```r
library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package = "systemPipeR")
read.delim(targetspath, comment.char = "#")[1:4, ]
```

```
##                      FileName SampleName Factor SampleLong Experiment
## 1 ./data/SRR446027_1.fastq.gz        M1A     M1  Mock.1h.A          1
## 2 ./data/SRR446028_1.fastq.gz        M1B     M1  Mock.1h.B          1
## 3 ./data/SRR446029_1.fastq.gz        A1A     A1   Avr.1h.A          1
## 4 ./data/SRR446030_1.fastq.gz        A1B     A1   Avr.1h.B          1
##          Date
## 1 23-Mar-2012
## 2 23-Mar-2012
## 3 23-Mar-2012
## 4 23-Mar-2012
```

To work with custom data, users need to generate a _`targets`_ file containing 
the paths to their own FASTQ files and then provide under _`targetspath`_ the
path to the corresponding _`targets`_ file. 

### Structure of _`targets`_ file for paired-end (PE) samples

For paired-end (PE) samples, the structure of the targets file is similar, where
users need to provide two FASTQ path columns: *`FileName1`* and *`FileName2`* 
with the paths to the PE FASTQ files. 


```r
targetspath <- system.file("extdata", "targetsPE.txt", package = "systemPipeR")
read.delim(targetspath, comment.char = "#")[1:2, 1:6]
```

```
##                     FileName1                   FileName2 SampleName Factor
## 1 ./data/SRR446027_1.fastq.gz ./data/SRR446027_2.fastq.gz        M1A     M1
## 2 ./data/SRR446028_1.fastq.gz ./data/SRR446028_2.fastq.gz        M1B     M1
##   SampleLong Experiment
## 1  Mock.1h.A          1
## 2  Mock.1h.B          1
```

### Sample comparisons

Sample comparisons are defined in the header lines of the _`targets`_ file 
starting with '``# <CMP>``'. 


```r
readLines(targetspath)[1:4]
```

```
## [1] "# Project ID: Arabidopsis - Pseudomonas alternative splicing study (SRA: SRP010938; PMID: 24098335)"                                                                              
## [2] "# The following line(s) allow to specify the contrasts needed for comparative analyses, such as DEG identification. All possible comparisons can be specified with 'CMPset: ALL'."
## [3] "# <CMP> CMPset1: M1-A1, M1-V1, A1-V1, M6-A6, M6-V6, A6-V6, M12-A12, M12-V12, A12-V12"                                                                                             
## [4] "# <CMP> CMPset2: ALL"
```

The function _`readComp`_ imports the comparison information and stores it in a 
_`list`_. Alternatively, _`readComp`_ can obtain the comparison information from 
the corresponding _`SYSargs`_ object (see below). Note, these header lines are 
optional. They are mainly useful for controlling comparative analyses according 
to certain biological expectations, such as identifying differentially expressed
genes in RNA-Seq experiments based on simple pair-wise comparisons.
 

```r
readComp(file = targetspath, format = "vector", delim = "-")
```

```
## $CMPset1
## [1] "M1-A1"   "M1-V1"   "A1-V1"   "M6-A6"   "M6-V6"   "A6-V6"   "M12-A12"
## [8] "M12-V12" "A12-V12"
## 
## $CMPset2
##  [1] "M1-A1"   "M1-V1"   "M1-M6"   "M1-A6"   "M1-V6"   "M1-M12"  "M1-A12" 
##  [8] "M1-V12"  "A1-V1"   "A1-M6"   "A1-A6"   "A1-V6"   "A1-M12"  "A1-A12" 
## [15] "A1-V12"  "V1-M6"   "V1-A6"   "V1-V6"   "V1-M12"  "V1-A12"  "V1-V12" 
## [22] "M6-A6"   "M6-V6"   "M6-M12"  "M6-A12"  "M6-V12"  "A6-V6"   "A6-M12" 
## [29] "A6-A12"  "A6-V12"  "V6-M12"  "V6-A12"  "V6-V12"  "M12-A12" "M12-V12"
## [36] "A12-V12"
```

## Structure of the new _`param`_ files and construct _`SYSargs2`_ container

_`SYSargs2`_ stores all the information and instructions needed for processing
a set of input files with a single or many command-line steps within a workflow
(*i.e.* several components of the software or several independent software tools).
The _`SYSargs2`_ object is created and fully populated with the *loadWF* 
and *renderWF* functions, respectively. 

In CWL, files with the extension *`.cwl`* define the parameters of a chosen 
command-line step or workflow, while files with the extension *`.yml`* define 
the input variables of command-line steps. Note, input variables provided
by a *targets* file can be passed on to a _`SYSargs2`_ instance via the *inputvars* 
argument of the *renderWF* function. 





The following imports a *`.cwl`* file (here *`hisat2-mapping-se.cwl`*) for running
the short read aligner HISAT2 [@Kim2015-ve]. The *loadWorkflow* and *renderWF* 
functions render the proper command-line strings for each sample and software tool.


```r
library(systemPipeR)
targets <- system.file("extdata", "targets.txt", package = "systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package = "systemPipeR")
WF <- loadWF(targets = targets, wf_file = "hisat2-mapping-se.cwl", input_file = "hisat2-mapping-se.yml", 
    dir_path = dir_path)

WF <- renderWF(WF, inputvars = c(FileName = "_FASTQ_PATH1_", SampleName = "_SampleName_"))
```

Several accessor methods are available that are named after the slot names of the _`SYSargs2`_ object. 

```r
names(WF)
```

```
##  [1] "targets"       "targetsheader" "modules"       "wf"           
##  [5] "clt"           "yamlinput"     "cmdlist"       "input"        
##  [9] "output"        "cwlfiles"      "inputvars"
```

Of particular interest is the *`cmdlist()`* method. It constructs the system
commands for running command-line software as specified by a given *`.cwl`*
file combined with the paths to the input samples (*e.g.* FASTQ files) provided
by a *`targets`* file. The example below shows the *`cmdlist()`* output for
running HISAT2 on the first SE read sample. Evaluating the output of
*`cmdlist()`* can be very helpful for designing and debugging *`.cwl`* files
of new command-line software or changing the parameter settings of existing
ones.  


```r
cmdlist(WF)[1]
```

```
## $M1A
## $M1A$`hisat2-mapping-se`
## [1] "hisat2 -S ./results/M1A.sam  -x ./data/tair10.fasta  -k 1  --min-intronlen 30  --max-intronlen 3000  -U ./data/SRR446027_1.fastq.gz --threads 4"
```

The output components of _`SYSargs2`_ define the expected output files for 
each step in the workflow; some of which are the input for the next workflow step, 
here next _`SYSargs2`_ instance (see Figure 2).


```r
output(WF)[1]
```

```
## $M1A
## $M1A$`hisat2-mapping-se`
## [1] "./results/M1A.sam"
```

```r
modules(WF)
```

```
##        module1 
## "hisat2/2.1.0"
```

```r
targets(WF)[1]
```

```
## $M1A
## $M1A$FileName
## [1] "./data/SRR446027_1.fastq.gz"
## 
## $M1A$SampleName
## [1] "M1A"
## 
## $M1A$Factor
## [1] "M1"
## 
## $M1A$SampleLong
## [1] "Mock.1h.A"
## 
## $M1A$Experiment
## [1] 1
## 
## $M1A$Date
## [1] "23-Mar-2012"
```

```r
targets.as.df(targets(WF))[1:4, 1:4]
```

```
##                      FileName SampleName Factor SampleLong
## 1 ./data/SRR446027_1.fastq.gz        M1A     M1  Mock.1h.A
## 2 ./data/SRR446028_1.fastq.gz        M1B     M1  Mock.1h.B
## 3 ./data/SRR446029_1.fastq.gz        A1A     A1   Avr.1h.A
## 4 ./data/SRR446030_1.fastq.gz        A1B     A1   Avr.1h.B
```

```r
output(WF)[1]
```

```
## $M1A
## $M1A$`hisat2-mapping-se`
## [1] "./results/M1A.sam"
```

```r
cwlfiles(WF)
```

```
## $cwl
## [1] "/home/tgirke/R/x86_64-pc-linux-gnu-library/4.0/systemPipeR/extdata/cwl/hisat2/hisat2-se/hisat2-mapping-se.cwl"
## 
## $yml
## [1] "/home/tgirke/R/x86_64-pc-linux-gnu-library/4.0/systemPipeR/extdata/cwl/hisat2/hisat2-se/hisat2-mapping-se.yml"
## 
## $steps
## [1] "hisat2-mapping-se"
```

```r
inputvars(WF)
```

```
## $FileName
## [1] "_FASTQ_PATH1_"
## 
## $SampleName
## [1] "_SampleName_"
```

In an 'R-centric' rather than a 'CWL-centric' workflow design the connectivity
among workflow steps is established by writing all relevant output with the
*writeTargetsout* function to a new targets file that serves as input to the
next *loadWorkflow* and *renderWF* call. By chaining several _`SYSargs2`_ steps
together one can construct complex workflows involving many sample-level
input/output file operations with any combination of command-line or R-based
software. Alternatively, a CWL-centric workflow design can be used that defines
all/most workflow steps with CWL workflow and parameter files. Due to time and
space restrictions, the CWL-centric approach is not covered by this tutorial. 

### Third-party software tools {#tools}

Current, *systemPipeR* provides the _`param`_ file templates for third-party software tools. Please check the listed software tools. 































































































































































