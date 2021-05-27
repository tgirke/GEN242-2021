---
title: "Automate Creation of CWL Instructions" 
author: "Author: Daniela Cassol, Le Zhang, Thomas Girke"
date: "Last update: 26 May, 2021" 
output:
  BiocStyle::html_document:
    toc_float: true
    code_folding: show
  BiocStyle::pdf_document: default
package: systemPipeR
vignette: |
  %\VignetteEncoding{UTF-8}
  %\VignetteIndexEntry{systemPipeR and CWL}
  %\VignetteEngine{knitr::rmarkdown}
fontsize: 14pt
weight: 18
bibliography: bibtex.bib
---

<!---
- Compile from command-line
Rscript -e "rmarkdown::render('systemPipeR_CWL_SYSargs2.Rmd', c('html_document'), clean=FALSE); knitr::knit('systemPipeR_CWL_SYSargs2.Rmd', tangle=TRUE)"
-->

<div style="text-align: right"> 
Source code downloads: &nbsp; &nbsp;
[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/cmdToCwl/cmdToCwl.Rmd) ] &nbsp; &nbsp; 
[ [.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/cmdToCwl/cmdToCwl.R) ]
</div>




## Introduction

A central concept for designing workflows within the `systemPipeR` environment is 
the use of workflow management containers. `systemPipeR` adopted the widely used community standard [Common Workflow Language](https://www.commonwl.org/) (CWL) 
[@Amstutz2016-ka] for describing analysis workflows in a generic and reproducible 
manner.
Using this community standard in `systemPipeR` has many advantages. For instance, 
the integration of CWL allows running `systemPipeR` workflows from a single 
specification instance either entirely from within R, from various command line 
wrappers (e.g., cwl-runner) or from other languages (, e.g., Bash or Python). 
`systemPipeR` includes support for both command line and R/Bioconductor software 
as well as resources for containerization, parallel evaluations on computer 
clusters along with the automated generation of interactive analysis reports.

An important feature of `systemPipeR's` CWL interface is that it provides two 
options to run command line tools and workflows based on CWL. 
First, one can run CWL in its native way via an R-based wrapper utility for 
`cwl-runner` or `cwl-tools` (CWL-based approach). Second, one can run workflows 
using CWL's command line and workflow instructions from within R (R-based approach). 
In the latter case the same CWL workflow definition files (e.g. *.cwl* and *.yml*) 
are used but rendered and executed entirely with R functions defined by `systemPipeR`, 
and thus use CWL mainly as a command line and workflow definition format rather 
than software to run workflows. In this regard `systemPipeR` also provides several 
convenience functions that are useful for designing and debugging workflows, 
such as a command line rendering function to retrieve the exact command line 
strings for each data set and processing step prior to running a command line.

This overview introduces how CWL describes command line tools and how to connect 
them to create workflows. In addition, we will demonstrate how the workflow can
be easily scalable with `systemPipeR.`


## Load package 

Important: this tutorial uses several new functions that are currently
only available in the development version of the `systemPipeR` package that is installed under R-4.0.3 and R-4.1.0
on the HPCC cluster at UCR.



## CWL command line specifications

CWL command line specifications are written in [YAML](http://yaml.org/) format.

In CWL, files with the extension `.cwl` define the parameters of a chosen 
command line step or workflow, while files with the extension `.yml` define 
the input variables of command line steps. 

Let's explore the `.cwl` file: 











































