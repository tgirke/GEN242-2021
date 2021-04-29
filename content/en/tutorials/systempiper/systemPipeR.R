## pre code {

## white-space: pre !important;

## overflow-x: scroll !important;

## word-break: keep-all !important;

## word-wrap: initial !important;

## }


## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()
options(width=80, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")), 
    tidy.opts=list(width.cutoff=80), tidy=TRUE)


## ----setup, echo=FALSE, message=FALSE, warning=FALSE--------------------------
suppressPackageStartupMessages({
    library(systemPipeR)
    library(BiocParallel)
    library(Biostrings)
    library(Rsamtools)
    library(GenomicRanges)
    library(ggplot2)
    library(GenomicAlignments)
    library(ShortRead)
    library(ape)
    library(batchtools)
    library(magrittr)
})


## ----install, eval=FALSE------------------------------------------------------
## if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
## BiocManager::install("systemPipeR")
## BiocManager::install("systemPipeRdata")


## ----documentation, eval=FALSE------------------------------------------------
## library("systemPipeR") # Loads the package
## library(help="systemPipeR") # Lists package info
## vignette("systemPipeR") # Opens vignette


## ----targetsSE, eval=TRUE-----------------------------------------------------
library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR") 
read.delim(targetspath, comment.char = "#")[1:4,]


## ----targetsPE, eval=TRUE-----------------------------------------------------
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")[1:2,1:6]


## ----comment_lines, echo=TRUE-------------------------------------------------
readLines(targetspath)[1:4]


## ----targetscomp, eval=TRUE---------------------------------------------------
readComp(file=targetspath, format="vector", delim="-")


## ----CWL_structure, echo = FALSE, eval=FALSE----------------------------------
## hisat2.cwl <- system.file("extdata", "cwl/hisat2/hisat2-se/hisat2-mapping-se.cwl", package="systemPipeR")
## yaml::read_yaml(hisat2.cwl)


## ----yaml_structure, echo = FALSE, eval=FALSE---------------------------------
## hisat2.yml <- system.file("extdata", "cwl/hisat2/hisat2-se/hisat2-mapping-se.yml", package="systemPipeR")
## yaml::read_yaml(hisat2.yml)


## ----SYSargs2_structure, eval=TRUE--------------------------------------------
library(systemPipeR)
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package="systemPipeR")
WF <- loadWF(targets=targets, wf_file="hisat2-mapping-se.cwl",
                   input_file="hisat2-mapping-se.yml",
                   dir_path=dir_path)

WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))


## ----names_WF, eval=TRUE------------------------------------------------------
names(WF)


## ----cmdlist, eval=TRUE-------------------------------------------------------
cmdlist(WF)[1]


## ----output_WF, eval=TRUE-----------------------------------------------------
output(WF)[1]
modules(WF)
targets(WF)[1]
targets.as.df(targets(WF))[1:4,1:4]
output(WF)[1]
cwlfiles(WF)
inputvars(WF)


## ----table_tools, echo=FALSE, message=FALSE-----------------------------------
library(magrittr)
# library(kableExtra)
SPR_software <- system.file("extdata", "SPR_software.csv", package="systemPipeR")
#SPR_software <- "SPR_software.csv"
software <- read.delim(SPR_software, sep=",", comment.char = "#")
colors <- colorRampPalette((c("darkseagreen", "indianred1")))(length(unique(software$Category)))
id <- as.numeric(c((unique(software$Category))))
software %>%
  dplyr::mutate(Step = kableExtra::cell_spec(Step, color = "white", bold = TRUE,
    background = factor(Category, id, colors)
  )) %>%
   dplyr::select(Tool, Description, Step) %>%
  dplyr::arrange(Tool) %>% 
  kableExtra::kable(escape = FALSE, align = "c", col.names = c("Tool Name",	"Description", "Step")) %>%
  kableExtra::kable_styling(c("striped", "hover", "condensed"), full_width = TRUE) %>%
  kableExtra::scroll_box(width = "100%", height = "500px")


## ----test_tool_path, eval=FALSE-----------------------------------------------
## tryCL(command="grep")


## ----workflow, eval=FALSE-----------------------------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="rnaseq")
## setwd("rnaseq")


## ----test_software, eval=FALSE------------------------------------------------
## tryCL(command="hisat2") ## "All set up, proceed!"


## ----initwf, eval=FALSE-------------------------------------------------------
## getwd() ## rnaseq
## script <- "systemPipeRNAseq.Rmd"
## targetspath <- "targets.txt"
## sysargslist <- initWF(script = script, targets = targetspath)


## ----config, eval=FALSE-------------------------------------------------------
## sysargslist <- configWF(x=sysargslist, input_steps = "1:3")
## sysargslist <- runWF(sysargslist = sysargslist, steps = "1:2")


## ----pipes, eval=FALSE--------------------------------------------------------
## sysargslist <- initWF(script ="systemPipeRNAseq.Rmd", overwrite = TRUE) %>%
##     configWF(input_steps = "1:3") %>%
##     runWF(steps = "1:2")


## ----closeR, eval=FALSE-------------------------------------------------------
## q("no") # closes R session on head node


## srun --x11 --partition=short --mem=2gb --cpus-per-task 4 --ntasks 1 --time 2:00:00 --pty bash -l

## module load R/4.0.3

## R


## ----r_environment, eval=FALSE------------------------------------------------
## system("hostname") # should return name of a compute node starting with i or c
## getwd() # checks current working directory of R session
## dir() # returns content of current working directory


## ----clusterRun, eval=FALSE---------------------------------------------------
## library(batchtools)
## targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
## dir_path <- system.file("extdata/cwl/hisat2/hisat2-pe", package="systemPipeR")
## args <- loadWorkflow(targets=targetspath, wf_file="hisat2-mapping-pe.cwl",
##                      input_file="hisat2-mapping-pe.yml", dir_path=dir_path)
## args <- renderWF(args, inputvars=c(FileName1="_FASTQ_PATH1_", FileName2="_FASTQ_PATH2_", SampleName="_SampleName_"))
## resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
## reg <- clusterRun(args, FUN = runCommandline, more.args = list(args=args, make_bam=TRUE, dir=FALSE),
##                   conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl",
##                   Njobs=18, runid="01", resourceList=resources)
## getStatus(reg=reg)
## waitForJobs(reg=reg)


## ----genRna_workflow_single, eval=FALSE---------------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="rnaseq")
## setwd("rnaseq")


## ----genChip_workflow_single, eval=FALSE--------------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="chipseq")
## setwd("chipseq")


## ----genVar_workflow_single, eval=FALSE---------------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="varseq")
## setwd("varseq")


## ----genRibo_workflow_single, eval=FALSE--------------------------------------
## library(systemPipeRdata)
## genWorkenvir(workflow="riboseq")
## setwd("riboseq")


## ----sessionInfo--------------------------------------------------------------
sessionInfo()

