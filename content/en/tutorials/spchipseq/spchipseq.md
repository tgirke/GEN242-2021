---
title: "ChIP-Seq Workflow Template" 
author: "Author: First Last"
date: "Last update: 06 June, 2021" 
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

<div style="text-align: right">

Source code downloads:    
\[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242/main/static/custom/spWFtemplates/spchipseq.Rmd) \]    
\[ [.html](https://girke.bioinformatics.ucr.edu/GEN242/custom/spWFtemplates/spchipseq.html) \]    
\[ [old version .Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/spchipseq/back_spchipseq_rmd) \]

</div>

## Introduction

The following analyzes the ChIP-Seq data from Kaufman et al. (2010) using
for peak calling MACS2 where the uninduced sample serves as input (reference).
The details about all download steps are provided [here](https://girke.bioinformatics.ucr.edu/GEN242/assignments/projects/project_data/).

Users want to extend this section to provide all background information relevant for this
ChIP-Seq project.

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
srun --x11 --partition=gen242 --mem=20gb --cpus-per-task 8 --ntasks 1 --time 20:00:00 --pty bash -l
module unload R; module load R/4.0.3_gcc-8.3.0
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

To apply workflows to custom data, the user needs to modify the *`targets`* file and if
necessary update the corresponding parameter (*`.cwl`* and *`.yml`*) files.
A collection of pre-generated *`.cwl`* and *`.yml`* files are provided in the *`param/cwl`* subdirectory
of each workflow template. They are also viewable in the GitHub repository of *`systemPipeRdata`* ([see
here](https://github.com/tgirke/systemPipeRdata/tree/master/inst/extdata/param/cwl)).
For more information of the structure of the *targets* file, please consult the documentation
[here](http://www.bioconductor.org/packages/release/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html#25_structure_of_targets_file). More details about the new parameter files from systemPipeR can be found [here](http://www.bioconductor.org/packages/release/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html#26_structure_of_the_new_param_files_and_construct_sysargs2_container).

### Import custom functions

Custem functions for the challenge projects can be imported with the source
command from a local R script (here [challengeProject\_Fct.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/spchipseq/challengeProject_Fct.R)). Skip this step if such a
script is not available. Alternatively, these functions can be loaded from a
custom R package.

``` r
source("challengeProject_Fct.R")
```

### Experiment definition provided by `targets` file

The `targets` file defines all FASTQ files and sample comparisons of the analysis workflow.
If needed the tab separated (TSV) version of this file can be downloaded from [here](https://github.com/tgirke/GEN242/tree/main/content/en/assignments/Projects/targets_files)
and the corresponding Google Sheet is [here](https://docs.google.com/spreadsheets/d/1w9V3JDOsXR8qW_qNqoXO0E0UR_Oev4k7-o2y6NMDTt8/edit#gid=472150521).

``` r
targetspath <- "targets_chipseq.txt"
targets <- read.delim(targetspath, comment.char = "#")
knitr::kable(targets)
```

| FileName                     | SampleName | Factor | SampleLong       | Experiment | Date      | SampleReference |
|:-----------------------------|:-----------|:-------|:-----------------|-----------:|:----------|:----------------|
| ./data/SRR038845\_1.fastq.gz | AP1\_1     | AP1    | APETALA1 Induced |          1 | 23-Mar-12 |                 |
| ./data/SRR038846\_1.fastq.gz | AP1\_2A    | AP1    | APETALA1 Induced |          1 | 23-Mar-12 |                 |
| ./data/SRR038847\_1.fastq.gz | AP1\_2B    | AP1    | APETALA1 Induced |          1 | 23-Mar-12 |                 |
| ./data/SRR038848\_1.fastq.gz | C\_1A      | C      | Control Mock     |          1 | 23-Mar-12 | AP1\_1          |
| ./data/SRR038849\_1.fastq.gz | C\_1B      | C      | Control Mock     |          1 | 23-Mar-12 | AP1\_1          |
| ./data/SRR038850\_1.fastq.gz | C\_2A      | C      | Control Mock     |          1 | 23-Mar-12 | AP1\_2A         |
| ./data/SRR038851\_1.fastq.gz | C\_2B      | C      | Control Mock     |          1 | 23-Mar-12 | AP1\_2B         |

## Read preprocessing

### Read quality filtering and trimming

The following example shows how one can design a custom read
preprocessing function using utilities provided by the `ShortRead` package, and then
apply it with `preprocessReads` in batch mode to all FASTQ samples referenced in the
corresponding `SYSargs2` instance (`trim` object below). More detailed information on
read preprocessing is provided in `systemPipeR's` main vignette.

First, we construct *`SYSargs2`* object from *`cwl`* and *`yml`* param and *`targets`* files.

``` r
trim <- loadWF(targets = targetspath, wf_file = "trim-se.cwl", 
    input_file = "trim-se.yml", dir_path = "param/cwl/preprocessReads/trim-se")
trim <- renderWF(trim, inputvars = c(FileName = "_FASTQ_PATH_", 
    SampleName = "_SampleName_"))
trim
output(trim)[1:2]
```

Next, we execute the code for trimming all the raw data. Note, the quality
settings are relatively relaxed in this step (Phred score of at least 10 and tolerating
two Ns per read), because this data is from a time when the quality of
Illumina sequencing was still low. Setting the quality parameter
more stringent would remove too many reads, which would negatively impact the
read coverage required for the downstream peak calling.

``` r
filterFct <- function(fq, cutoff = 10, Nexceptions = 2) {
    qcount <- rowSums(as(quality(fq), "matrix") <= cutoff, na.rm = TRUE)
    fq[qcount <= Nexceptions]
    # Retains reads where Phred scores are >= cutoff with N
    # exceptions
}
preprocessReads(args = trim, Fct = "filterFct(fq, cutoff=10, Nexceptions=2)", 
    batchsize = 1e+05)
writeTargetsout(x = trim, file = "targets_chip_trim.txt", step = 1, 
    new_col = c("FileName"), new_col_output_index = 1, overwrite = TRUE)
```

### FASTQ quality report

The following `seeFastq` and `seeFastqPlot` functions generate and plot a
series of useful quality statistics for a set of FASTQ files including per
cycle quality box plots, base proportions, base-level quality trends,
relative k-mer diversity, length and occurrence distribution of reads, number
of reads above quality cutoffs and mean quality distribution. The results are
written to a PDF file named `fastqReport.png`. Parallelization of FASTQ
quality report via scheduler (*e.g.* Slurm) across several compute nodes.

``` r
library(BiocParallel)
library(batchtools)
f <- function(x) {
    library(systemPipeR)
    targets <- "targets_chip_trim.txt"
    dir_path <- "param/cwl/preprocessReads/trim-se"
    trim <- loadWorkflow(targets = targets, wf_file = "trim-se.cwl", 
        input_file = "trim-se.yml", dir_path = dir_path)
    trim <- renderWF(trim, inputvars = c(FileName = "_FASTQ_PATH_", 
        SampleName = "_SampleName_"))
    outfile <- subsetWF(trim, slot = "output", subset = 1, index = 1)
    test = seeFastq(fastq = outfile[x], batchsize = 1e+05, klength = 8)
}

resources <- list(walltime = 120, ntasks = 1, ncpus = 4, memory = 1024)
param <- BatchtoolsParam(workers = 4, cluster = "slurm", template = "batchtools.slurm.tmpl", 
    resources = resources)
fqlist <- bplapply(seq(along = trim), f, BPPARAM = param)

png("./results/fastqReport.png", height = 18 * 96, width = 4 * 
    96 * length(fqlist))
seeFastqPlot(unlist(fqlist, recursive = FALSE))
dev.off()
```

![](../results/fastqReport.png)

<div align="center">

Figure 1: FASTQ quality report for 7 samples.

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
idx <- loadWorkflow(targets = NULL, wf_file = "bowtie2-index.cwl", 
    input_file = "bowtie2-index.yml", dir_path = "param/cwl/bowtie2/bowtie2-idx")
idx <- renderWF(idx)
idx
cmdlist(idx)

## Run in single machine
runCommandline(idx, make_bam = FALSE)
```

The following submits 7 alignment jobs via a scheduler to a computer cluster.

``` r
targets <- "targets_chip_trim.txt"
dir_path <- "param/cwl/bowtie2/bowtie2-se"
args <- loadWF(targets = targets, wf_file = "bowtie2-mapping-se.cwl", 
    input_file = "bowtie2-mapping-se.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))
args
cmdlist(args)[1:2]
output(args)[1:2]
```

``` r
moduleload(modules(args))  # Skip if a module system is not used
resources <- list(walltime = 120, ntasks = 1, ncpus = 4, memory = 1024)
reg <- clusterRun(args, FUN = runCommandline, more.args = list(args = args, 
    dir = FALSE), conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl", 
    Njobs = 7, runid = "01", resourceList = resources)
getStatus(reg = reg)
waitForJobs(reg = reg)
args <- output_update(args, dir = FALSE, replace = TRUE, extension = c(".sam", 
    ".bam"))  ## Updates the output(args) to the right location in the subfolders
output(args)
```

Alternatively, one can run the alignments sequentially on a single system.
Note: this step is not used here!

``` r
# args <- runCommandline(args, force=FALSE)
```

Check whether all BAM files and the corresponding new targets have been
created.

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
read_statsDF <- alignStats(args = args, output_index = 1, subset = "FileName")
write.table(read_statsDF, "results/alignStats.xls", row.names = FALSE, 
    quote = FALSE, sep = "\t")
read.delim("results/alignStats.xls")
```

### Create symbolic links for viewing BAM files in IGV

The `symLink2bam` function creates symbolic links to view the BAM alignment
files in a genome browser such as IGV without moving these large files to a
local system. The corresponding URLs are written to a file with a path
specified under `urlfile`, here `IGVurl.txt`. Please replace the directory
and the user name. The following parameter settings will create a
subdirectory under `~/.html` called `somedir` of the user account. The user
name under `urlbase`, here `ttest`, needs to be changed to the corresponding
user name of the person running this function.

``` r
symLink2bam(sysargs = args, htmldir = c("~/.html/", "somedir/"), 
    urlbase = "http://cluster.hpcc.ucr.edu/~<username>/", urlfile = "./results/IGVurl.txt")
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
dir_path <- "param/cwl/mergeBamByFactor"
args <- loadWF(targets = "targets_bam.txt", wf_file = "merge-bam.cwl", 
    input_file = "merge-bam.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName = "_BAM_PATH_", 
    SampleName = "_SampleName_"))
args_merge <- mergeBamByFactor(args = args, overwrite = TRUE)
writeTargetsout(x = args_merge, file = "targets_mergeBamByFactor.txt", 
    step = 1, new_col = "FileName", new_col_output_index = 1, 
    overwrite = TRUE, remove = TRUE)
```

### Peak calling with input/reference sample

MACS2 can perform peak calling on ChIP-Seq data with and without input
samples (Zhang et al. 2008).

The following performs peak calling with input sample. The input sample
can be most conveniently specified in the `SampleReference` column of the
initial `targets` file. The `writeTargetsRef` function uses this
information to create a `targets` file intermediate for running MACS2
with the corresponding input sample(s).

``` r
writeTargetsRef(infile = "targets_mergeBamByFactor.txt", outfile = "targets_bam_ref.txt", 
    silent = FALSE, overwrite = TRUE)
dir_path <- "param/cwl/MACS2/MACS2-input"
args_input <- loadWF(targets = "targets_bam_ref.txt", wf_file = "macs2-input.cwl", 
    input_file = "macs2.yml", dir_path = dir_path)
args_input <- renderWF(args_input, inputvars = c(FileName1 = "_FASTQ_PATH2_", 
    FileName2 = "_FASTQ_PATH1_", SampleName = "_SampleName_"))
cmdlist(args_input)[1]
## Run MACS2
args_input <- runCommandline(args_input, make_bam = FALSE, force = TRUE)
outpaths_input <- subsetWF(args_input, slot = "output", subset = 1, 
    index = 1)
file.exists(outpaths_input)
writeTargetsout(x = args_input, file = "targets_macs_input.txt", 
    step = 1, new_col = "FileName", new_col_output_index = 1, 
    overwrite = TRUE)
```

The peak calling results from MACS2 are written for each sample to the
`results` directory. They are named after the corresponding reference sample
with extensions used by MACS2.

## Annotate peaks with genomic context

### Annotation with `ChIPseeker` package

To annotate the identified peaks with genomic context information
one can use the `ChIPpeakAnno` or `ChIPseeker` package (Zhu et al. 2010; Yu, Wang, and He 2015).
The following code uses the `ChIPseeker` package for annotating the peaks.

``` r
library(ChIPseeker)
library(GenomicFeatures)
dir_path <- "param/cwl/annotate_peaks"
args <- loadWF(targets = "targets_macs_input.txt", wf_file = "annotate-peaks.cwl", 
    input_file = "annotate-peaks.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))

txdb <- makeTxDbFromGFF(file = "data/tair10.gff", format = "gff", 
    dataSource = "TAIR", organism = "Arabidopsis thaliana")
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

The peak annotation results are written to the `results` directory.
The files are named after the corresponding peak files with extensions
specified in the `annotate_peaks.param` file, here `*.peaks.annotated.xls`.

## Count reads overlapping peaks

The `countRangeset` function is a convenience wrapper to perform read counting
iteratively over serveral range sets, here peak range sets. Internally,
the read counting is performed with the `summarizeOverlaps` function from the
`GenomicAlignments` package. The resulting count tables are directly saved to
files, one for each peak set.

``` r
library(GenomicRanges)
dir_path <- "param/cwl/count_rangesets"
args <- loadWF(targets = "targets_macs_input.txt", wf_file = "count_rangesets.cwl", 
    input_file = "count_rangesets.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))

## Bam Files
targets <- "targets_chip_trim.txt"
dir_path <- "param/cwl/bowtie2/bowtie2-se"
args_bam <- loadWF(targets = targets, wf_file = "bowtie2-mapping-se.cwl", 
    input_file = "bowtie2-mapping-se.yml", dir_path = dir_path)
args_bam <- renderWF(args_bam, inputvars = c(FileName = "_FASTQ_PATH1_", 
    SampleName = "_SampleName_"))
args_bam <- output_update(args_bam, dir = FALSE, replace = TRUE, 
    extension = c(".sam", ".bam"))
outpaths <- subsetWF(args_bam, slot = "output", subset = 1, index = 1)

register(MulticoreParam(workers = 3))
bfl <- BamFileList(outpaths, yieldSize = 50000, index = character())
countDFnames <- countRangeset(bfl, args, mode = "Union", ignore.strand = TRUE)
writeTargetsout(x = args, file = "targets_countDF.txt", step = 1, 
    new_col = "FileName", new_col_output_index = 1, overwrite = TRUE)
```

Shows count table generated in previous step (`C_1A_peaks.countDF.xls`).
To avoid slowdowns of the load time of this page, ony 200 rows of the source
table are imported into the below `datatable` view .

``` r
countDF <- read.delim("results/C_1A_peaks.countDF.xls")[1:200, 
    ]
colnames(countDF)[1] <- "PeakIDs"
library(DT)
datatable(countDF)
```

<div id="htmlwidget-1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-1">{"x":{"filter":"none","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199","200"],["Chr1_10000123-10001504","Chr1_10015050-10015329","Chr1_10044084-10044292","Chr1_10053060-10053403","Chr1_10065376-10065801","Chr1_10077918-10078299","Chr1_10111557-10111830","Chr1_10137133-10137852","Chr1_10139865-10140426","Chr1_10141346-10141599","Chr1_10142995-10143536","Chr1_10143619-10144022","Chr1_10145791-10146202","Chr1_10146333-10147178","Chr1_10186868-10187102","Chr1_10187659-10187893","Chr1_10203888-10204202","Chr1_10219688-10219985","Chr1_10228577-10229094","Chr1_10229262-10229710","Chr1_10242962-10243288","Chr1_10244590-10244954","Chr1_10247304-10247528","Chr1_10251234-10251447","Chr1_10263957-10264248","Chr1_102694-103029","Chr1_10294478-10294800","Chr1_10295972-10296231","Chr1_10353709-10354318","Chr1_10364919-10365439","Chr1_10369960-10370330","Chr1_10375634-10375976","Chr1_10376847-10377070","Chr1_10489119-10489357","Chr1_10489629-10489865","Chr1_1049818-1050057","Chr1_10504608-10505267","Chr1_10505986-10506492","Chr1_10519205-10519487","Chr1_10520174-10520649","Chr1_10520971-10521455","Chr1_10530766-10530996","Chr1_10535179-10535644","Chr1_10537198-10537646","Chr1_10538872-10539144","Chr1_10541237-10541831","Chr1_10542202-10542449","Chr1_10553478-10553874","Chr1_10562570-10563068","Chr1_10563816-10565000","Chr1_10565793-10566457","Chr1_10567369-10567945","Chr1_10569384-10569620","Chr1_105715-106037","Chr1_10581383-10581685","Chr1_10622431-10622757","Chr1_10631994-10632273","Chr1_10649137-10649367","Chr1_10668660-10668958","Chr1_10669009-10669489","Chr1_10676832-10677414","Chr1_10678636-10679743","Chr1_10679884-10680612","Chr1_10683747-10684382","Chr1_10684727-10685161","Chr1_10685376-10685989","Chr1_10691341-10691562","Chr1_10691767-10692130","Chr1_10692583-10693433","Chr1_10693954-10694617","Chr1_10695463-10696565","Chr1_10701975-10702395","Chr1_10704093-10704316","Chr1_10713567-10714255","Chr1_10715301-10715933","Chr1_1072348-1072586","Chr1_10726169-10726568","Chr1_1072802-1073152","Chr1_10737905-10738222","Chr1_10793236-10793520","Chr1_10794884-10795869","Chr1_10799998-10800389","Chr1_10800942-10801268","Chr1_10801874-10802349","Chr1_10803701-10803926","Chr1_1084359-1084961","Chr1_10845312-10845874","Chr1_10864966-10865601","Chr1_10871974-10872308","Chr1_10874129-10874328","Chr1_10874476-10875185","Chr1_10891639-10891997","Chr1_10909961-10910163","Chr1_10953463-10953804","Chr1_10954240-10954628","Chr1_10956271-10956623","Chr1_10958841-10959244","Chr1_10959439-10960076","Chr1_10972812-10973048","Chr1_1103200-1103416","Chr1_11033716-11034097","Chr1_11059998-11060940","Chr1_11061272-11061800","Chr1_11062322-11062750","Chr1_11067587-11068014","Chr1_11071124-11071500","Chr1_11084496-11084695","Chr1_11128-11443","Chr1_11155158-11155410","Chr1_11196849-11197415","Chr1_11200138-11200385","Chr1_11210773-11211568","Chr1_11212068-11212348","Chr1_11212560-11212978","Chr1_1121262-1121569","Chr1_11222722-11222952","Chr1_11227689-11228070","Chr1_11231069-11231415","Chr1_1129384-1129772","Chr1_1134217-1134971","Chr1_11343316-11343530","Chr1_11350683-11350893","Chr1_11374190-11375132","Chr1_11377612-11378008","Chr1_11385968-11386241","Chr1_11386635-11386872","Chr1_11388323-11388599","Chr1_11421622-11422073","Chr1_11428546-11428802","Chr1_11429167-11429366","Chr1_11429941-11430188","Chr1_11430639-11430926","Chr1_11441925-11442147","Chr1_11446002-11446822","Chr1_11447290-11447887","Chr1_115641-116419","Chr1_1159091-1159510","Chr1_11594952-11595190","Chr1_11598295-11598860","Chr1_11600614-11600983","Chr1_11601139-11601823","Chr1_11605622-11605926","Chr1_11620054-11620600","Chr1_11621059-11621447","Chr1_11621523-11623089","Chr1_11624914-11625625","Chr1_11626261-11626946","Chr1_11628379-11628622","Chr1_11628853-11629244","Chr1_11630475-11631803","Chr1_11662337-11662575","Chr1_1172920-1173129","Chr1_1173195-1173526","Chr1_11763531-11763800","Chr1_11770202-11770807","Chr1_11771164-11771564","Chr1_11798751-11798970","Chr1_11801046-11801322","Chr1_11804660-11805026","Chr1_11807363-11807831","Chr1_11808375-11808791","Chr1_11809718-11810175","Chr1_11829516-11829982","Chr1_11830901-11831419","Chr1_11874425-11874646","Chr1_1188104-1188901","Chr1_11925063-11925373","Chr1_11927715-11928115","Chr1_11928466-11928909","Chr1_11978528-11979276","Chr1_1197874-1198387","Chr1_12004740-12004982","Chr1_12005493-12006453","Chr1_12026542-12027136","Chr1_12050647-12051422","Chr1_12051618-12051822","Chr1_12054702-12055258","Chr1_12055957-12056342","Chr1_12056485-12057547","Chr1_12059598-12059861","Chr1_12064129-12064328","Chr1_12070355-12070698","Chr1_12123473-12123764","Chr1_12132602-12132955","Chr1_12177510-12177803","Chr1_12180624-12180825","Chr1_12254744-12255260","Chr1_1233120-1233354","Chr1_1235877-1236728","Chr1_12458546-12458848","Chr1_12474908-12475186","Chr1_128465-128823","Chr1_12916443-12916688","Chr1_12986436-12986729","Chr1_12986791-12987139","Chr1_13029530-13030312","Chr1_1303036-1303264","Chr1_130507-131207","Chr1_1306355-1306663","Chr1_1307140-1307396"],[384,19,54,31,89,82,36,95,135,29,44,44,59,134,27,14,50,49,193,75,36,29,22,24,33,21,19,16,108,125,38,17,20,20,28,36,135,53,23,168,101,19,67,26,31,55,40,39,40,290,229,40,21,39,28,13,36,16,30,40,96,271,81,63,48,106,18,40,207,86,121,85,19,158,100,38,72,40,30,16,277,51,54,97,22,88,120,77,19,18,137,64,23,44,53,47,33,150,22,29,39,107,161,49,94,23,8,60,53,91,27,58,32,69,33,26,26,47,49,222,15,32,191,35,17,26,43,71,17,16,13,39,6,120,80,93,116,17,86,41,178,36,125,48,178,59,92,23,41,284,80,17,38,23,104,104,12,13,35,72,80,82,49,98,15,100,16,36,32,251,72,25,143,127,177,16,94,37,219,22,12,80,55,37,30,19,102,9,268,26,13,74,31,23,37,304,30,93,40,16],[571,50,31,71,122,115,59,201,191,38,102,72,120,296,32,44,59,56,357,129,54,62,36,35,47,64,66,43,213,149,69,48,33,45,37,19,281,95,63,301,180,39,109,83,58,86,50,80,100,625,371,100,35,69,46,78,46,44,68,106,148,447,152,152,123,205,38,71,347,159,264,96,31,196,144,37,150,42,47,37,434,93,79,135,37,132,186,137,71,38,236,85,33,54,86,63,76,243,53,29,114,319,248,135,277,73,30,92,37,165,62,157,91,138,61,34,81,82,63,330,31,32,295,78,58,48,52,131,46,26,43,54,36,167,112,154,128,31,89,49,199,60,165,90,399,162,177,29,84,485,141,32,70,41,152,139,39,46,75,127,120,169,157,141,35,199,74,82,114,343,164,44,237,228,284,28,164,61,433,37,43,99,44,75,69,33,160,34,285,59,62,153,37,52,45,481,32,281,53,38],[282,26,19,32,89,46,21,90,95,26,58,30,58,149,25,20,35,24,158,80,26,22,22,21,23,26,49,17,96,83,24,28,13,18,16,15,107,40,32,152,74,21,49,36,27,36,19,31,45,348,151,54,14,31,15,32,31,18,21,41,83,222,82,62,40,102,10,30,143,98,119,38,12,92,66,24,94,21,20,27,189,41,52,70,19,41,94,59,33,19,132,26,11,24,29,25,39,119,13,14,53,155,110,55,110,50,14,42,15,77,21,91,27,58,37,24,42,44,31,151,15,25,142,32,22,20,18,67,28,15,22,21,19,75,46,58,43,14,45,27,90,25,51,44,191,76,83,19,48,248,67,21,35,16,68,81,14,19,42,57,78,67,74,58,9,76,24,23,45,169,72,28,117,96,141,14,67,29,195,20,17,57,25,38,29,13,78,26,118,16,25,71,15,28,16,200,16,129,27,23],[56,1,4,6,1,12,3,15,7,5,10,34,1,1,2,0,4,13,1,17,1,0,0,0,4,3,1,1,15,15,10,9,2,5,3,1,6,19,0,3,13,6,11,0,8,2,18,4,5,9,7,7,2,10,0,4,5,0,10,7,0,11,3,5,1,17,0,10,12,0,25,3,0,4,19,10,4,15,5,0,12,2,1,2,2,60,16,2,5,14,2,3,2,12,12,11,3,6,4,1,2,4,9,0,24,3,14,1,7,0,6,8,2,9,4,20,1,10,6,1,1,7,10,34,7,4,3,4,6,2,0,22,1,6,13,1,23,0,33,9,1,17,61,2,21,1,19,0,9,3,5,5,10,5,0,16,3,34,6,14,0,0,5,25,0,8,11,1,1,43,9,0,10,15,10,0,22,16,14,4,0,1,14,2,19,21,3,1,34,5,11,67,29,5,15,8,1,9,1,3],[43,3,3,6,6,4,8,8,8,8,10,12,1,4,3,1,9,5,2,5,7,1,6,1,6,3,4,1,9,10,8,6,2,8,1,0,6,22,2,3,18,4,10,3,8,7,9,8,5,32,17,7,7,10,0,3,4,0,8,3,0,21,23,9,1,20,1,11,12,1,33,7,2,9,13,8,9,16,3,2,9,4,8,8,2,19,14,11,5,5,17,1,0,7,6,5,5,9,2,4,6,17,12,3,16,1,8,8,5,1,7,2,3,3,9,16,0,7,7,19,9,2,23,19,4,0,7,14,2,3,0,15,0,6,24,4,19,2,10,7,8,12,17,0,39,14,11,4,3,6,22,4,12,7,1,11,5,10,11,15,0,1,9,17,2,13,11,0,4,26,12,0,21,18,17,0,15,5,19,1,1,5,6,8,11,6,1,1,22,9,5,27,7,1,9,14,0,6,1,7],[104,12,14,25,18,10,17,20,20,5,21,11,11,15,8,10,13,16,32,16,12,10,9,7,18,11,7,4,21,27,12,8,2,6,10,4,22,16,13,30,11,17,13,15,13,16,6,14,14,33,36,11,7,2,6,11,12,12,16,19,13,41,21,20,10,25,6,22,30,15,33,15,14,17,19,9,35,10,14,9,32,22,8,8,8,15,18,17,20,14,22,15,5,7,13,5,12,11,6,7,14,41,22,20,46,8,11,14,7,16,5,9,12,10,9,5,12,19,9,31,8,12,48,15,11,4,5,36,4,4,12,5,13,20,20,18,28,6,13,15,31,10,12,7,41,23,26,5,22,36,43,5,13,2,28,29,7,5,16,18,18,18,24,12,10,30,9,9,13,55,26,3,46,37,30,6,29,5,21,11,7,21,7,13,9,3,24,5,36,23,16,25,4,6,11,31,12,43,18,4],[55,5,3,4,9,4,9,15,17,6,4,8,10,14,4,4,7,8,19,5,4,12,5,5,2,7,5,4,12,17,6,6,4,1,3,5,16,5,2,19,10,4,7,6,2,5,8,7,9,14,24,3,7,3,3,9,7,5,7,11,9,22,11,2,8,16,4,10,14,7,13,10,2,10,11,4,30,5,8,3,20,3,6,7,5,9,5,8,10,4,12,9,10,0,9,2,7,18,2,6,11,26,11,15,23,9,4,10,6,9,6,9,4,9,4,1,3,8,7,13,3,8,20,10,5,5,6,12,3,1,5,9,3,8,11,10,17,1,4,5,18,4,12,6,25,13,11,0,10,19,22,5,10,4,19,17,2,8,6,18,7,9,7,7,4,8,3,2,4,23,11,1,24,13,14,4,16,7,12,5,2,9,3,11,7,4,20,6,13,6,14,14,2,4,9,17,4,32,3,4]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>PeakIDs<\/th>\n      <th>AP1_1<\/th>\n      <th>AP1_2A<\/th>\n      <th>AP1_2B<\/th>\n      <th>C_1A<\/th>\n      <th>C_1B<\/th>\n      <th>C_2A<\/th>\n      <th>C_2B<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

## Differential binding analysis

The `runDiff` function performs differential binding analysis in batch mode for
several count tables using `edgeR` or `DESeq2` (Robinson, McCarthy, and Smyth 2010; Love, Huber, and Anders 2014).
Internally, it calls the functions `run_edgeR` and `run_DESeq2`. It also returns
the filtering results and plots from the downstream `filterDEGs` function using
the fold change and FDR cutoffs provided under the `dbrfilter` argument.

``` r
dir_path <- "param/cwl/rundiff"
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

The following performs GO term enrichment analysis for each annotated peak
set. Note: the following assumes that the GO annotation data exists under
`data/GO/catdb.RData`. If this is not the case then it can be generated with
the instructions from [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/sprnaseq/sprnaseq/#obtain-gene-to-go-mappings).

``` r
dir_path <- "param/cwl/annotate_peaks"
args <- loadWF(targets = "targets_bam_ref.txt", wf_file = "annotate-peaks.cwl", 
    input_file = "annotate-peaks.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
    FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))

args_anno <- loadWF(targets = "targets_macs_input.txt", wf_file = "annotate-peaks.cwl", 
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
write.table(BatchResult, "results/GOBatchAll.xls", row.names = FALSE, 
    quote = FALSE, sep = "\t")
```

Shows GO term enrichment results from previous step. The last gene identifier column (10)
of this table has been excluded in this viewing instance to minimze the complexity of the
result.
To avoid slowdowns of the load time of this page, ony 200 rows of the source
table are imported into the below `datatable` view .

``` r
BatchResult <- read.delim("results/GOBatchAll.xls")[1:200, ]
library(DT)
datatable(BatchResult[, -10], options = list(scrollX = TRUE, 
    autoWidth = TRUE))
```

<div id="htmlwidget-2" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-2">{"x":{"filter":"none","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199","200"],["C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A","C_1A"],[4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903,4903],["GO:0003700","GO:0140110","GO:0098772","GO:0043565","GO:0000976","GO:0001067","GO:1990837","GO:0003677","GO:0003690","GO:0042802","GO:0000981","GO:0046983","GO:0000987","GO:0003676","GO:0000978","GO:1901363","GO:0033612","GO:0097159","GO:0042562","GO:0005515","GO:0000977","GO:0060089","GO:0004672","GO:0019840","GO:0016773","GO:0043178","GO:0004860","GO:0019210","GO:0016705","GO:0015291","GO:0038023","GO:0001653","GO:0016301","GO:0022857","GO:0004861","GO:0030291","GO:0022804","GO:0042803","GO:0000156","GO:0005215","GO:0004674","GO:0019900","GO:0015293","GO:0005102","GO:0004805","GO:0005516","GO:0010427","GO:0080161","GO:0008509","GO:0019207","GO:0010011","GO:0004497","GO:0016757","GO:0016702","GO:0102336","GO:0102337","GO:0102338","GO:0102756","GO:0102682","GO:0019887","GO:0016165","GO:0046423","GO:0016837","GO:0030570","GO:0015103","GO:0001216","no_annot_MF","GO:0010033","GO:0009725","GO:0009719","GO:0042221","GO:0006355","GO:1903506","GO:2001141","GO:0051252","GO:0006351","GO:0097659","GO:0009889","GO:0031326","GO:0032774","GO:0019219","GO:2000112","GO:0010556","GO:0050789","GO:0031323","GO:0050794","GO:1901700","GO:0009755","GO:0080090","GO:0051171","GO:0032870","GO:0065007","GO:0034654","GO:0010468","GO:0050896","GO:0071495","GO:0048367","GO:0070887","GO:1901362","GO:0018130","GO:0019222","GO:0019438","GO:0071310","GO:0060255","GO:0048856","GO:0032502","GO:0007275","GO:0048731","GO:0009628","GO:0032501","GO:0099402","GO:0009791","GO:0009058","GO:0033993","GO:0090567","GO:0009908","GO:0044249","GO:1901576","GO:0050793","GO:0023052","GO:0048827","GO:0009059","GO:0034645","GO:0044271","GO:0007165","GO:0007154","GO:0048608","GO:0061458","GO:0051716","GO:0016070","GO:0051239","GO:2000026","GO:0009888","GO:0048366","GO:0009723","GO:0010016","GO:1901701","GO:0003006","GO:0097305","GO:0009737","GO:0000160","GO:0090698","GO:0009733","GO:0071396","GO:0044260","GO:0048437","GO:0007389","GO:0009653","GO:0010200","GO:0009873","GO:0071369","GO:0071704","GO:0048580","GO:0009416","GO:0048831","GO:0009753","GO:0009314","GO:0090304","GO:0003002","GO:0048507","GO:0071365","GO:0010467","GO:0036294","GO:0071453","GO:0009734","GO:0009739","GO:0008152","GO:0071456","GO:1905392","GO:0090696","GO:0044238","GO:0023051","GO:0048438","GO:1901360","GO:0046483","GO:0010646","GO:0001666","GO:0036293","GO:0048583","GO:0006725","GO:0070482","GO:0010817","GO:0006139","GO:0009611","GO:0009966","GO:0030154","GO:0044237","GO:0009965","GO:0009415","GO:0009605","GO:0009414","GO:1905393","GO:0022414","GO:0000003","GO:0080167","GO:0009738","GO:0048869","GO:0048364","GO:0022622","GO:0014070"],[1667,1837,2497,1298,1007,1007,1032,2566,1142,394,336,664,306,4598,284,7859,46,7889,32,8016,324,236,1133,23,1296,39,27,30,400,313,207,17,1498,1200,8,8,515,208,19,1273,880,138,101,126,12,235,19,19,240,111,8,325,643,29,23,23,23,23,11,105,4,4,34,27,107,86,null,2061,1592,1630,3149,2447,2447,2447,2545,2614,2615,2830,2804,2642,2605,2710,2719,5920,3405,5294,1690,943,3263,3186,998,6674,3012,3008,6430,1037,964,1714,3434,3239,3733,3312,1188,3482,3202,3329,2827,2007,2212,3177,1122,1602,5767,813,560,547,5519,5628,811,1984,553,4146,4065,4039,1942,2283,1263,1266,3279,3919,632,574,705,392,281,198,741,1522,610,603,260,196,414,458,7585,249,212,1017,138,203,210,11761,397,745,262,196,769,4451,174,303,240,4540,241,241,223,135,13241,239,462,195,11090,345,194,5423,5170,349,265,269,772,5313,270,273,4846,221,340,759,11348,111,411,1915,402,125,1926,1935,127,262,841,583,584,323],[682,698,807,506,412,412,416,777,432,133,117,195,102,966,94,1550,26,1552,20,1566,95,73,260,14,289,19,15,16,104,85,61,11,325,266,7,7,127,60,11,276,198,42,33,39,8,62,10,10,63,34,6,81,146,13,11,11,11,11,7,32,4,4,14,12,32,27,1151,741,608,614,968,798,798,798,810,821,821,869,863,824,815,839,841,1486,967,1359,569,378,927,910,392,1617,870,869,1569,398,376,565,950,907,1011,921,429,949,883,910,797,610,655,851,380,486,1341,291,221,216,1277,1295,282,552,213,996,979,969,532,604,378,378,799,925,218,203,232,152,120,95,238,415,205,202,111,92,153,163,1587,105,94,293,70,90,92,2322,143,228,106,86,229,977,78,114,96,989,96,96,91,65,2558,95,152,82,2179,122,81,1149,1101,122,100,101,222,1126,101,101,1035,87,118,217,2214,54,134,459,131,57,457,458,57,94,229,171,171,109],[3.49341724975865e-120,1.67823421235582e-104,2.18308574281622e-80,2.03153636203605e-78,2.16539339765293e-70,2.16539339765293e-70,7.55925494113875e-69,8.11545366519395e-63,5.44122511323004e-62,8.45596025352572e-15,3.02911820739919e-14,5.20591490791628e-14,2.60120228948827e-11,4.80717206533446e-11,2.37684179178145e-10,3.18096489951193e-09,3.70334206206444e-09,6.71010417736292e-09,2.23753662482306e-08,4.36325237409659e-08,1.70493636893304e-07,4.830138137427e-07,2.6496064818211e-06,4.80082717610204e-06,8.57859203202179e-06,9.09042267661985e-06,1.04457669068979e-05,1.06626043165579e-05,1.95557875198908e-05,2.00638856116351e-05,2.16501006655692e-05,2.30169720521419e-05,2.58444261599642e-05,3.08098675362986e-05,3.69557186521355e-05,3.69557186521355e-05,3.83718364858986e-05,5.00968427978002e-05,9.90476514846966e-05,0.000111473574255863,0.000134651553167959,0.000179143121967622,0.000200017440650394,0.000203270191184441,0.000242118052007588,0.000569854299398944,0.000584020859104921,0.000584020859104921,0.000596646365491369,0.000620806400430897,0.000622639080657412,0.000638313285153159,0.000651389109767473,0.000661481208436504,0.000882972477787645,0.000882972477787645,0.000882972477787645,0.000882972477787645,0.000923273167503436,0.000971151812072339,0.000985176859023356,0.000985176859023356,0.00116073621246116,0.0011663254960488,0.00137619404670869,0.0014365184733804,null,2.08030559215139e-98,8.84445753188716e-92,7.99044687862247e-90,1.00136196817579e-85,2.84586426469234e-82,2.84586426469234e-82,2.84586426469234e-82,7.00756126270748e-78,5.95376139055854e-76,7.30906580000395e-76,1.81242752422528e-75,1.99261400683763e-75,1.24136249790299e-74,1.96460409423546e-74,1.97533553385046e-74,2.06279821613432e-74,7.67538383919188e-66,2.89284068409029e-65,3.82543694084247e-65,2.67856771191651e-62,2.80641515897491e-62,3.18419709301977e-62,3.97164317598723e-62,8.23525676222737e-62,1.57900938424425e-61,1.85359397506143e-61,2.03757284488313e-61,2.6534631374365e-61,1.39468209887232e-59,2.60254817546187e-58,2.72038670909883e-58,7.48099772262661e-58,1.7871662950701e-57,3.86452182214192e-57,5.98728191189957e-57,1.15942827321177e-55,1.92111016478512e-54,1.62722245927251e-52,1.72054922933218e-52,6.25109322747183e-51,1.87298563989032e-49,5.98915618143388e-49,1.00035968014342e-44,1.17513090330965e-41,1.04094451961807e-38,1.81418354805277e-37,8.8174668813003e-37,3.37375828508898e-35,1.82553090880975e-34,4.28521139480153e-34,1.76229289056109e-33,3.87278130094404e-33,1.39372135984633e-32,3.53058878516096e-32,1.41310986323543e-31,2.17428141064555e-31,2.36528878189112e-30,1.61612349825563e-29,3.41027282659704e-29,1.10534314184679e-28,1.84797161784541e-28,1.34341062809817e-26,3.52582121396515e-26,2.92183011311569e-25,4.00379761301029e-25,1.15978372803612e-23,1.27607305822323e-23,4.68010170104501e-23,6.21140915965268e-23,1.18377074898836e-22,1.58929187462153e-22,2.98563745757507e-22,9.3415090572645e-22,2.0362142235376e-21,2.04669139904858e-21,2.44206102514965e-21,1.02185290935601e-20,1.60293642271604e-20,7.46594782843918e-20,1.03820276922678e-19,1.17981797952921e-19,5.93169430745261e-19,6.27876335335262e-19,6.80092440869772e-19,7.17365466830399e-19,7.24324615762554e-19,7.53027717301326e-19,1.92985572105862e-18,8.31416865902976e-18,2.25608620757649e-17,4.47840002781269e-17,6.20485402861349e-17,6.61912498339823e-17,1.88537469071654e-16,2.57092935272483e-16,2.59815436809989e-16,2.59815436809989e-16,2.59950335203861e-16,3.05617165431272e-16,3.44602524485475e-16,4.36746530503724e-16,6.67126864497589e-16,9.81946273544576e-16,1.22495210283574e-15,1.32863383466433e-15,2.38797420668098e-15,2.51610912956324e-15,2.91703871897955e-15,3.54126861820981e-15,4.21515603778389e-15,4.51337314889426e-15,4.82521330186237e-15,5.27317350616465e-15,5.97853531428908e-15,1.37142664594468e-14,1.38572967871244e-14,1.45790151577778e-14,1.80024013013447e-14,1.91920602110596e-14,2.0543623294913e-14,5.64837011578974e-14,7.6350405522318e-14,9.19678310958977e-14,1.51253480118034e-13,3.7530095309677e-13,5.67622965819194e-13,8.01680345745693e-13,8.55867558970148e-13,9.09167342406123e-13,9.64628756886268e-13,1.15426640418311e-12,1.35346160577817e-12,1.39828784646005e-12],[2.12399768785326e-117,1.02036640111234e-101,1.32731613163226e-77,1.23517410811792e-75,1.31655918577298e-67,1.31655918577298e-67,4.59602700421236e-66,4.93419582843792e-60,3.30826486884386e-59,5.14122383414364e-12,1.8417038700987e-11,3.1651962640131e-11,1.58153099200887e-08,2.92276061572335e-08,1.44511980940312e-07,1.93402665890325e-06,2.25163197373518e-06,4.07974333983666e-06,1.36042226789242e-05,2.65285744345073e-05,0.000103660131231129,0.000293672398755561,0.00161096074094723,0.00291890292307004,0.00521578395546925,0.00552697698738487,0.00635102627939394,0.00648286342446718,0.0118899188120936,0.0121988424518741,0.0131632612046661,0.0139943190077023,0.0157134111052582,0.0187323994620695,0.0224690769404984,0.0224690769404984,0.0233300765834263,0.0304588804210625,0.0602209721026955,0.0677759331475647,0.0818681443261192,0.108919018156314,0.121610603915439,0.12358827624014,0.147207775620613,0.346471414034558,0.355084682335792,0.355084682335792,0.362760990218752,0.377450291461985,0.378564561039707,0.388094477373121,0.396044578738623,0.402180574729394,0.536847266494888,0.536847266494888,0.536847266494888,0.536847266494888,0.561350085842089,0.590460301739982,0.5989875302862,0.5989875302862,0.705727617176388,0.709125901597669,0.836725980398886,0.873403231815286,null,2.6066229069657e-95,1.10821052874546e-88,1.00120299389139e-86,1.25470654612427e-82,3.5658679236595e-79,3.5658679236595e-79,3.5658679236595e-79,8.78047426217248e-75,7.46006302236985e-73,9.15825944740495e-73,2.27097168785427e-72,2.49674535056755e-72,1.55542720987245e-71,2.46164893007703e-71,2.47509542391462e-71,2.5846861648163e-71,9.61725595050743e-63,3.62472937716513e-62,4.79327248687561e-62,3.35624534303139e-59,3.51643819419556e-59,3.98979895755377e-59,4.97646889951199e-59,1.03187767230709e-58,1.97849875845805e-58,2.32255325075197e-58,2.55307877463856e-58,3.32478931120793e-58,1.74753666988702e-56,3.26099286385372e-55,3.40864454650083e-55,9.37369014645114e-55,2.23931936772284e-54,4.84224584314383e-54,7.50206423561016e-54,1.45276362633435e-52,2.40715103647576e-51,2.03890974146846e-49,2.15584818435322e-49,7.8326198140222e-48,2.34685100678257e-46,7.50441269533666e-46,1.25345067921971e-41,1.47243902184699e-38,1.30430348308145e-35,2.27317198571012e-34,1.10482860022693e-33,4.2273191312165e-32,2.28739022873862e-31,5.36936987768632e-31,2.20815299187305e-30,4.85259497008289e-30,1.74633286388745e-29,4.42382774780669e-29,1.77062665863399e-28,2.72437460753887e-28,2.96370684370957e-27,2.02500274331431e-26,4.27307185172609e-26,1.38499495673403e-25,2.31550843716031e-25,1.68329351700701e-23,4.41785398109834e-23,3.66105313173395e-22,5.01675840910189e-22,1.45320901122925e-20,1.5989195419537e-20,5.8641674314094e-20,7.7828956770448e-20,1.48326474848241e-19,1.99138271890078e-19,3.74100373434156e-19,1.17049108487524e-18,2.55137642209262e-18,2.56450432300788e-18,3.05990246451251e-18,1.28038169542308e-17,2.00847933766319e-17,9.3548326290343e-17,1.30086806984116e-16,1.4783119283501e-16,7.43241296723812e-16,7.86729048175083e-16,8.52155828409824e-16,8.9885892993849e-16,9.0757874355048e-16,9.43543729778562e-16,2.41810921848646e-15,1.04176533297643e-14,2.82687601809335e-14,5.6114352348493e-14,7.7746820978527e-14,8.29376360419798e-14,2.36237448746782e-13,3.22137447896421e-13,3.25548742322916e-13,3.25548742322916e-13,3.25717770010438e-13,3.82938308285384e-13,4.31786963180301e-13,5.47243402721166e-13,8.35909961215479e-13,1.23037868075135e-12,1.53486498485318e-12,1.66477819483441e-12,2.99213168097127e-12,3.15268473934274e-12,3.65504951488138e-12,4.43720957861689e-12,5.28159051534322e-12,5.65525655556451e-12,6.04599226723355e-12,6.60728640322431e-12,7.49110474880422e-12,1.71839758736868e-11,1.73631928742668e-11,1.82675059926956e-11,2.25570088305849e-11,2.40476514444576e-11,2.5741159988526e-11,7.07740775508454e-11,9.56670581194644e-11,1.1523569236316e-10,1.89520610587896e-10,4.70252094230252e-10,7.1123157617145e-10,1.00450547321935e-09,1.0724020513896e-09,1.13918668003487e-09,1.20867983237849e-09,1.44629580444144e-09,1.69588739204005e-09,1.75205467161445e-09],["DNA-binding transcription factor activity","transcription regulator activity","molecular function regulator","sequence-specific DNA binding","transcription regulatory region sequence-specific DNA binding","regulatory region nucleic acid binding","sequence-specific double-stranded DNA binding","DNA binding","double-stranded DNA binding","identical protein binding","DNA-binding transcription factor activity, RNA polymerase II-specific","protein dimerization activity","cis-regulatory region sequence-specific DNA binding","nucleic acid binding","RNA polymerase II cis-regulatory region sequence-specific DNA binding","heterocyclic compound binding","receptor serine/threonine kinase binding","organic cyclic compound binding","hormone binding","protein binding","RNA polymerase II transcription regulatory region sequence-specific DNA binding","molecular transducer activity","protein kinase activity","isoprenoid binding","phosphotransferase activity, alcohol group as acceptor","alcohol binding","protein kinase inhibitor activity","kinase inhibitor activity","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen","secondary active transmembrane transporter activity","signaling receptor activity","peptide receptor activity","kinase activity","transmembrane transporter activity","cyclin-dependent protein serine/threonine kinase inhibitor activity","protein serine/threonine kinase inhibitor activity","active transmembrane transporter activity","protein homodimerization activity","phosphorelay response regulator activity","transporter activity","protein serine/threonine kinase activity","kinase binding","symporter activity","signaling receptor binding","trehalose-phosphatase activity","calmodulin binding","abscisic acid binding","auxin transmembrane transporter activity","anion transmembrane transporter activity","kinase regulator activity","auxin binding","monooxygenase activity","transferase activity, transferring glycosyl groups","oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, incorporation of two atoms of oxygen","3-oxo-arachidoyl-CoA synthase activity","3-oxo-cerotoyl-CoA synthase activity","3-oxo-lignoceronyl-CoA synthase activity","very-long-chain 3-ketoacyl-CoA synthase activity","N6-(Delta2-isopentenyl)-adenosine 5'-monophosphate phosphoribohydrolase activity","protein kinase regulator activity","linoleate 13S-lipoxygenase activity","allene-oxide cyclase activity","carbon-oxygen lyase activity, acting on polysaccharides","pectate lyase activity","inorganic anion transmembrane transporter activity","DNA-binding transcription activator activity",null,"response to organic substance","response to hormone","response to endogenous stimulus","response to chemical","regulation of transcription, DNA-templated","regulation of nucleic acid-templated transcription","regulation of RNA biosynthetic process","regulation of RNA metabolic process","transcription, DNA-templated","nucleic acid-templated transcription","regulation of biosynthetic process","regulation of cellular biosynthetic process","RNA biosynthetic process","regulation of nucleobase-containing compound metabolic process","regulation of cellular macromolecule biosynthetic process","regulation of macromolecule biosynthetic process","regulation of biological process","regulation of cellular metabolic process","regulation of cellular process","response to oxygen-containing compound","hormone-mediated signaling pathway","regulation of primary metabolic process","regulation of nitrogen compound metabolic process","cellular response to hormone stimulus","biological regulation","nucleobase-containing compound biosynthetic process","regulation of gene expression","response to stimulus","cellular response to endogenous stimulus","shoot system development","cellular response to chemical stimulus","organic cyclic compound biosynthetic process","heterocycle biosynthetic process","regulation of metabolic process","aromatic compound biosynthetic process","cellular response to organic substance","regulation of macromolecule metabolic process","anatomical structure development","developmental process","multicellular organism development","system development","response to abiotic stimulus","multicellular organismal process","plant organ development","post-embryonic development","biosynthetic process","response to lipid","reproductive shoot system development","flower development","cellular biosynthetic process","organic substance biosynthetic process","regulation of developmental process","signaling","phyllome development","macromolecule biosynthetic process","cellular macromolecule biosynthetic process","cellular nitrogen compound biosynthetic process","signal transduction","cell communication","reproductive structure development","reproductive system development","cellular response to stimulus","RNA metabolic process","regulation of multicellular organismal process","regulation of multicellular organismal development","tissue development","leaf development","response to ethylene","shoot system morphogenesis","cellular response to oxygen-containing compound","developmental process involved in reproduction","response to alcohol","response to abscisic acid","phosphorelay signal transduction system","post-embryonic plant morphogenesis","response to auxin","cellular response to lipid","cellular macromolecule metabolic process","floral organ development","pattern specification process","anatomical structure morphogenesis","response to chitin","ethylene-activated signaling pathway","cellular response to ethylene stimulus","organic substance metabolic process","regulation of post-embryonic development","response to light stimulus","regulation of shoot system development","response to jasmonic acid","response to radiation","nucleic acid metabolic process","regionalization","meristem development","cellular response to auxin stimulus","gene expression","cellular response to decreased oxygen levels","cellular response to oxygen levels","auxin-activated signaling pathway","response to gibberellin","metabolic process","cellular response to hypoxia","plant organ morphogenesis","post-embryonic plant organ development","primary metabolic process","regulation of signaling","floral whorl development","organic cyclic compound metabolic process","heterocycle metabolic process","regulation of cell communication","response to hypoxia","response to decreased oxygen levels","regulation of response to stimulus","cellular aromatic compound metabolic process","response to oxygen levels","regulation of hormone levels","nucleobase-containing compound metabolic process","response to wounding","regulation of signal transduction","cell differentiation","cellular metabolic process","leaf morphogenesis","response to water","response to external stimulus","response to water deprivation","plant organ formation","reproductive process","reproduction","response to karrikin","abscisic acid-activated signaling pathway","cellular developmental process","root development","root system development","response to organic cyclic compound"],["MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP","BP"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>CLID<\/th>\n      <th>CLSZ<\/th>\n      <th>GOID<\/th>\n      <th>NodeSize<\/th>\n      <th>SampleMatch<\/th>\n      <th>Phyper<\/th>\n      <th>Padj<\/th>\n      <th>Term<\/th>\n      <th>Ont<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scrollX":true,"autoWidth":true,"columnDefs":[{"className":"dt-right","targets":[2,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"orderClasses":false}},"evals":[],"jsHooks":[]}</script>

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
dir_path <- "param/cwl/annotate_peaks"
args <- loadWF(targets = "targets_macs_input.txt", wf_file = "annotate-peaks.cwl", 
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
png("results/seqlogo.png")
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

This project was supported by funds from the National Institutes of
Health (NIH) and the National Science Foundation (NSF).

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-H_Backman2016-bt" class="csl-entry">

H Backman, Tyler W, and Thomas Girke. 2016. “<span class="nocase">systemPipeR: NGS workflow and report generation environment</span>.” *BMC Bioinformatics* 17 (1): 388. <https://doi.org/10.1186/s12859-016-1241-0>.

</div>

<div id="ref-Kaufmann2010-me" class="csl-entry">

Kaufmann, Kerstin, Frank Wellmer, Jose M Muiño, Thilia Ferrier, Samuel E Wuest, Vijaya Kumar, Antonio Serrano-Mislata, et al. 2010. “<span class="nocase">Orchestration of floral initiation by APETALA1</span>.” *Science* 328 (5974): 85–89. <https://doi.org/10.1126/science.1185244>.

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
