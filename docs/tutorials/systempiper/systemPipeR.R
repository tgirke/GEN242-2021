#######################################################################
## systemPipeR: Workflow design and reporting generation environment ##
#######################################################################


## Loading package and documentation

library("systemPipeR") # Loads the package
library(help="systemPipeR") # Lists package info
vignette("systemPipeR") # Opens vignette

## Structure of targets file for single-end (SE) samples

library(systemPipeR)
targetspath <- system.file("extdata", "targets.txt", package="systemPipeR") 
read.delim(targetspath, comment.char = "#")[1:4,]

## Structure of targets file for paired-end (PE) samples

targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")[1:2,1:6]

## Sample comparisons

## Sample comparisons are defined in the header lines of the _`targets`_ file 
## starting with '# <CMP>'. 

readLines(targetspath)[1:4]

The function readComp imports the comparison information and stores it in a 
list. 
 
readComp(file=targetspath, format="vector", delim="-")

## Structure and initialization of SYSargs2

## SYSargs2 stores all the information and instructions needed for processing
## a set of input files with a single or many command-line steps within a workflow
## (e.g. several components of the software or several independent software tools).
## The SYSargs2 object is created and fully populated with the loadWF 
## and renderWF functions, respectively. 

## In CWL, files with the extension *.cwl define the parameters of a chosen 
## command-line step or workflow, while files with the extension *.yml define 
## the input variables of command-line steps. Note, input variables provided
## by a targets file can be passed on to a SYSargs2 instance via the inputvars
## argument of the renderWF function. 

hisat2.cwl <- system.file("extdata", "cwl/hisat2/hisat2-se/hisat2-mapping-se.cwl", package="systemPipeR")
yaml::read_yaml(hisat2.cwl)
hisat2.yml <- system.file("extdata", "cwl/hisat2/hisat2-se/hisat2-mapping-se.yml", package="systemPipeR")
yaml::read_yaml(hisat2.yml)

## The following imports a *.cwl file (here hisat2-mapping-se.cwl) for running
## the short read aligner HISAT2. The loadWF and renderWF 
## functions render the proper command-line strings for each sample and software tool.

library(systemPipeR)
targets <- system.file("extdata", "targets.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-se", package="systemPipeR")
WF <- loadWF(targets=targets, wf_file="hisat2-mapping-se.cwl",
                   input_file="hisat2-mapping-se.yml",
                   dir_path=dir_path)

WF <- renderWF(WF, inputvars=c(FileName="_FASTQ_PATH1_", SampleName="_SampleName_"))

## Several accessor methods are available that are named after the slot names of the _`SYSargs2`_ object. 

names(WF)

## Of particular interest is the cmdlist() method. It constructs the system
## commands for running command-line software as specified by a given *.cwl
## file combined with the paths to the input samples (e.g. FASTQ files) provided
## by a targets file. The example below shows the cmdlist() output for
## running HISAT2 on the first SE read sample. Evaluating the output of
## cmdlist() can be very helpful for designing and debugging *.cwl files
## of new command-line software or changing the parameter settings of existing
## ones.  

cmdlist(WF)[1]

## The output components of SYSargs2 define the expected output files for 
## each step in the workflow; some of which are the input for the next workflow step, 
## here next SYSargs2 instance.

output(WF)[1]
modules(WF)
targets(WF)[1]
targets.as.df(targets(WF))[1:4,1:4]
output(WF)[1]
cwlfiles(WF)
inputvars(WF)

## How to run the workflow on a cluster

## This section of the tutorial provides an introduction to the usage of the
## systemPipeR features on a cluster. 

## Now open the R markdown script `*.Rmd`in your R IDE (_e.g._vim-r or RStudio)
## and run the workflow as outlined below. If you work under Vim-R-Tmux, the
## following command sequence will connect the user in an interactive session with
## a node on the cluster. The code of the `Rmd`
## script can then be sent from Vim on the login (head) node to an open R session running
## on the corresponding computer node. This is important since Tmux sessions
## should not be run on the computer nodes. 

q("no") # closes R session on head node
srun --x11 --partition=short --mem=2gb --cpus-per-task 4 --ntasks 1 --time 2:00:00 --pty bash -l
module load R/4.0.3
R

## Now check whether your R session is running on a computer node of the
## cluster and not on a head node.

system("hostname") # should return name of a compute node starting with i or c 
getwd() # checks current working directory of R session
dir() # returns content of current working directory

## Parallelization on clusters

## Alternatively, the computation can be greatly accelerated by processing many files 
## in parallel using several compute nodes of a cluster, where a scheduling/queuing
## system is used for load balancing. For this the clusterRun function submits 
## the computing requests to the scheduler using the run specifications
## defined by runCommandline. 

## The following example will run 18 processes in parallel using for each 4 CPU
## cores. If the resources available on a cluster allow running all 18 processes
## at the same time then the shown sample submission will utilize in total 72 CPU
## cores. Note, clusterRun can be used with most queueing systems as it is based
## on utilities from the batchtools package which supports the use of template
## files (*.tmpl) for defining the run parameters of different schedulers. The
## following example uses the sample conf and template files for the Slurm
## scheduler provided by this package.  

library(batchtools)
targetspath <- system.file("extdata", "targetsPE.txt", package="systemPipeR")
dir_path <- system.file("extdata/cwl/hisat2/hisat2-pe", package="systemPipeR")
args <- loadWorkflow(targets=targetspath, wf_file="hisat2-mapping-pe.cwl", 
                     input_file="hisat2-mapping-pe.yml", dir_path=dir_path)
args <- renderWF(args, inputvars=c(FileName1="_FASTQ_PATH1_", FileName2="_FASTQ_PATH2_", SampleName="_SampleName_"))
resources <- list(walltime=120, ntasks=1, ncpus=4, memory=1024)
reg <- clusterRun(args, FUN = runCommandline, more.args = list(args=args, make_bam=TRUE, dir=FALSE), 
                  conffile = ".batchtools.conf.R", template = "batchtools.slurm.tmpl", 
                  Njobs=18, runid="01", resourceList=resources)
getStatus(reg=reg)
waitForJobs(reg=reg)

## Workflow templates

## Sample workflow templates are provided via systemPipeRdata and GitHub. Instances of these 
## workflows can be created with a single command. 

### RNA-Seq sample

## Load the RNA-Seq sample workflow into your current working directory.

library(systemPipeRdata)
genWorkenvir(workflow="rnaseq")
setwd("rnaseq")

## Run RNA-Seq workflow interactively

## Open under workflow directory (here 'rnaseq') systemPipeRNAseq.Rmd and continue there
