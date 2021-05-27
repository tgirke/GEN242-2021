## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----load_library, eval=TRUE, include=FALSE-----------------------------------
library(systemPipeR)


## -----------------------------------------------------------------------------
dir_path <- system.file("extdata/cwl/example/", package="systemPipeR")
cwl <- yaml::read_yaml(file.path(dir_path, "example.cwl"))


## -----------------------------------------------------------------------------
cwl[1:2]


## -----------------------------------------------------------------------------
cwl[3]


## -----------------------------------------------------------------------------
cwl[4]


## -----------------------------------------------------------------------------
cwl[5]


## -----------------------------------------------------------------------------
cwl[6]


## -----------------------------------------------------------------------------
yaml::read_yaml(file.path(dir_path, "example_single.yml"))


## ----fromFile, eval=TRUE------------------------------------------------------
HW <- loadWF( wf_file="example.cwl", input_file="example_single.yml",
              dir_path = dir_path)
HW <- renderWF(HW)
HW
cmdlist(HW)


## -----------------------------------------------------------------------------
yml <- yaml::read_yaml(file.path(dir_path, "example.yml"))
yml


## -----------------------------------------------------------------------------
targetspath <- system.file("extdata/cwl/example/targets_example.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")


## ----fromFile_example, eval=TRUE----------------------------------------------
HW_mul <- loadWorkflow(targets = targetspath, wf_file="example.cwl",
                   input_file="example.yml", dir_path = dir_path)
HW_mul <- renderWF(HW_mul, inputvars = c(Message = "_STRING_", SampleName = "_SAMPLE_"))
HW_mul
cmdlist(HW_mul)


## ----cmd, eval=TRUE-----------------------------------------------------------
command <- "
hisat2 \
    -S <F, out: ./results/M1A.sam> \
    -x <F: ./data/tair10.fasta> \
    -k <int: 1> \
    -min-intronlen <int: 30> \
    -max-intronlen <int: 3000> \
    -threads <int: 4> \
    -U <F: ./data/SRR446027_1.fastq.gz>
"


## -----------------------------------------------------------------------------
cmd <- createParamFiles(command, writeParamFiles = FALSE) 


## ----saving, eval=FALSE-------------------------------------------------------
## writeParamFiles(cmd, overwrite = TRUE)


## -----------------------------------------------------------------------------
printParam(cmd, position = "baseCommand") ## Print a baseCommand section
printParam(cmd, position = "outputs")
printParam(cmd, position = "inputs", index = 1:2) ## Print by index
printParam(cmd, position = "inputs", index = -1:-2) ## Negative indexing printing to exclude certain indices in a position


## -----------------------------------------------------------------------------
cmd2 <- subsetParam(cmd, position = "inputs", index = 1:2, trim = TRUE)
cmdlist(cmd2)

cmd2 <- subsetParam(cmd, position = "inputs", index = c("S", "x"), trim = TRUE)
cmdlist(cmd2)



## -----------------------------------------------------------------------------
cmd3 <- replaceParam(cmd, "base", index = 1, replace = list(baseCommand = "bwa"))
cmdlist(cmd3)


## -----------------------------------------------------------------------------
new_inputs <- new_inputs <- list(
    "new_input1" = list(type = "File", preF="-b", yml ="myfile"),
    "new_input2" = "-L <int: 4>"
)
cmd4 <- replaceParam(cmd, "inputs", index = 1:2, replace = new_inputs)
cmdlist(cmd4)


## -----------------------------------------------------------------------------
newIn <- new_inputs <- list(
    "new_input1" = list(type = "File", preF="-b1", yml ="myfile1"),
    "new_input2" = list(type = "File", preF="-b2", yml ="myfile2"),
    "new_input3" = "-b3 <F: myfile3>"
)
cmd5 <- appendParam(cmd, "inputs", index = 1:2, append = new_inputs)
cmdlist(cmd5)

cmd6 <- appendParam(cmd, "inputs", index = 1:2, after=0, append = new_inputs)
cmdlist(cmd6)


## -----------------------------------------------------------------------------
new_outs <- list(
    "sam_out" = "<F: $(inputs.results_path)/test.sam>"
) 
cmd7 <- replaceParam(cmd, "outputs", index = 1, replace = new_outs)
output(cmd7) 


## ----sessionInfo--------------------------------------------------------------
sessionInfo()

