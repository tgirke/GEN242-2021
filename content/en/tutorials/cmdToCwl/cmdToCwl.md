---
title: "Automate Creation of CWL Instructions" 
author: "Author: Daniela Cassol, Le Zhang, Thomas Girke"
date: "Last update: 26 May, 2021" 
output:
  html_document:
    toc: true
    toc_float:
        collapsed: true
        smooth_scroll: true
    toc_depth: 3
    fig_caption: yes
    code_folding: show
    number_sections: true
fontsize: 14pt
bibliography: bibtex.bib
weight: 18
type: docs
---

<!---
- Compile from command-line
Rscript -e "rmarkdown::render('cmdToCwl.Rmd', c('html_document'), clean=FALSE); knitr::knit('cmdToCwl.Rmd', tangle=TRUE)"
-->

<div style="text-align: right">

Source code downloads:    
\[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/cmdToCwl/cmdToCwl.Rmd) \]    
\[ [.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/cmdToCwl/cmdToCwl.R) \]

</div>

## Introduction

A central concept for designing workflows within the `systemPipeR` environment is
the use of workflow management containers. `systemPipeR` adopted the widely used community standard [Common Workflow Language](https://www.commonwl.org/) (CWL)
(Amstutz et al. 2016) for describing analysis workflows in a generic and reproducible
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
using CWL’s command line and workflow instructions from within R (R-based approach).
In the latter case the same CWL workflow definition files (e.g. *.cwl* and *.yml*)
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

Let’s explore the `.cwl` file:

``` r
dir_path <- system.file("extdata/cwl/example/", package="systemPipeR")
cwl <- yaml::read_yaml(file.path(dir_path, "example.cwl"))
```

-   The `cwlVersion` component shows the CWL specification version used by the document.
-   The `class` component shows this document describes a command line tool.
    Note that CWL has another `class`, called `Workflow` which represents a union of one
    or more command line tools together.

``` r
cwl[1:2]
```

    ## $cwlVersion
    ## [1] "v1.0"
    ## 
    ## $class
    ## [1] "CommandLineTool"

-   `baseCommand` component provides the name of the software that we desire to execute.

``` r
cwl[3]
```

    ## $baseCommand
    ## [1] "echo"

-   The `inputs` section provides the input information to run the tool. Important
    components of this section are:
    -   `id`: each input has an id describing the input name;
    -   `type`: describe the type of input value (string, int, long, float, double,
        File, Directory or Any);
    -   `inputBinding`: It is optional. This component indicates if the input
        parameter should appear on the command line. If this component is missing
        when describing an input parameter, it will not appear in the command line
        but can be used to build the command line.

``` r
cwl[4]
```

    ## $inputs
    ## $inputs$message
    ## $inputs$message$type
    ## [1] "string"
    ## 
    ## $inputs$message$inputBinding
    ## $inputs$message$inputBinding$position
    ## [1] 1
    ## 
    ## 
    ## 
    ## $inputs$SampleName
    ## $inputs$SampleName$type
    ## [1] "string"
    ## 
    ## 
    ## $inputs$results_path
    ## $inputs$results_path$type
    ## [1] "Directory"

-   The `outputs` section should provide a list of the expected outputs after running the command line tools. Important
    components of this section are:
    -   `id`: each input has an id describing the output name;
    -   `type`: describe the type of output value (string, int, long, float, double,
        File, Directory, Any or `stdout`);
    -   `outputBinding`: This component defines how to set the outputs values. The `glob` component will define the name of the output value.

``` r
cwl[5]
```

    ## $outputs
    ## $outputs$string
    ## $outputs$string$type
    ## [1] "stdout"

-   `stdout`: component to specify a `filename` to capture standard output.
    Note here we are using a syntax that takes advantage of the inputs section,
    using results\_path parameter and also the `SampleName` to construct the output `filename.`

``` r
cwl[6]
```

    ## $stdout
    ## [1] "$(inputs.results_path.basename)/$(inputs.SampleName).txt"

Next, let’s explore the *.yml* file, which provide the parameter values for all
the components we describe above.

For this simple example, we have three parameters defined:

``` r
yaml::read_yaml(file.path(dir_path, "example_single.yml"))
```

    ## $message
    ## [1] "Hello World!"
    ## 
    ## $SampleName
    ## [1] "M1"
    ## 
    ## $results_path
    ## $results_path$class
    ## [1] "Directory"
    ## 
    ## $results_path$path
    ## [1] "./results"

Note that if we define an input component in the *.cwl* file, this value needs
to be also defined here in the *.yml* file.

### How to connect CWL description files within `systemPipeR`

`SYSargs2` container stores all the information and instructions needed for processing
a set of input files with a single or many command-line steps within a workflow
(i.e. several components of the software or several independent software tools).
The `SYSargs2` object is created and fully populated with `loadWF` and \`renderWF\`\` construct
functions.

The following imports a `.cwl` file (here `example.cwl`) for running the `echo Hello World`
example.

``` r
HW <- loadWF( wf_file="example.cwl", input_file="example_single.yml",
              dir_path = dir_path)
HW <- renderWF(HW)
HW
```

    ## Instance of 'SYSargs2':
    ##    Slot names/accessors: 
    ##       targets: 0 (...), targetsheader: 0 (lines)
    ##       modules: 0
    ##       wf: 0, clt: 1, yamlinput: 3 (components)
    ##       input: 1, output: 1
    ##       cmdlist: 1
    ##    WF Steps:
    ##       1. example (rendered: TRUE)

``` r
cmdlist(HW)
```

    ## $defaultid
    ## $defaultid$example
    ## [1] "echo Hello World! > results/M1.txt"

However, we are limited to run just one command line or one sample in this example.
To scale the command line over many samples, a simple solution offered by `systemPipeR`
is to provide a `variable` for each of the parameters that we want to run with multiple samples.

Let’s explore the example:

``` r
yml <- yaml::read_yaml(file.path(dir_path, "example.yml"))
yml
```

    ## $message
    ## [1] "_STRING_"
    ## 
    ## $SampleName
    ## [1] "_SAMPLE_"
    ## 
    ## $results_path
    ## $results_path$class
    ## [1] "Directory"
    ## 
    ## $results_path$path
    ## [1] "./results"

For the `message` and `SampleName` parameter, we are passing a variable connecting
with a third file called `targets.`

Now, let’s explore the `targets` file structure:

``` r
targetspath <- system.file("extdata/cwl/example/targets_example.txt", package="systemPipeR")
read.delim(targetspath, comment.char = "#")
```

    ##                Message SampleName
    ## 1         Hello World!         M1
    ## 2           Hello USA!         M2
    ## 3 Hello Bioconcudctor!         M3

The `targets` file defines all input files or values and sample ids of an analysis workflow.
For this example, we have defined a string message for the `echo` command line tool,
in the first column that will be evaluated,
and the second column is the `SampleName` id for each one of the messages.
Any number of additional columns can be added as needed.

Users should note here, the usage of `targets` files is optional when using
`systemPipeR's` new CWL interface. Since for organizing experimental variables
targets files are extremely useful and user-friendly. Thus, we encourage users to keep using them.

#### How to connect the parameter files and `targets` file information?

The constructor functions create an `SYSargs2` S4 class object connecting three input files:

    - CWL command line specification file (`wf_file` argument);
    - Input variables (`input_file` argument);
    - Targets file (`targets` argument).

As demonstrated above, the latter is optional for workflow steps lacking input files.
The connection between input variables (here defined by `input_file` argument)
and the `targets` file are defined under the `inputvars` argument.
A named vector is required, where each element name needs to match with column
names in the `targets` file, and the value must match the names of the *.yml*
variables. This is used to replace the CWL variable and construct all the command-line
for that particular step.

The variable pattern `_XXXX_` is used to distinguish CWL variables that target
columns will replace. This pattern is recommended for consistency and easy identification
but not enforced.

The following imports a `.cwl` file (same example demonstrated above) for running
the `echo Hello World` example. However, now we are connecting the variable defined
on the `.yml` file with the `targets` file inputs.

``` r
HW_mul <- loadWorkflow(targets = targetspath, wf_file="example.cwl",
                   input_file="example.yml", dir_path = dir_path)
HW_mul <- renderWF(HW_mul, inputvars = c(Message = "_STRING_", SampleName = "_SAMPLE_"))
HW_mul
```

    ## Instance of 'SYSargs2':
    ##    Slot names/accessors: 
    ##       targets: 3 (M1...M3), targetsheader: 1 (lines)
    ##       modules: 0
    ##       wf: 0, clt: 1, yamlinput: 3 (components)
    ##       input: 3, output: 3
    ##       cmdlist: 3
    ##    WF Steps:
    ##       1. example (rendered: TRUE)

``` r
cmdlist(HW_mul)
```

    ## $M1
    ## $M1$example
    ## [1] "echo Hello World! > results/M1.txt"
    ## 
    ## 
    ## $M2
    ## $M2$example
    ## [1] "echo Hello USA! > results/M2.txt"
    ## 
    ## 
    ## $M3
    ## $M3$example
    ## [1] "echo Hello Bioconcudctor! > results/M3.txt"

![**Figure 1:** Connectivity between CWL param files and targets files.](../images/SPR_CWL_hello.png)

## Creating the CWL param files from the command line

Users need to define the command line in a pseudo-bash script format:

``` r
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
```

### Define prefix and defaults

-   First line is the base command. Each line is an argument with its default value.

-   For argument lines (starting from the second line), any word before the first
    space with leading `-` or `--` in each will be treated as a prefix, like `-S` or
    `--min`. Any line without this first word will be treated as no prefix.

-   All defaults are placed inside `<...>`.

-   First argument is the input argument type. `F` for “File,” “int,” “string” are unchanged.

-   Optional: use the keyword `out` followed the type with a `,` comma separation to
    indicate if this argument is also an CWL output.

-   Then, use `:` to separate keywords and default values, any non-space value after the `:`
    will be treated as the default value.

-   If any argument has no default value, just a flag, like `--verbose`, there is no need to add any `<...>`

### `createParamFiles` Function

`createParamFiles` function requires the `string` as defined above as an input.

First of all, the function will print the three components of the `cwl` file:
- `BaseCommand`: Specifies the program to execute.
- `Inputs`: Defines the input parameters of the process.
- `Outputs`: Defines the parameters representing the output of the process.

The four component is the original command line.

If in interactive mode, the function will verify that everything is correct and
will ask you to proceed. Here, the user can answer “no” and provide more
information at the string level. Another question is to save the param created here.

If running the workflow in non-interactive mode, the `createParamFiles` function will
consider “yes” and returning the container.

``` r
cmd <- createParamFiles(command, writeParamFiles = FALSE) 
```

    ## *****BaseCommand*****
    ## hisat2 
    ## *****Inputs*****
    ## S:
    ##     type: File
    ##     preF: -S
    ##     yml: ./results/M1A.sam
    ## x:
    ##     type: File
    ##     preF: -x
    ##     yml: ./data/tair10.fasta
    ## k:
    ##     type: int
    ##     preF: -k
    ##     yml: 1
    ## min-intronlen:
    ##     type: int
    ##     preF: -min-intronlen
    ##     yml: 30
    ## max-intronlen:
    ##     type: int
    ##     preF: -max-intronlen
    ##     yml: 3000
    ## threads:
    ##     type: int
    ##     preF: -threads
    ##     yml: 4
    ## U:
    ##     type: File
    ##     preF: -U
    ##     yml: ./data/SRR446027_1.fastq.gz
    ## *****Outputs*****
    ## output1:
    ##     type: File
    ##     value: ./results/M1A.sam
    ## *****Parsed raw command line*****
    ## hisat2 -S ./results/M1A.sam -x ./data/tair10.fasta -k 1 -min-intronlen 30 -max-intronlen 3000 -threads 4 -U ./data/SRR446027_1.fastq.gz

If the user chooses not to save the `param` files on the above operation,
it can use the `writeParamFiles` function.

``` r
writeParamFiles(cmd, overwrite = TRUE)
```

### How to access and edit param files

#### Print a component

``` r
printParam(cmd, position = "baseCommand") ## Print a baseCommand section
```

    ## *****BaseCommand*****
    ## hisat2

``` r
printParam(cmd, position = "outputs")
```

    ## *****Outputs*****
    ## output1:
    ##     type: File
    ##     value: ./results/M1A.sam

``` r
printParam(cmd, position = "inputs", index = 1:2) ## Print by index
```

    ## *****Inputs*****
    ## S:
    ##     type: File
    ##     preF: -S
    ##     yml: ./results/M1A.sam
    ## x:
    ##     type: File
    ##     preF: -x
    ##     yml: ./data/tair10.fasta

``` r
printParam(cmd, position = "inputs", index = -1:-2) ## Negative indexing printing to exclude certain indices in a position
```

    ## *****Inputs*****
    ## k:
    ##     type: int
    ##     preF: -k
    ##     yml: 1
    ## min-intronlen:
    ##     type: int
    ##     preF: -min-intronlen
    ##     yml: 30
    ## max-intronlen:
    ##     type: int
    ##     preF: -max-intronlen
    ##     yml: 3000
    ## threads:
    ##     type: int
    ##     preF: -threads
    ##     yml: 4
    ## U:
    ##     type: File
    ##     preF: -U
    ##     yml: ./data/SRR446027_1.fastq.gz

#### Subsetting the command line

``` r
cmd2 <- subsetParam(cmd, position = "inputs", index = 1:2, trim = TRUE)
```

    ## *****Inputs*****
    ## S:
    ##     type: File
    ##     preF: -S
    ##     yml: ./results/M1A.sam
    ## x:
    ##     type: File
    ##     preF: -x
    ##     yml: ./data/tair10.fasta
    ## *****Parsed raw command line*****
    ## hisat2 -S ./results/M1A.sam -x ./data/tair10.fasta

``` r
cmdlist(cmd2)
```

    ## $defaultid
    ## $defaultid$hisat2
    ## [1] "hisat2 -S ./results/M1A.sam -x ./data/tair10.fasta"

``` r
cmd2 <- subsetParam(cmd, position = "inputs", index = c("S", "x"), trim = TRUE)
```

    ## *****Inputs*****
    ## S:
    ##     type: File
    ##     preF: -S
    ##     yml: ./results/M1A.sam
    ## x:
    ##     type: File
    ##     preF: -x
    ##     yml: ./data/tair10.fasta
    ## *****Parsed raw command line*****
    ## hisat2 -S ./results/M1A.sam -x ./data/tair10.fasta

``` r
cmdlist(cmd2)
```

    ## $defaultid
    ## $defaultid$hisat2
    ## [1] "hisat2 -S ./results/M1A.sam -x ./data/tair10.fasta"

#### Replacing a existing argument in the command line

``` r
cmd3 <- replaceParam(cmd, "base", index = 1, replace = list(baseCommand = "bwa"))
```

    ## Replacing baseCommand
    ## *****BaseCommand*****
    ## bwa 
    ## *****Parsed raw command line*****
    ## bwa -S ./results/M1A.sam -x ./data/tair10.fasta -k 1 -min-intronlen 30 -max-intronlen 3000 -threads 4 -U ./data/SRR446027_1.fastq.gz

``` r
cmdlist(cmd3)
```

    ## $defaultid
    ## $defaultid$hisat2
    ## [1] "bwa -S ./results/M1A.sam -x ./data/tair10.fasta -k 1 -min-intronlen 30 -max-intronlen 3000 -threads 4 -U ./data/SRR446027_1.fastq.gz"

``` r
new_inputs <- new_inputs <- list(
    "new_input1" = list(type = "File", preF="-b", yml ="myfile"),
    "new_input2" = "-L <int: 4>"
)
cmd4 <- replaceParam(cmd, "inputs", index = 1:2, replace = new_inputs)
```

    ## Replacing inputs
    ## *****Inputs*****
    ## new_input1:
    ##     type: File
    ##     preF: -b
    ##     yml: myfile
    ## new_input2:
    ##     type: int
    ##     preF: -L
    ##     yml: 4
    ## k:
    ##     type: int
    ##     preF: -k
    ##     yml: 1
    ## min-intronlen:
    ##     type: int
    ##     preF: -min-intronlen
    ##     yml: 30
    ## max-intronlen:
    ##     type: int
    ##     preF: -max-intronlen
    ##     yml: 3000
    ## threads:
    ##     type: int
    ##     preF: -threads
    ##     yml: 4
    ## U:
    ##     type: File
    ##     preF: -U
    ##     yml: ./data/SRR446027_1.fastq.gz
    ## *****Parsed raw command line*****
    ## hisat2 -b myfile -L 4 -k 1 -min-intronlen 30 -max-intronlen 3000 -threads 4 -U ./data/SRR446027_1.fastq.gz

``` r
cmdlist(cmd4)
```

    ## $defaultid
    ## $defaultid$hisat2
    ## [1] "hisat2 -b myfile -L 4 -k 1 -min-intronlen 30 -max-intronlen 3000 -threads 4 -U ./data/SRR446027_1.fastq.gz"

#### Adding new arguments

``` r
newIn <- new_inputs <- list(
    "new_input1" = list(type = "File", preF="-b1", yml ="myfile1"),
    "new_input2" = list(type = "File", preF="-b2", yml ="myfile2"),
    "new_input3" = "-b3 <F: myfile3>"
)
cmd5 <- appendParam(cmd, "inputs", index = 1:2, append = new_inputs)
```

    ## Replacing inputs
    ## *****Inputs*****
    ## S:
    ##     type: File
    ##     preF: -S
    ##     yml: ./results/M1A.sam
    ## x:
    ##     type: File
    ##     preF: -x
    ##     yml: ./data/tair10.fasta
    ## k:
    ##     type: int
    ##     preF: -k
    ##     yml: 1
    ## min-intronlen:
    ##     type: int
    ##     preF: -min-intronlen
    ##     yml: 30
    ## max-intronlen:
    ##     type: int
    ##     preF: -max-intronlen
    ##     yml: 3000
    ## threads:
    ##     type: int
    ##     preF: -threads
    ##     yml: 4
    ## U:
    ##     type: File
    ##     preF: -U
    ##     yml: ./data/SRR446027_1.fastq.gz
    ## new_input1:
    ##     type: File
    ##     preF: -b1
    ##     yml: myfile1
    ## new_input2:
    ##     type: File
    ##     preF: -b2
    ##     yml: myfile2
    ## new_input3:
    ##     type: File
    ##     preF: -b3
    ##     yml: myfile3
    ## *****Parsed raw command line*****
    ## hisat2 -S ./results/M1A.sam -x ./data/tair10.fasta -k 1 -min-intronlen 30 -max-intronlen 3000 -threads 4 -U ./data/SRR446027_1.fastq.gz -b1 myfile1 -b2 myfile2 -b3 myfile3

``` r
cmdlist(cmd5)
```

    ## $defaultid
    ## $defaultid$hisat2
    ## [1] "hisat2 -S ./results/M1A.sam -x ./data/tair10.fasta -k 1 -min-intronlen 30 -max-intronlen 3000 -threads 4 -U ./data/SRR446027_1.fastq.gz -b1 myfile1 -b2 myfile2 -b3 myfile3"

``` r
cmd6 <- appendParam(cmd, "inputs", index = 1:2, after=0, append = new_inputs)
```

    ## Replacing inputs
    ## *****Inputs*****
    ## new_input1:
    ##     type: File
    ##     preF: -b1
    ##     yml: myfile1
    ## new_input2:
    ##     type: File
    ##     preF: -b2
    ##     yml: myfile2
    ## new_input3:
    ##     type: File
    ##     preF: -b3
    ##     yml: myfile3
    ## S:
    ##     type: File
    ##     preF: -S
    ##     yml: ./results/M1A.sam
    ## x:
    ##     type: File
    ##     preF: -x
    ##     yml: ./data/tair10.fasta
    ## k:
    ##     type: int
    ##     preF: -k
    ##     yml: 1
    ## min-intronlen:
    ##     type: int
    ##     preF: -min-intronlen
    ##     yml: 30
    ## max-intronlen:
    ##     type: int
    ##     preF: -max-intronlen
    ##     yml: 3000
    ## threads:
    ##     type: int
    ##     preF: -threads
    ##     yml: 4
    ## U:
    ##     type: File
    ##     preF: -U
    ##     yml: ./data/SRR446027_1.fastq.gz
    ## *****Parsed raw command line*****
    ## hisat2 -b1 myfile1 -b2 myfile2 -b3 myfile3 -S ./results/M1A.sam -x ./data/tair10.fasta -k 1 -min-intronlen 30 -max-intronlen 3000 -threads 4 -U ./data/SRR446027_1.fastq.gz

``` r
cmdlist(cmd6)
```

    ## $defaultid
    ## $defaultid$hisat2
    ## [1] "hisat2 -b1 myfile1 -b2 myfile2 -b3 myfile3 -S ./results/M1A.sam -x ./data/tair10.fasta -k 1 -min-intronlen 30 -max-intronlen 3000 -threads 4 -U ./data/SRR446027_1.fastq.gz"

#### Editing `output` param

``` r
new_outs <- list(
    "sam_out" = "<F: $(inputs.results_path)/test.sam>"
) 
cmd7 <- replaceParam(cmd, "outputs", index = 1, replace = new_outs)
```

    ## Replacing outputs
    ## *****Outputs*****
    ## sam_out:
    ##     type: File
    ##     value: $(inputs.results_path)/test.sam
    ## *****Parsed raw command line*****
    ## hisat2 -S ./results/M1A.sam -x ./data/tair10.fasta -k 1 -min-intronlen 30 -max-intronlen 3000 -threads 4 -U ./data/SRR446027_1.fastq.gz

``` r
output(cmd7) 
```

    ## $defaultid
    ## $defaultid$hisat2
    ## [1] "./results/test.sam"

## Version information

``` r
sessionInfo()
```

    ## R version 4.1.0 (2021-05-18)
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
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] systemPipeR_1.27.3          ShortRead_1.50.0           
    ##  [3] GenomicAlignments_1.28.0    SummarizedExperiment_1.22.0
    ##  [5] Biobase_2.52.0              MatrixGenerics_1.4.0       
    ##  [7] matrixStats_0.58.0          BiocParallel_1.26.0        
    ##  [9] Rsamtools_2.8.0             Biostrings_2.60.0          
    ## [11] XVector_0.32.0              GenomicRanges_1.44.0       
    ## [13] GenomeInfoDb_1.28.0         IRanges_2.26.0             
    ## [15] S4Vectors_0.30.0            BiocGenerics_0.38.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] colorspace_2.0-1         rjson_0.2.20             hwriter_1.3.2           
    ##   [4] ellipsis_0.3.2           bit64_4.0.5              AnnotationDbi_1.54.0    
    ##   [7] fansi_0.4.2              splines_4.1.0            cachem_1.0.5            
    ##  [10] knitr_1.33               jsonlite_1.7.2           annotate_1.70.0         
    ##  [13] GO.db_3.13.0             dbplyr_2.1.1             png_0.1-7               
    ##  [16] pheatmap_1.0.12          graph_1.70.0             compiler_4.1.0          
    ##  [19] httr_1.4.2               GOstats_2.58.0           backports_1.2.1         
    ##  [22] assertthat_0.2.1         Matrix_1.3-3             fastmap_1.1.0           
    ##  [25] limma_3.48.0             htmltools_0.5.1.1        prettyunits_1.1.1       
    ##  [28] tools_4.1.0              gtable_0.3.0             glue_1.4.2              
    ##  [31] GenomeInfoDbData_1.2.6   Category_2.58.0          dplyr_1.0.6             
    ##  [34] rsvg_2.1.2               batchtools_0.9.15        rappdirs_0.3.3          
    ##  [37] V8_3.4.2                 Rcpp_1.0.6               jquerylib_0.1.4         
    ##  [40] vctrs_0.3.8              blogdown_1.3             rtracklayer_1.52.0      
    ##  [43] xfun_0.23                stringr_1.4.0            lifecycle_1.0.0         
    ##  [46] restfulr_0.0.13          XML_3.99-0.6             edgeR_3.34.0            
    ##  [49] zlibbioc_1.38.0          scales_1.1.1             BSgenome_1.60.0         
    ##  [52] VariantAnnotation_1.38.0 hms_1.1.0                RBGL_1.68.0             
    ##  [55] RColorBrewer_1.1-2       yaml_2.2.1               curl_4.3.1              
    ##  [58] memoise_2.0.0            ggplot2_3.3.3            sass_0.4.0              
    ##  [61] biomaRt_2.48.0           latticeExtra_0.6-29      stringi_1.6.2           
    ##  [64] RSQLite_2.2.7            genefilter_1.74.0        BiocIO_1.2.0            
    ##  [67] checkmate_2.0.0          GenomicFeatures_1.44.0   filelock_1.0.2          
    ##  [70] DOT_0.1                  rlang_0.4.11             pkgconfig_2.0.3         
    ##  [73] bitops_1.0-7             evaluate_0.14            lattice_0.20-44         
    ##  [76] purrr_0.3.4              bit_4.0.4                tidyselect_1.1.1        
    ##  [79] GSEABase_1.54.0          AnnotationForge_1.34.0   magrittr_2.0.1          
    ##  [82] bookdown_0.22            R6_2.5.0                 generics_0.1.0          
    ##  [85] base64url_1.4            DelayedArray_0.18.0      DBI_1.1.1               
    ##  [88] withr_2.4.2              pillar_1.6.1             survival_3.2-11         
    ##  [91] KEGGREST_1.32.0          RCurl_1.98-1.3           tibble_3.1.2            
    ##  [94] crayon_1.4.1             utf8_1.2.1               BiocFileCache_2.0.0     
    ##  [97] rmarkdown_2.8            jpeg_0.1-8.1             progress_1.2.2          
    ## [100] locfit_1.5-9.4           grid_4.1.0               data.table_1.14.0       
    ## [103] blob_1.2.1               Rgraphviz_2.36.0         digest_0.6.27           
    ## [106] xtable_1.8-4             brew_1.0-6               munsell_0.5.0           
    ## [109] bslib_0.2.5.1

## Funding

This project is funded by NSF award [ABI-1661152](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1661152).

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Amstutz2016-ka" class="csl-entry">

Amstutz, Peter, Michael R Crusoe, Nebojša Tijanić, Brad Chapman, John Chilton, Michael Heuer, Andrey Kartashov, et al. 2016. “Common Workflow Language, V1.0,” July. <https://doi.org/10.6084/m9.figshare.3115156.v2>.

</div>

</div>
