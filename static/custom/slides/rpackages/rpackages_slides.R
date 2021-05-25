## pre {

##   max-height: 300px;

##   overflow-y: auto;

## }

## 
## pre[class] {

##   max-height: 300px;

## }


## .scroll-300 {

##   max-height: 300px;

##   overflow-y: auto;

##   background-color: inherit;

## }


## ----package_skeleton1, eval=FALSE--------------------------------------------
## ## Download R script (here pkg_build_fct.R) containing two sample functions
## download.file("https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/manuals/rprogramming/helper_functions/pkg_build_fct.R", "pkg_build_fct.R")
## ## Build package skeleton based on functions in pkg_build_fct.R
## package.skeleton(name="mypackage", code_files=c("pkg_build_fct.R"))


## ----r_build_package, eval=FALSE----------------------------------------------
## system("R CMD build mypackage")
## system("R CMD check mypackage_1.0.tar.gz")


## ----install_package, eval=FALSE----------------------------------------------
## install.packages("mypackage_1.0.tar.gz", repos=NULL)


## ----package_devtools_skeleton1, eval=FALSE-----------------------------------
## library("devtools"); library("roxygen2"); library("usethis"); library(sinew) # If not availble install these packages with 'install.packages(...)'
## create("myfirstpkg") # Creates package skeleton. The chosen name (here myfirstpkg) will be the name of the package.
## setwd("myfirstpkg") # Set working directory of R session to package directory 'myfirstpkg'
## use_mit_license() # Add license information to description file (here MIT). To look up alternatives, do ?use_mit_license


## ----package_devtools_skeleton2, eval=FALSE-----------------------------------
## download.file("https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/manuals/rprogramming/helper_functions/pkg_build_fct.R", "R/pkg_build_fct.R")


## ----package_devtools_skeleton3, eval=FALSE-----------------------------------
## load_all() # Loads package in a simulated way without installing it.
## writeLines(makeOxygen(myMAcomp), "myroxylines") # This creates a 'myroxylines' file in current directory. Delete this file after adding its content to the corresponding functions.


## ----package_devtools_skeleton4, eval=FALSE-----------------------------------
## document() # Auto-generates/updates *.Rd files under man directory (here: myMAcomp.Rd and talkToMe.Rd)
## tools::Rd2txt("man/myMAcomp.Rd") # Renders Rd file from source
## tools::checkRd("man/myMAcomp.Rd") # Checks Rd file for problems


## ----package_devtools_skeleton5, eval=FALSE-----------------------------------
## use_vignette("introduction", title="Introduction to this package")


## ----package_devtools_skeleton6, eval=FALSE-----------------------------------
## setwd("..") # Redirect R session to parent directory
## check("myfirstpkg") # Check package for problems, when in pkg dir one can just use check()
## # remove.packages("myfirstpkg") # Optional. Removes test package if already installed
## install("myfirstpkg", build_vignettes=TRUE) # Installs package
## build("myfirstpkg") # Creates *.tar.gz file for package required to for submission to CRAN/Bioc


## ----package_devtools_skeleton7, eval=FALSE-----------------------------------
## library("myfirstpkg")
## library(help="myfirstpkg")
## ?myMAcomp
## vignette("introduction", "myfirstpkg")


## ----package_devtools_skeleton8, eval=FALSE-----------------------------------
## devtools::install_github("<github_user_name>/<mypkg_repos>", subdir="myfirstpkg") # If the package is in the root directory of the repos, then the 'subdir' argument can be dropped.


## ----sessionInfo--------------------------------------------------------------
sessionInfo()

