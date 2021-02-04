# setup
options(repos=structure(c(CRAN="http://cran.rstudio.com/")))
FORCE_UPDATE <- FALSE
if(!requireNamespace("remotes") | FORCE_UPDATE) install.packages("remotes")
if(!requireNamespace("yaml") | FORCE_UPDATE) install.packages("yaml")
pkgs <- yaml::yaml.load_file("deps.yaml")

# define func
install_deps <- function(pkgs, type, FORCE_UPDATE){
  if(length(pkgs) < 1) return(cat("No ", type, " pkgs to install\n"))
  switch(type,
    'CRAN' = {
      for(pkg in pkgs){
        if(!requireNamespace(pkg) | FORCE_UPDATE) install.packages(pkg)
        else cat(pkg, " found skip\n")
      }
    },
    "Bioc" = {
      if(!requireNamespace("BiocManager") | FORCE_UPDATE) install.packages("BiocManager")
      for(pkg in pkgs){
        if(!requireNamespace(pkg) | FORCE_UPDATE) BiocManager::install(pkg, update = FORCE_UPDATE)
        else cat(pkg, " found skip\n")
      }
    },
    "Github" = {
      for(pkg in pkgs){
        if(sum(gregexpr("/", pkg, fixed = TRUE)[[1]] > 0) != 1) stop("Invalid Github pkg name: ", pkg, ". e.g. user/repo")
        if(!requireNamespace(strsplit(pkg, "/")[[1]][2]) | FORCE_UPDATE) remotes::install_github(pkg)
        else cat(pkg, " found skip\n")
      }
    }
  )
  cat("Install ", type, " pkgs done.\n")
}

# install
install_deps(pkgs$CRAN, "CRAN", FORCE_UPDATE)
install_deps(pkgs$Bioc, "Bioc", FORCE_UPDATE)
install_deps(pkgs$Github, "Github", FORCE_UPDATE)
dir.create("public", showWarnings = FALSE, recursive = TRUE)
cat("Deps installation done.\n")
