#################################################
## Helper Functions used by scRNA-Seq Projects ##
#################################################
## Author: Anonymous Genius
## Last update: 08-Jun-21

## (A) Function to Arrange Single Cell Embedding Results in Single Plot
arrangePlots <- function(x, method_list, layout=list(nrow=2, ncol=2), plot_title, ...){
    require(gridExtra)
    for(i in seq_along(method_list[["run"]])) {
            x <- suppressWarnings(method_list[["run"]][[i]](x)) 
    }
    plotList <- lapply(seq_along(method_list[["plot"]]), 
                       function(y) method_list[["plot"]][[y]](x, ...))
    ml <- marrangeGrob(plotList, nrow=layout$nrow, ncol=layout$ncol, top=plot_title)
    return(ml)
}
## Usage:
embedMethodList <- list(run=c(runTSNE, runMDS, runUMAP, runPCA), plot=c(plotTSNE, plotMDS, plotUMAP, plotPCA))
arrangePlots(x=sce, method_list=embedMethodList, layout=list(nrow=2, ncol=2), plot_title="My Test Plot", colour_by="label", text_by="label")
