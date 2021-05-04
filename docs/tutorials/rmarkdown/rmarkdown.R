## ----install_rmarkdown, eval=FALSE--------------------------------------------
## install.packages("rmarkdown")


## ----render_rmarkdown, eval=FALSE, message=FALSE------------------------------
## rmarkdown::render("sample.Rmd", clean=TRUE, output_format="html_document")


## $ Rscript -e "rmarkdown::render('sample.Rmd', output_format='html_document', clean=TRUE)"


## $ make -B


## ----kable--------------------------------------------------------------------
library(knitr)
kable(iris[1:12,])


## ----dt-----------------------------------------------------------------------
library(DT)
datatable(iris, filter = 'top', options = list(
  pageLength = 100, scrollX = TRUE, scrollY = "600px", autoWidth = TRUE
))


## ----some_jitter_plot, eval=TRUE----------------------------------------------
library(ggplot2)
dsmall <- diamonds[sample(nrow(diamonds), 1000), ]
ggplot(dsmall, aes(color, price/carat)) + geom_jitter(alpha = I(1 / 2), aes(color=color))


## ----some_custom_inserted_plot, eval=TRUE, warning=FALSE, message=FALSE-------
png("myplot.png")
ggplot(dsmall, aes(color, price/carat)) + geom_jitter(alpha = I(1 / 2), aes(color=color))
dev.off()


## ----rmarkdown_symbolic_link, eval=FALSE--------------------------------------
## cd ~/.html
## ln -s ~/bigdata/today/rmarkdown/sample.html sample.html


## ----sessionInfo--------------------------------------------------------------
sessionInfo()

