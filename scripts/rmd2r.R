# Convert a .Rmd file to a .R script
library(knitr)
setwd("~/projects/thesis/noisy-microbes/scripts")
purl("dada2_pipeline.Rmd")
q()