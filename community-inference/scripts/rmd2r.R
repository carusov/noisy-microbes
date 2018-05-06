#!/usr/bin/Rscript

# Convert a .Rmd file to a .R script
library(knitr)
library(optparse)
library(stringr)

option_list = list(
  make_option(c("-i", "--in_file"), type = "character", default = NULL, 
              help = "dataset file name", metavar = "character"),
  make_option(c("-d", "--out_dir"), type = "character", 
              default = "~/thesis/noisy-microbes/community-inference/scripts")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input <- opt$in_file
fn <- str_split(input, "/")[[1]]
fn <- fn[length(fn)]
bn <- str_split(fn, "\\.")[[1]]
if (length(bn) > 2){
  bn <- paste(bn[-length(bn)], collapse = ".")
} else {
  bn <- bn[1]
}
output <- paste0(bn, ".R")

#setwd(opt$out_dir)
purl(input = input, output = file.path(opt$out_dir, output))
q()
