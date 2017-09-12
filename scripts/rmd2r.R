#!/usr/bin/Rscript

# Convert a .Rmd file to a .R script
library(knitr)
library(optparse)

option_list = list(
  make_option(c("-i", "--in_file"), type = "character", default = NULL, 
              help = "dataset file name", metavar = "character")
  # make_option(c("-d", "--out_dir"), type = "character", 
  #             default = "~projects/thesis/noisy-microbes/scripts")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#setwd(opt$out_dir)
purl(opt$in_file)
q()
