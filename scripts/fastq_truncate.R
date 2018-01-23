## ----setup, include=FALSE------------------------------------------------

library("knitr")
knitr::opts_chunk$set(echo = TRUE)
#opts_knit$set(root.dir = "~/projects/thesis/data")

## ------------------------------------------------------------------------

library(optparse)
 
option_list = list(
  make_option(c("-d", "--data"), type = "character", default = NULL, 
              help = "data directory"),
  make_option(c("-f", "--ftrunc"), type = "integer", default = 230,
              help = "forward read truncate position"),
  make_option(c("-b", "--rtrunc"), type = "integer", default = 210,
              help = "reverse read truncate position")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


## ------------------------------------------------------------------------
#source("https://bioconductor.org/biocLite.R")
#biocLite("dada2")
library("dada2"); packageVersion("dada2")
library("ggplot2")
library("stringr")


## ----rename files--------------------------------------------------------

#data_path <- "~/thesis/data/zymo_neat"
data_path <- opt$data
raw_path <- file.path(data_path, "raw")

file_names <- list.files(raw_path)
fastqs <- str_subset(file_names, ".fastq$")

# rename files in MockCommunities folder

if (str_detect(fastqs[1], "^lane1-")){
  for (i in seq_along(fastqs)){
    new_name <- fastqs[i]
    new_name <- str_replace(new_name, "lane1-", "")
    new_name <- str_replace(new_name, "index-[ACGT]+-", "")
    new_name <- str_replace(new_name, "_S\\d{3}_L001", "")
    new_name <- str_replace(new_name, "_001", "")
    #new_name <- str_replace(new_name, "-", ".")
    file.rename(file.path(raw_path, fastqs[i]), file.path(raw_path, new_name))
  }
}


## ------------------------------------------------------------------------

file_names <- list.files(raw_path)
fastqs <- str_subset(file_names, ".fastq$")
fastq_Fs <- str_subset(fastqs, "_R1")
fastq_Rs <- str_subset(fastqs, "_R2")

#get the sample names
sample_names <- sapply(str_split(fastq_Fs, "_R\\d"), `[`, 1)

qual_path <- file.path(raw_path, "quality")
if (!file_test("-d", qual_path)) dir.create(qual_path)


#plot quality
if (length(list.files(qual_path)) == 0){
  for (fq in fastqs){
    plotQualityProfile(file.path(raw_path, fq)) +
    scale_y_continuous(limits = c(10, 40), breaks = seq(10, 40, 5)) +
    scale_x_continuous(limits = c(0, 250), breaks = seq(0, 250, 10)) +
    theme(panel.grid.major = element_line(colour="grey", size=0.5)) +
    ggtitle(str_replace(fq, ".fastq", ""))
  
    fname <- str_replace(fq, "fastq", "png")
    ggsave(file.path(qual_path, fname), width = 10, height = 7)
  }
}


## ------------------------------------------------------------------------

# Set up a directory and file names for the truncated reads
trunc_path <- file.path(data_path, "truncated")
if (!file_test("-d", trunc_path)) dir.create(trunc_path)

trunc_Fs <- paste0(sample_names, "_trunc_R1.fastq")
trunc_Rs <- paste0(sample_names, "_trunc_R2.fastq")

# Truncate and filter paired read sets
for (i in seq_along(fastq_Fs)){
  fastqPairedFilter(file.path(raw_path, c(fastq_Fs[i], fastq_Rs[i])), file.path(trunc_path, c(trunc_Fs[i], trunc_Rs[i])),
                    truncLen = c(opt$ftrunc, opt$rtrunc), trimLeft = c(15, 15),
                    #truncLen = c(240, 220), trimLeft = c(15, 15),
                    maxN = Inf, maxEE = c(Inf, Inf), truncQ = 2, rm.phix = TRUE,
                    compress = FALSE, verbose = TRUE)
}


