## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----libraries-----------------------------------------------------------

library("dada2")
library("stringr")
library("ggplot2")


## ------------------------------------------------------------------------

library(optparse)
 
option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "data directory", metavar = NULL),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "results directory", metavar = NULL),
  make_option(c("-f", "--ftrunc"), type = "integer", default = 230,
              help = "forward read truncate position"),
  make_option(c("-b", "--rtrunc"), type = "integer", default = 210,
              help = "reverse read truncate position"),
  make_option(c("-s", "--min_len"), type = "integer", default = 230, 
              help = "minimum merged length"),
  make_option(c("-l", "--max_len"), type = "integer", default = 235,
              help = "maximum merged length")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


## ----paths---------------------------------------------------------------

#data_path <- "~/projects/thesis/data/test_bal"
data_path <- opt$input     # parent directory for raw and filtered data
#result_path <- "~/projects/thesis/results/test_bal"
dada2_path <- opt$output    # directory for outputs of DADA2 read processsing
#ref_path <- "~/projects/thesis/references"    # directory containing reference databases
raw_path <- file.path(data_path, "raw")     # directory containing raw read files
#dada2_path <- file.path(result_path, "dada2")     # directory where DADA2 processing results will be stored
filt_path <- file.path(dada2_path, "filtered")     # directory where filtered reads will be stored

if (!file_test("-d", dada2_path)) dir.create(dada2_path)
if (!file_test("-d", filt_path)) dir.create(filt_path)


## ----plot qualities------------------------------------------------------

#setwd(raw_path)
file_names <- list.files(raw_path)
fastqs <- str_subset(file_names, ".fastq$")
fastq_Fs <- str_subset(fastqs, "_R1")
fastq_Rs <- str_subset(fastqs, "_R2")

#get the sample names
sample_names <- sapply(str_split(fastq_Fs, "_R\\d"), `[`, 1)


## ----filter--------------------------------------------------------------

# Define file names for the filtered reads
filt_Fs <- paste0(sample_names, "_R1_filt.fastq")
filt_Rs <- paste0(sample_names, "_R2_filt.fastq")

# Filter paired read sets
filt_stats <- filterAndTrim(fwd = file.path(raw_path, fastq_Fs), filt = file.path(filt_path, filt_Fs),
                            rev = file.path(raw_path, fastq_Rs), filt.rev = file.path(filt_path, filt_Rs),
                            #truncLen = c(240, 220), trimLeft = 15, maxEE = c(3, 4), truncQ = 2, rm.phix = TRUE, 
                            truncLen = c(opt$ftrunc, opt$rtrunc), trimLeft = 15, maxEE = c(3, 4), truncQ = 2, rm.phix = TRUE, 
                            compress = TRUE, verbose = TRUE, multithread = TRUE)


## ----errors--------------------------------------------------------------

err_F <- learnErrors(file.path(filt_path, filt_Fs), multithread = TRUE)
err_R <- learnErrors(file.path(filt_path, filt_Rs), multithread = TRUE)


## ----dereplicate---------------------------------------------------------

derep_Fs <- derepFastq(file.path(filt_path, filt_Fs), verbose = TRUE)
derep_Rs <- derepFastq(file.path(filt_path, filt_Rs), verbose = TRUE)

#names(derep_Fs) <- sample_names
#names(derep_Rs) <- sample_names


## ----SV inference--------------------------------------------------------

dada_Fs <- dada(derep_Fs, err = err_F, multithread = TRUE)
dada_Rs <- dada(derep_Rs, err = err_R, multithread = TRUE)

# Save the dada objects
save(err_F, err_R, derep_Fs, derep_Rs, dada_Fs, dada_Rs, file = file.path(dada2_path, "dada2.RData"))


## ----merge SVs-----------------------------------------------------------

#load(file = file.path(dada2_path, "dada2.RData"))
mergers <- mergePairs(dada_Fs, derep_Fs, dada_Rs, derep_Rs, 
                     propagateCol = c("n0", "n1", "birth_fold", "birth_ham"), 
                     verbose = TRUE)


## ----sequence table------------------------------------------------------

sv_table <- makeSequenceTable(mergers)
row.names(sv_table) <- sample_names
table(nchar(getSequences(sv_table)))


## ----remove bad lengths--------------------------------------------------

min_len <- opt$min_len
max_len <- opt$max_len
sv_table <- sv_table[, nchar(getSequences(sv_table)) %in% seq(min_len, max_len)]

table(nchar(getSequences(sv_table)))


## ----remove chimeras-----------------------------------------------------

sv_table.no_chim <- removeBimeraDenovo(sv_table, method = "consensus", verbose = TRUE)

#check what percentage of reads remain
sum(sv_table.no_chim) / sum(sv_table)


## ----track reads---------------------------------------------------------

getN <- function(x) sum(getUniques(x))

if (length(sample_names) > 1){
  track_table <- cbind(filt_stats, sapply(dada_Fs, getN), sapply(mergers, getN), rowSums(sv_table), rowSums(sv_table.no_chim))
} else {
  track_table <- cbind(filt_stats, getN(dada_Fs), getN(mergers), sum(sv_table), sum(sv_table.no_chim))
}

colnames(track_table) <- c("raw", "filtered", "denoised", "merged", "tabled", "non_chim")
rownames(track_table) <- sample_names
track_table

save(mergers, sv_table, sv_table.no_chim, track_table, file = file.path(dada2_path, "tables.RData"))
write.table(sv_table.no_chim, file = file.path(dada2_path, "sv_table.no_chim.txt"), quote = FALSE, sep = "\t")


