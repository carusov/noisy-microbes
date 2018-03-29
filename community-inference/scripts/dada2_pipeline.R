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
  # make_option(c("-f", "--ftrunc"), type = "integer", default = 230,
  #             help = "forward read truncate position"),
  # make_option(c("-b", "--rtrunc"), type = "integer", default = 210,
  #             help = "reverse read truncate position"),
  make_option(c("-s", "--min_len"), type = "integer", default = 221, 
              help = "minimum merged length"),
  make_option(c("-l", "--max_len"), type = "integer", default = 225,
              help = "maximum merged length"),
  make_option(c("-F", "--maxee_F"), type = "numeric", default = 2.5,
              help = "max EE forward reads"),
  make_option(c("-R", "--maxee_R"), type = "numeric", default = 2.5,
              help = "max EE reverse reads")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cat("DADA2 parameters:\n")
cat("Input directory:", opt$input, "\n")
cat("Output directory:", opt$output, "\n")
cat("Minimum sequence length:", opt$min_len, "\n")
cat("Maximum sequence length:", opt$max_len, "\n")
cat("Max EE forward reads:", opt$maxee_F, "\n")
cat("Max EE reverse reads:", opt$maxee_R, "\n")


## ----paths---------------------------------------------------------------

data_path <- opt$input     # parent directory for raw and filtered data
dada2_path <- opt$output    # directory for outputs of DADA2 read processsing
# data_path <- "~/thesis/data/test"
# dada2_path <- file.path(data_path, "dada2")     # directory where DADA2 processing results will be stored
trunc_path <- file.path(data_path, "truncated")
filt_path <- file.path(dada2_path, "filtered")     # directory where filtered reads will be stored
# ref_path <- "~/thesis/references"    # directory containing reference databases
# raw_path <- file.path(data_path, "raw")     # directory containing raw read files

if (!file_test("-d", dada2_path)) dir.create(dada2_path)
if (!file_test("-d", filt_path)) dir.create(filt_path)


## ----plot qualities------------------------------------------------------

#setwd(raw_path)
file_names <- list.files(trunc_path)
fastqs <- str_subset(file_names, ".fastq$")
trunc_Fs <- str_subset(fastqs, "_R1")
trunc_Rs <- str_subset(fastqs, "_R2")

#get the sample names
sample_names <- sapply(str_split(trunc_Fs, "_trunc_R\\d"), `[`, 1)


## ----filter--------------------------------------------------------------

# Define file names for the filtered reads
filt_Fs <- paste0(sample_names, "_filt_R1.fastq")
filt_Rs <- paste0(sample_names, "_filt_R2.fastq")

# Filter paired read sets
filt_stats <- filterAndTrim(fwd = file.path(trunc_path, trunc_Fs), filt = file.path(filt_path, filt_Fs),
                            rev = file.path(trunc_path, trunc_Rs), filt.rev = file.path(filt_path, filt_Rs),
                            #truncLen = c(240, 220), trimLeft = 15, maxEE = c(3, 4), truncQ = 2, rm.phix = TRUE, 
                            maxEE = c(opt$maxee_F, opt$maxee_R), 
                            #truncLen = c(opt$ftrunc, opt$rtrunc), trimLeft = 15, truncQ = 2, rm.phix = TRUE, 
                            compress = FALSE, verbose = TRUE, multithread = TRUE)


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

print("Sequence lengths before length filtering:")
table(nchar(getSequences(sv_table)))


## ----remove bad lengths--------------------------------------------------

min_len <- opt$min_len
max_len <- opt$max_len
sv_table <- sv_table[, nchar(getSequences(sv_table)) %in% seq(min_len, max_len)]

print("Sequence lengths after length filtering:")
table(nchar(getSequences(sv_table)))


## ----remove chimeras-----------------------------------------------------

sv_table.no_chim <- removeBimeraDenovo(sv_table, method = "consensus", verbose = TRUE)

#check what percentage of reads remain
print("Percentage of reads remaining after bimera removal:")
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
print("Read counts at each stage of the DADA2 pipeline:")
track_table

save(mergers, sv_table, sv_table.no_chim, track_table, file = file.path(dada2_path, "tables.RData"))
write.table(sv_table.no_chim, file = file.path(dada2_path, "sv_table.no_chim.txt"), quote = FALSE, sep = "\t")


