## Author: Vincent Caruso
## Date written: 1/19/18
## Last modified: 1/19/18
## Purpose: This script contains many data processing functions for use in 
##      analyzing the results output by various 16S rRNA sequence inference 
##      methods.

library(tidyverse)
library(stringr)
library(Biostrings)
library(ShortRead)

options(stringsAsFactors = FALSE)

# function to load UCLUST OTU table
load_uclust <- function(otu_file_path, seq_file_path, sample_names = NULL){
  uclust_table <- as.tibble(read.table(file = otu_file_path, header = TRUE, sep = "\t", 
                                       skip = 1, comment.char = ""))
  colnames(uclust_table)[1] <- "id"
  
  if (is.null(sample_names)){
    sample_names <- colnames(uclust_table)[-1]  
    sample_names <- sample_names[order(sample_names)]
    uclust_table <- uclust_table %>% select(id, sample_names)  # put columns in order of sample names
  } else {
    tmp_names <- colnames(uparse_table)[-1][order(colnames(uparse_table)[-1])]  
    uparse_table <- uparse_table %>% select(id, tmp_names)  # put sample columns in order
    colnames(uparse_table)[-1] <- sample_names  # rename columns with provided sample names
  }
  
  # get sequences and add them to the table
  uclust_seqs <- readDNAStringSet(seq_file_path)
  uclust_seqs <- as.character(uclust_seqs)
  names(uclust_seqs) <- sapply(str_split(names(uclust_seqs), " "), `[`, 1)  # remove annotation from names
  uclust_table$sequence <- uclust_seqs[uclust_table$id]  
  uclust_table$id <- str_replace(uclust_table$id, "denovo", "denovo_")  # make ids more readable
  
  return(uclust_table)
}


# function to load UPARSE OTU table
load_uparse <- function(otu_file_path, seq_file_path, sample_names = NULL){
  uparse_table <- as.tibble(read.table(file = otu_file_path, header = TRUE, sep = "\t", 
                                       comment.char = ""))
  colnames(uparse_table)[1] <- "id"
  
  if (is.null(sample_names)){
    sample_names <- colnames(uparse_table)[-1]  
    sample_names <- sample_names[order(sample_names)]
  } else {
    tmp_names <- colnames(uparse_table)[-1][order(colnames(uparse_table)[-1])]  
    uparse_table <- uparse_table %>% select(id, tmp_names)  # put sample columns in order
    colnames(uparse_table)[-1] <- sample_names  # rename columns with provided sample names
  }
  
  # get sequences and add them to the table
  uparse_seqs <- readDNAStringSet(seq_file_path)
  uparse_seqs <- as.character(uparse_seqs)
  uparse_table$sequence <- uparse_seqs[uparse_table$id]  # add sequences to table
  uparse_table$id <- str_replace(uparse_table$id, "OTU", "OTU_")
  
  return(uparse_table)
}


# function to load mothur OTU table
load_mothur <- function(otu_file_path, seq_file_path, sample_names = NULL){
  mothur_table <- read.table(file = otu_file_path, header = TRUE, sep = "\t", 
                                       comment.char = "")
  tmp_names <- mothur_table[[2]][order(mothur_table[2])]
  mothur_table <- data.frame(t(mothur_table[, -c(1:3)]))
  mothur_table <- as.tibble(rownames_to_column(mothur_table, var = "id"))
  colnames(mothur_table)[-1] <- tmp_names
  
  if (is.null(sample_names)){
    sample_names <- tmp_names
    mothur_table <- mothur_table %>% select(id, sample_names)  # put columns in order of sample names
  } else {
    tmp_names <- colnames(mothur_table)[-1][order(colnames(mothur_table)[-1])]  
    mothur_table <- mothur_table %>% select(id, tmp_names)  # put sample columns in order
    colnames(mothur_table)[-1] <- sample_names  # rename columns with provided sample names
  }
  
  # get sequences and add them to the table
  mothur_seqs <- readDNAStringSet(seq_file_path)
  mothur_seqs <- as.character(mothur_seqs)
  mothur_table$sequence <- mothur_seqs[mothur_table$id]  # add sequences to table
  mothur_table$id <- str_replace(mothur_table$id, "Otu", "Otu_")
  
  return(mothur_table)
}

# function to load UNOISE ZOTU table
load_unoise <- function(otu_file_path, seq_file_path, sample_names = NULL){
  unoise_table <- as.tibble(read.table(file = otu_file_path, header = TRUE, sep = "\t", 
                                       comment.char = ""))
  colnames(unoise_table)[1] <- "id"
  
  if (is.null(sample_names)){
    sample_names <- colnames(unoise_table)[-1]  
    sample_names <- sample_names[order(sample_names)]
  } else {
    tmp_names <- colnames(unoise_table)[-1][order(colnames(unoise_table)[-1])]  
    unoise_table <- unoise_table %>% select(id, tmp_names)  # put sample columns in order
    colnames(unoise_table)[-1] <- sample_names
  }
  
  # get sequences and add them to the table
  unoise_seqs <- readDNAStringSet(seq_file_path)
  unoise_seqs <- as.character(unoise_seqs)
  unoise_table$sequence <- unoise_seqs[unoise_table$id]  # add sequences to table
  unoise_table$id <- str_replace(unoise_table$id, "OTU", "ZOTU_")
  
  return(unoise_table)
}


# function to load MED Node table
load_med <- function(otu_file_path, seq_file_path, chimera_file_path, sample_names = NULL){
  med_table <- read.table(file = otu_file_path, header = TRUE, sep = "\t")
  row.names(med_table) <- med_table[, 1]
  med_table <- med_table[, -1]  # remove sample name column
  med_table <- data.frame(t(med_table))   # samples as columns
  med_table <- as.tibble(rownames_to_column(med_table, var = "id"))
  med_table$id <- str_replace(med_table$id, "X", "Node_")
  
  if (is.null(sample_names)){
    sample_names <- colnames(med_table)[-1]  
    sample_names <- sample_names[order(sample_names)]
  } else {
    tmp_names <- colnames(med_table)[-1][order(colnames(med_table)[-1])]  
    med_table <- med_table %>% select(id, tmp_names)  # put sample columns in order
    colnames(med_table)[-1] <- sample_names
  }
  
  # get sequences and add them to the table
  med_seqs <- readDNAStringSet(seq_file_path)
  med_seqs <- as.character(med_seqs)
  med_seqs <- sapply(med_seqs, FUN = str_replace, "-+", "")
  names(med_seqs) <- paste0("Node_", sapply(str_split(names(med_seqs), "\\|"), FUN = `[`, 1))  # remove size annotation
  med_table$sequence <- med_seqs[med_table$id]  # add chimera-free sequences to table
  
  # remove chimeras detected with UCHIME method
  med_seqs.chimeras <- readDNAStringSet(chimera_file_path)
  if (length(med_seqs.chimeras) > 0){
    med_seqs.chimeras <- as.character(med_seqs.chimeras)
    names(med_seqs.chimeras) <- paste0("Node_", sapply(str_split(names(med_seqs.chimeras), ";"), FUN = `[`, 1))
    med_table <- med_table[!med_table$sequence %in% med_seqs.chimeras, ]
  }
  
  return(med_table)
}


# function to load Deblur sOTU table
load_deblur <- function(otu_file_path, sample_names = NULL){
  deblur_table <- as.tibble(read.table(file = otu_file_path, header = TRUE, sep = "\t", 
                                       skip = 1, comment.char = ""))
  colnames(deblur_table)[1] <- "sequence"
  deblur_table$sequence <- toupper(deblur_table$sequence)
  # deblur_table$sequence <- as.character(deblur_table$sequence)
  
  if (is.null(sample_names)){
    sample_names <- colnames(deblur_table)[-1]  
    sample_names <- sample_names[order(sample_names)]
  } else {
    tmp_names <- colnames(deblur_table)[-1][order(colnames(deblur_table)[-1])]  
    deblur_table <- deblur_table %>% select(sequence, tmp_names)  # put sample columns in order
    colnames(deblur_table)[-1] <- sample_names
  }  
  
  deblur_table$id <- paste0("sOTU_", 1:nrow(deblur_table))
  deblur_table <- deblur_table %>% select(id, sample_names, sequence)
  
  return(deblur_table)  
}


# function to load DADA2 ASV table
load_dada2 <- function(otu_file_path, sample_names = NULL){
  dada2_table <- read.table(file = otu_file_path, header = TRUE, sep = "\t")
  
  if (length(dada2_table) == 1){
    colnames(dada2_table) <- ifelse(is.null(sample_names), "sample_1", sample_names)
  } else {
    dada2_table <- data.frame(t(dada2_table))
  }
  
  dada2_table <- as.tibble(rownames_to_column(dada2_table, var = "sequence"))
  
  if (is.null(sample_names)){
    sample_names <- colnames(dada2_table)[-1]  
    sample_names <- sample_names[order(sample_names)]
  } else {
    tmp_names <- colnames(dada2_table)[-1][order(colnames(dada2_table)[-1])]  
    dada2_table <- dada2_table %>% select(sequence, tmp_names)  # put sample columns in order
    colnames(dada2_table)[-1] <- sample_names
  }  
  
  dada2_table <- dada2_table[order(dada2_table[[ sample_names[1] ]], decreasing = TRUE), ]
  dada2_table$id <- paste0("ASV_", 1:nrow(dada2_table))
  dada2_table <- dada2_table %>% select(id, sample_names, sequence)
  
  return(dada2_table)
}


# define a function to remove low abundance sequences from any sample
remove_rare <- function(seq_table, sample_names, min_abund = 2, discard_empty = TRUE){
  seq_table[, sample_names][seq_table[, sample_names] < min_abund]  <- 0  # set rare abundances to zero
  
  if (discard_empty){
    seq_table <- seq_table[rowSums(seq_table[, sample_names]) > 0, ]  # remove seqs with all zeros
  }
  
  return(seq_table)
}


# define a function to remove sequences below a minimum abundance across all samples
remove_rare_across <- function(seq_table, sample_names, min_abund = 2){
  seq_table <- seq_table[rowSums(seq_table[, sample_names]) >= min_abund, ]  # remove sequences with less than min_abund reads across samples
  return(seq_table)
}


# define a function to detect if a string is a substring of any in a vector of strings
is_substring <- function(pattern, string){
  return(any(grepl(pattern, string, fixed = TRUE)))
}


# define a function to combine tables from several methods for a single sample into one
merge_tables <- function(table_list, sample_names, collapse = FALSE, id = NULL){
  all_seqs <- sapply(table_list, function(t) return(t$sequence)) %>% unlist() %>% unique()
  all_table <- matrix(0, nrow = length(table_list), ncol = length(all_seqs))
  row.names(all_table) <- names(table_list)
  colnames(all_table) <- all_seqs
  
  # populate the new table with sequence abundances 
  for (n in names(table_list)){
    cols <- colnames(table_list[[n]])  # get the column names of an input table
    # scol <- cols[cols %in% sample_names]
    scol <- cols[sapply(cols, is_substring, sample_names)]  # find the column name that corresponds to the sample
    all_table[n, table_list[[n]]$sequence] <- table_list[[n]][[scol]]  # populate the appropriate merged table row with sample abundances
  }
  
  if (collapse) all_table <- collapseNoMismatch(all_table)
  
  if (is.null(id)){
    prefix <- "Seq"
  } else if (id == "keep"){
    prefix <- str_split(table_list[[1]]$id[1], "_")[[1]][1]
  } else{
    prefix <- id
  }
  
  all_table <- data.frame(t(all_table))
  all_table <- rownames_to_column(all_table, var = "sequence") %>% as.tibble()
  all_table <- all_table[order(all_table[[ names(table_list)[1] ]], decreasing = TRUE), ]
  all_table <- all_table %>% mutate(id = paste(prefix, 1:nrow(all_table), sep = "_")) %>% 
    select(id, names(table_list), sequence)
  
  return(all_table)
}


# define a function to compute Levenshtein distances between two sequences
levDist <- Vectorize(function(query, ref, ...) {
  mmi <- dada2:::nweval(query, ref, ...)    # this gives matches, mismatches, and indels
  ldist <- mmi[2] + mmi[3]
  if (nchar(query) > sum(mmi)) {    # query not fully overlapped by reference
    ldist <- nchar(query) - mmi[1]  # include non-overlapped nucleotides in distance
  }
  return(ldist)
})


# function to compute distances between inferred and reference sequences
compute_ref_dist <- function(seq_table, ref_seqs){
  dist_mat <- outer(seq_table$sequence, ref_seqs, levDist, band = -1)
  row.names(dist_mat) <- seq_table$id
  return(dist_mat)
}


# function to collapse distances according to a grouping variable
collapse_group_dist <- function(dist_table, groups){
  group_to_ref <- t(aggregate(t(dist_table), list(groups), min))
  colnames(group_to_ref) <- group_to_ref[1, ]
  group_to_ref <- group_to_ref[-1, ]
  class(group_to_ref) <- "integer"
  return(group_to_ref)
}


# function to find sequences (or any named vector, really) with names matching any of a set of names
# get_matching_seqs <- function(seqs, patterns){
#   ref_names <- names(seqs)
#   m_ind <- sapply(patterns, function(p) str_detect(ref_names, p))
#   
#   if (length(m_ind) > 0){
#     m_ind <- apply(m_ind, 1, any)
#     return(seqs[m_ind])
#   } else return(NULL)
# }


# function to determine reference strains completely missing from a set of reads
get_missing_strains <- function(data_seqs, refs){
  if (length(refs) > 0){
    # result <- list()
    dist_mat <- outer(data_seqs, refs, levDist, band = -1)
    missing <- apply(dist_mat, 2, min) > 0
    missing <- tapply(missing, as.factor(names(missing)), all)
    missing <- names(missing[missing])
    # result[["missing"]] <- missing
    # result[["distances"]] <- dist_mat
    # return(result)
    return(missing)
  } else return(character(0))
}


# function to generate a table of "noisy"" sequences: those that are very similar to the reference sequences
countNoisySeqs <- function(seq_tab, dist_mat, min_dist = 1, max_dist = 3) {
  if (max(colSums(dist_mat == 0)) > 1) {
    print("WARNING: Your distance matrix contains seqeunces that are identical except for length differences.")
    print("This may result in inaccurate counts of noisy seuqences.")
  }
  
  if (sum(!sapply(seq_tab, class) %in% c("integer", "numeric")) > 0){
    stop("The sequence table must only contain column vectors of type 'numeric' or 'integer'.")
  }
  
  rmat <- matrix(0, nrow = ncol(seq_tab), ncol = (max_dist - min_dist) + 1, 
                 dimnames = list(colnames(seq_tab), paste0(min_dist:max_dist, "_off")))
  
  for (n in min_dist:max_dist){
    mins <- apply(dist_mat, 1, min)   
    noisy_ind <- which(dist_mat == n, arr.ind = TRUE)
    noisy_ind <- noisy_ind[dist_mat[noisy_ind] == mins[noisy_ind[, "row"]], , drop = FALSE]   # make sure it's a minimum distance
    has_zero <- apply(dist_mat[, noisy_ind[, "col"], drop = FALSE], 2, min) == 0  # is there a hit to the nearest reference?
    noisy_ind <- noisy_ind[has_zero, , drop = FALSE]
    
    if (length(noisy_ind) > 0){
      ref_row <- apply(dist_mat[, noisy_ind[, "col"] , drop = FALSE] == 0, 2, which)  # get nearest sequence that matches reference
      noisy_ind <- cbind(noisy_ind, ref_row)
      is_n_off <- ( (seq_tab[noisy_ind[, "row"], ] > 0) & 
                      (seq_tab[noisy_ind[, "row"], ] < seq_tab[noisy_ind[, "ref_row"], ]) ) * 1  # is abundance less than reference?
      row.names(is_n_off) <- row.names(noisy_ind)
      is_n_off <- aggregate(is_n_off, by = list(row.names(is_n_off)), FUN = max)
      row.names(is_n_off) <- is_n_off$Group.1
      is_n_off <- is_n_off[, -1]
      rmat[, n] <- colSums(is_n_off)
    }
    else rmat[, n] <- integer(nrow(rmat))
  }
  rmat <- cbind(rmat, "ref_noisy" =  rowSums(rmat))
  return(rmat)
}


# function that returns whether sequences of a sequence table are "noisy"
isNoisy <- function(seq_tab, dist_mat, min_dist = 1, max_dist = 3){
  if (max(colSums(dist_mat == 0)) > 1) {
    print("WARNING: Your distance matrix contains seqeunces that are identical except for length differences.")
    print("This may result in inaccurate counts of noisy seuqences.")
  }
  
  if (sum(!sapply(seq_tab, class) %in% c("integer", "numeric")) > 0){
    stop("The sequence table must only contain column vectors of type 'integer' or 'numeric'.")
  }
  
  noisy_ind <- which(dist_mat >= min_dist & dist_mat <= max_dist, arr.ind = TRUE) 
  if (length(noisy_ind) == 0){
    return(rep(FALSE, nrow(dist_mat)))
  } 

  mins <- apply(dist_mat, 1, min)   # get minimum distance for each sequence
  noisy_ind <- noisy_ind[dist_mat[noisy_ind] == mins[noisy_ind[, "row"]], , drop = FALSE]   # make sure it's a minimum distance
  if (length(noisy_ind) == 0){
    return(rep(FALSE, nrow(dist_mat)))
  }
  
  has_zero <- apply(dist_mat[, noisy_ind[, "col"], drop = FALSE], 2, min) == 0  # is there a hit to the nearest reference?
  noisy_ind <- noisy_ind[has_zero, , drop = FALSE]
  if (length(noisy_ind) == 0){
    return(rep(FALSE, nrow(dist_mat)))
  }
  
  ref_row <- apply(dist_mat[, noisy_ind[, "col"] , drop = FALSE] == 0, 2, which)   # find the nearest sequence that matches the reference
  noisy_ind <- cbind(noisy_ind, ref = ref_row)
  is_noisy <- ( (seq_tab[noisy_ind[, "row"], ] > 0) & 
                  (seq_tab[noisy_ind[, "row"], ] < seq_tab[noisy_ind[, "ref"], ]) )   # is abundance less than nearest reference abundance?
  is_noisy <- apply(is_noisy, 1, any)   # collapse to vector
  noisy_ind <- noisy_ind[is_noisy, , drop = FALSE]    # remove any with abundance higher than reference
  if (length(noisy_ind) == 0){
    return(rep(FALSE, nrow(dist_mat)))
  }
  
  noisy <- 1:nrow(dist_mat) %in% noisy_ind[, "row"]  
  return(noisy)

}


# function to annotate sequence table with 'Reference' and 'Ref_Noisy'
annotate_ref <- function(seq_table, dist_mat, sample_names, max_dist = 3){
  seq_table$ref_dist <- apply(dist_mat, 1, min)
  seq_table$reference <- seq_table$ref_dist == 0
  seq_table$ref_noisy <- isNoisy(seq_table[sample_names], dist_mat, max_dist = max_dist)
  seq_table$ref_like <- seq_table$reference | seq_table$ref_noisy
  return(seq_table)
}


# function to write sequences to a fasta file
write_fasta <- function(seq_table, fasta_path){
  seqs <- DNAStringSet(seq_table$sequence)
  names(seqs) <- seq_table$id
  writeFasta(seqs, fasta_path)
}


# function to read in a table of BLAST results
load_blast <- function(blast_path){
  blast_table <- read.table(blast_path, sep = "\t", comment.char = "#",
                            col.names = c("seqID", "seq_len", "subjectID", "subject_len", "kingdom", "sci_name", 
                                          "identity", "aln_len", "matches", "mismatches", "gapopens", "gaps", 
                                          "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
  blast_table <- as.tibble(blast_table)
  return(blast_table)
}


# function to annotate BLAST table with distance between query and subject sequence
annotate_blast_table <- function(blast_tab){
  blast_tab$lev_dist <- blast_tab$mismatches + blast_tab$gaps  # compute Levenshtein distances, not counting terminal gaps
  blast_tab$cov_dist <- pmax(blast_tab$seq_len - blast_tab$aln_len, 0)  # compute number of terminal gaps, if any
  blast_tab$tot_dist <- blast_tab$lev_dist + blast_tab$cov_dist  # sum the two to get the total distance between sequences
  return(blast_tab)
}


# function to annotate sequence table with distance to nearest NT sequence
annotate_nt_dist <- function(seq_tab, blast_tab){
  nt_dist <- tapply(blast_tab$tot_dist, INDEX = list(blast_tab$seqID), FUN = min)
  seq_tab$nt_dist <- nt_dist[seq_tab$id]
  return(seq_tab)
}


# function to identify exact matches to BLAST results
isBlastHit <- function(seq_tab, blast_tab){
  perfect_hits <- blast_tab %>% filter(tot_dist == 0)
  # perfect_hits <- blast_tab[ (blast_tab$identity == 100) & 
  #                              (blast_tab$seq_len == blast_tab$aln_len), ]
  return(seq_tab$id %in% perfect_hits$seqID)
}


# function to identify sequences that are similar, but not exact matches, to BLAST results
isBlastNoisy <- function(seq_tab, blast_tab, min_dist = 1, max_dist = 3){
  # annotate blast table to facilitate computation of differences between query and subject (nt) sequences
  # blast_tab$lev_dist <- blast_tab$mismatches + blast_tab$gaps  # compute Levenshtein distances, not counting terminal gaps
  # blast_tab$cov_dist <- pmax(blast_tab$seq_len - blast_tab$aln_len, 0)  # compute number of terminal gaps, if any
  # blast_tab$tot_dist <- blast_tab$lev_dist + blast_tab$cov_dist  # sum the two to get the total distance between sequences
  
  noisies <- blast_tab %>% group_by(seqID) %>% filter(tot_dist == min(tot_dist)) %>% ungroup() # only consider nearest hit for each sequence
  noisies <- noisies[between(noisies$tot_dist, min_dist, max_dist), ]
  is_noisy <- seq_tab$id %in% noisies$seqID
  
  hits <- isBlastHit(seq_tab, blast_tab)
  is_noisy[hits] <- FALSE
  
  return(is_noisy)
}


# function to annotate sequence table with 'Contaminant' and 'Contam_Mut'
annotate_contam <- function(seq_table, blast_table, sample_names, max_dist = 3){
  seq_table$contaminant <- isBlastHit(seq_table, blast_table) & !seq_table$ref_like
  # seq_table$contam_mut <- isBlastNoisy(seq_table, blast_table, 1, max_dist) & !seq_table$ref_like
  # seq_table$contam_noisy <- isNoisy(seq_table, contam_dist, max_dist = max_dist) & !seq_table$ref_like
  # seq_table$contam_like <- seq_table$contaminant | seq_table$contam_mut
  # seq_table$other <- !seq_table$ref_like & !seq_table$contam_like
  # seq_table$consensus <- apply(seq_table[, sample_names], 1, min) > 0
  return(seq_table)
}


# function to annotate sequence table with 'Contam_Noisy', 'Contam_Like', and 'Other'
annotate_contam_like <- function(seq_table, dist_mat, blast_table, sample_names, noisy_dist = 3, mut_dist = 1){
  seq_table$contam_dist <- apply(dist_mat, 1, min)
  seq_table$contam_noisy <- isNoisy(seq_table[sample_names], dist_mat, max_dist = noisy_dist) & !seq_table$ref_like
  seq_table$contam_mut <- isBlastNoisy(seq_table, blast_table, max_dist = mut_dist) & !seq_table$ref_like & !seq_table$contam_noisy
  seq_table$contam_like <- seq_table$contaminant | seq_table$contam_noisy | seq_table$contam_mut 
  seq_table$other <- !seq_table$ref_like & !seq_table$contaminant & !seq_table$contam_noisy
  return(seq_table)
}


# function to compute distances between sequences in the same table
compute_inter_dist <- function(seq_tab, class1, class2){
  class1_seqs <- seq_tab %>% filter(!!as.name(class1)) %>% .$sequence
  names(class1_seqs) <- seq_tab %>% filter(!!as.name(class1)) %>% .$id
  class2_seqs <- seq_tab %>% filter(!!as.name(class2)) %>% .$sequence
  names(class2_seqs) <- seq_tab %>% filter(!!as.name(class2)) %>% .$id
  class1_to_class2 <- outer(class1_seqs, class2_seqs, levDist, band = -1)
  return(class1_to_class2)
}


# function to annotate sequence table with minimum inter-sequence distance, 
# as computed above
annotate_inter_dist <- function(seq_tab, inter_dist, column_name){
  column_name <- quo_name(enquo(column_name))
  if (!column_name %in% colnames(seq_tab)){
    seq_tab <- seq_tab %>% mutate(!!column_name := rep(NA, nrow(seq_tab)))
  }
  if (length(inter_dist) == 0){
    return(seq_tab)
  }
  seq_tab[seq_tab$id %in% row.names(inter_dist), column_name] <- apply(inter_dist, 1, function(d) min(d[d > 0]))
  return(seq_tab)
}


# function to add a class label column to an annotated sequence table
annotate_class <- function(seq_table){
  seq_table$class <- factor(rep("Other", nrow(seq_table)), 
                            levels = c("Reference", "Ref Noisy", "Contaminant", "Contam Noisy", "Other"))
  seq_table$class[seq_table$reference] <- "Reference"
  seq_table$class[seq_table$ref_noisy] <- "Ref Noisy"
  seq_table$class[seq_table$contaminant] <- "Contaminant"
  seq_table$class[seq_table$contam_noisy] <- "Contam Noisy"
  # seq_table$class[seq_table$contam_mut] <- "Contam Mut"
  return(seq_table)
}


# function to assign taxonomy to sequences and add taxonomic labels to the sequence table
assign_taxonomy <- function(seq_table, genus_db_path, species_db_path){
  seqs <- seq_table[["sequence"]]
  taxa <- assignTaxonomy(seqs, refFasta = genus_db_path, multithread = TRUE, verbose = TRUE)
  taxa <- taxa %>% as.data.frame() %>% rownames_to_column(var = "sequence") %>% as.tibble()
  species <- assignSpecies(seqs, refFasta = species_db_path, verbose = TRUE)
  species <- species %>% as.data.frame() %>% rownames_to_column(var = "sequence") %>% as.tibble()
  
  taxa <- inner_join(taxa, species, by = "sequence")
  taxa <- taxa %>% mutate(Genus = ifelse(!is.na(Genus.x), Genus.x, Genus.y)) %>%
    mutate(Genus = ifelse(str_detect(Genus, Genus.y) & !is.na(Genus.y), Genus.y, Genus)) %>%
    select(-Genus.x, -Genus.y)
  taxa$taxonomy <- taxa %>% select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    apply(., 1, paste, collapse = " ")
  new_table <- inner_join(seq_table, taxa, by = "sequence")
  return(new_table)
}


# function to count the number of strains' sequences that match inferred sequences
count_strains <- function(seq_vector, strain_dist){
  seqs_to_strains <- apply(strain_dist[seq_vector > 0, ], 2, min)
  n_match <- sum(seqs_to_strains == 0)
  return(n_match)
}


# function to create a sequence summary table from a sequence table
summarize_seqs <- function(seq_table, dist_mat, refs, sample_names, max_dist){
  # summarize counts of various classes
  sample_names <- as.character(sample_names)
  exp_strains <- rep(length(unique(names(refs))), length(sample_names))
  total_count <- colSums(seq_table[, sample_names] > 0)
  strain_dist <- collapse_group_dist(dist_mat, names(refs))
  strain_count <- sapply(seq_table[, sample_names], count_strains, strain_dist)
  # strain_count <- sapply(seq_table[, sample_names], function(m) sum(colSums(strain_dist[m > 0,]) > 0))
  ref_count <- colSums(seq_table[, sample_names] > 0 & seq_table$reference)
  ref_noisy_count <- countNoisySeqs(seq_table[, sample_names], dist_mat, 1, max_dist)
  ref_like_count <- ref_count + ref_noisy_count[, "ref_noisy"]
  contam_count <- colSums(seq_table[, sample_names] > 0 & seq_table$contaminant)
  contam_noisy_count <- colSums(seq_table[, sample_names] > 0 & seq_table$contam_noisy)
  contam_like_count <- colSums(seq_table[, sample_names] > 0 & seq_table$contam_like)
  contam_mut_count <- colSums(seq_table[sample_names] > 0 & seq_table$contam_mut)
  other_count <- colSums(seq_table[, sample_names] > 0 & seq_table$other)
  
  # summarize percentage of primary classes
  ref_pct <- 100 * colSums(seq_table[seq_table$reference, sample_names]) / colSums(seq_table[, sample_names])
  ref_noisy_pct <- 100 * colSums(seq_table[seq_table$ref_noisy, sample_names]) / colSums(seq_table[, sample_names])
  contam_pct <- 100 * colSums(seq_table[seq_table$contaminant, sample_names]) / colSums(seq_table[, sample_names])
  contam_noisy_pct <- 100 * colSums(seq_table[seq_table$contam_noisy, sample_names]) / colSums(seq_table[, sample_names])
  other_pct <- 100 * colSums(seq_table[seq_table$other, sample_names]) / colSums(seq_table[, sample_names])
  
  summary_table <- data.frame(sample = sample_names, exp_strains = exp_strains, total = total_count, 
                              strains = strain_count,
                              reference = ref_count, ref_noisy_count, ref_like = ref_like_count,
                              contaminant = contam_count, contam_noisy = contam_noisy_count, 
                              contam_mut = contam_mut_count, contam_like = contam_like_count, 
                              other = other_count, 
                              pct_ref = ref_pct, pct_ref_noisy = ref_noisy_pct, 
                              pct_contam = contam_pct, pct_contam_noisy = contam_noisy_pct,
                              pct_other = other_pct,
                              row.names = NULL, check.names = FALSE) #%>% as.tibble
  return(summary_table)
}


# function to gather sample counts into a single column
gather_samples <- function(seq_table, sample_names, sample_colname){
  sample_colname <- enquo(sample_colname)
  gg_table <- gather(seq_table, !!sample_colname, "count", one_of(sample_names)) %>%
    select(id, !!sample_colname, count, everything())
  return(gg_table)
}


# function to create a taxa summary table from a sequence table
summarize_taxa <- function(seq_table, sample_names, group){
  group <- enquo(group)
  seq_table_gg <- gather_samples(seq_table , sample_names, !!group)
  # seq_table_gg$method <- factor(seq_table_gg$method, 
  #                               levels = c("uclust", "uparse", "med", "unoise", "deblur", "dada2"),
  #                               labels = c("UCLUST", "UPARSE", "MED", "UNOISE", "Deblur", "DADA2"))
  taxa_summary <- seq_table_gg %>% 
    group_by(!!group, class) %>% 
    filter(count > 0) %>% 
    summarise(n_taxa = length(unique(taxonomy))) %>%
    complete(class, fill = list(n_taxa = 0)) %>% ungroup()
  taxa_summary <- spread(taxa_summary, class, n_taxa)
  totals <- seq_table_gg %>% 
    group_by(!!group) %>%
    filter(count > 0) %>%
    summarise(n_taxa = length(unique(taxonomy)))
  taxa_summary <- bind_cols(Total = totals$n_taxa, taxa_summary) %>%
    select(!!group, everything())
  return(taxa_summary)
}


# function to add a sanity check to the summary table
sanity_check_summary <- function(sum_table){
  sum_table <- sum_table %>% 
    mutate(check_sum = reference + ref_noisy + contaminant + contam_noisy + other,
           pct_check = pct_ref + pct_ref_noisy + pct_contam+ pct_contam_noisy + pct_other) %>%
    select(-starts_with("pct"), everything())
  return(sum_table)
}


# function to convert list of method summaries to list of sample summaries
# method_to_sample <- function(method_table_list, mt_colname, st_colname){
#   sample_names <- method_table_list[[1]][[mt_colname]]
#   sample_summary <- sapply(sample_names, function(m) NULL)
#   
#   sample_summary <- mapply(function(ssum, sname, msum){
#     ssum <- lapply(msum, function(ms, ss, sn) {
#       ss <- rbind(ss, ms[ms[[mt_colname]] == sn, ])
#       return(ss)
#     }, ss = ssum, sn = sname)
#     
#     ssum <- do.call("rbind", ssum) %>% rownames_to_column(var = st_colname)
#     ssum <- ssum %>% select(-one_of(mt_colname))
#     return(ssum)
#   }, ssum = sample_summary, sname = sample_names, MoreArgs = list(msum = method_table_list), 
#   SIMPLIFY = FALSE)
#  
#   return(sample_summary) 
# }


# function to convert list of tables listed by one variable to a list of tables listed by a different variable
transpose_table_list <- function(old_list, old_id_col, new_id_col){
  new_list_names <- as.character(old_list[[1]][[old_id_col]])
  new_row_ids <- names(old_list)
  new_list <- list()
  
  for (nn in new_list_names){
    new_table <- data.frame()
    for (nid in new_row_ids){
      temp <- old_list[[nid]]
      new_row <- data.frame(nid, temp[temp[old_id_col] == nn, ], check.names = FALSE)
      new_table <- rbind(new_table, new_row)
    }
    colnames(new_table)[1] <- new_id_col
    new_table <- new_table %>% select(-one_of(old_id_col)) %>% as.tibble()
    new_list[[nn]] <- new_table
  }
  return(new_list)
}


# function to write a list of tables to a single file
write_tables <- function(table_list, file_path, compact = FALSE, overwrite = TRUE){
  if (file.exists(file_path) & overwrite == FALSE){
    cat("File", file_path, "exists, I will not overwrite it.\n")
    return(FALSE)
  } else if (file.exists(file_path) & overwrite == TRUE){
    file.remove(file_path)
  }
  
  header = paste(colnames(table_list[[1]]), collapse = "\t")
  if (compact) write(header, file_path, append = TRUE)
  tab_names <- names(table_list)
  for (i in seq_along(table_list)){
    if (!compact){
      write(c("\n", tab_names[i]), file_path, append = TRUE)
      write(header, file_path, append = TRUE)
    }
    write_tsv(table_list[[i]], file_path, append = TRUE)
  }
  
  return(TRUE)
}


# function to compute precision and recall at the sequence level
compute_pr_seqs <- function(sum_table){
  seq_stats <- sum_table %>% 
    mutate(TP = reference,
           FN = pmax(exp_strains - reference, 0),
           FP = ref_noisy + contaminant + contam_noisy + other,
           FP_NC = ref_noisy) %>%
    select(1, TP, FN, FP, FP_NC)
  
  seq_stats <- seq_stats %>%
    mutate(recall = 100 * TP / (TP + FN),
           precision = 100 * TP / (TP + FP),
           precision_NC = 100 * TP / (TP + FP_NC)) %>%
    as.tibble()
  
  return(seq_stats)
}


# function to compute precision and recall at the read level
compute_pr_reads <- function(seq_tab, sample_names, sample_colname = "sample"){
  # count the number of true positive, false negative, and false positive reads
  TP <- seq_tab %>% filter(reference | contaminant) %>% select(sample_names) %>% colSums()
  FN <- seq_tab %>% filter(ref_noisy | contam_noisy) %>% select(sample_names) %>% colSums()
  FP <- seq_tab %>% filter(other) %>% select(sample_names) %>% colSums()
  TP_ref <- seq_tab %>% filter(reference) %>% select(sample_names) %>% colSums()
  FN_ref <- seq_tab %>% filter(ref_noisy) %>% select(sample_names) %>% colSums()
  
  read_stats <- tibble(sample = sample_names,
                       recall = 100 * TP / (TP + FN),
                       precision = 100 * TP / (TP + FP),
                       recall_ref = 100 * TP_ref / (TP_ref + FN_ref),
                       precision_ref = 100 * TP_ref / (TP_ref + FP))
  colnames(read_stats)[1] <- sample_colname
 
  return(read_stats)
}


# function to compute percentage of sample reads inferred as reference
compute_ref_perc <- function(seq_tab, sample_names){
  ref_reads <- seq_tab %>% filter(reference) %>% select(sample_names) %>% colSums()
  all_reads <- seq_tab %>% select(sample_names) %>% colSums()
  ref_percent <- 100 * ref_reads / all_reads
  names(ref_percent) <- sample_names
  return(ref_percent)
}


# function to annotate normalized sequence counts
annotate_norms <- function(all_seq_table, group){
  all_seq_table <- all_seq_table %>%
    mutate(log10_count = log10(count)) %>%
    group_by_at(group) %>%
    mutate(total_reads = sum(count),
           rel_count = count / total_reads) %>% 
    ungroup()

  if (identical(group, "sample")){
    all_seq_table <- all_seq_table %>%
      mutate(norm_median = round(rel_count * median(total_reads)),
             log10_norm_med = log10(norm_median))
  }
  
  all_seq_table <- all_seq_table %>% select(-sequence, -total_reads, sequence)
  return(all_seq_table)
}

###############################################################################
#
# Plotting functions
#
###############################################################################

theme_set(theme_bw())

# function to define general theme parameters
big_labels <- function(title = 10, text = 9, angle = 45, hjust = 1, vjust = 1, legend.position = "top", key.size = 1){
  theme(axis.text.x = element_text(size = text, angle = angle, hjust = hjust, vjust = vjust),
        axis.text.y = element_text(size = text),
        axis.title.x = element_text(margin = margin(t = 6)),
        axis.title.y = element_text(margin = margin(r = 6)),
        plot.title = element_text(face = "bold", size = title, hjust = 0.5),
        plot.subtitle = element_text(size = title - 2, hjust = 0.5),
        axis.title = element_text(face = "bold", size = title),
        strip.text = element_text(size = title),
        legend.position = legend.position,
        legend.justification = "center",
        legend.title = element_text(face = "bold", size = title, margin = margin(r = 7)),
        legend.text = element_text(size = text, margin = margin(l = 0, r = 7), hjust = 0),
        # legend.spacing.x = unit(1.5,"cm"),
        legend.key.size = unit(key.size, "lines"),
        legend.margin = margin(b = 0, t = 0),
        legend.box.margin = margin(0, 0, 0, 0),
        panel.spacing = unit(0.25, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        plot.margin = margin(2, 2, 2, 2)) 
}

# ref vs. non-ref beeswarm plots, log10 raw counts
snr_beeswarm_plot <- function(gg_table, x_var, x_label, point_colors = c("blue", "orange"), title){
  
  snr_bees <- ggplot(data = gg_table %>% filter(count > 0), aes_string(x = x_var, y = "log10_count")) +
    geom_beeswarm(aes(color = factor(reference, levels = c(T, F))), size = 2.5, priority = "ascending", cex = 0.15, dodge.width = 0.5) +
    labs(title = title, x = x_label, y = "log10(count)") +
    scale_color_manual(name = "Sequence type", labels = c("Reference", "non-Reference"), values = point_colors) +
    theme(legend.position = "right",
          legend.justification = c(1, 0.5), 
          axis.title = element_text(face = "bold"))
    
  return(snr_bees)
}


# create dodged bar plot of sequence counts for a given label
dodged_bar_plot <- function(gg_table, x_axis, seq_label = NULL, bar_colors, x_labels, psub = NULL){
  if (is.null(seq_label)){
    gg_data <- gg_table %>% filter(count > 0)
    ptitle <- "Total inferred sequences"
  }
  else {
    gg_data <- gg_table %>% filter(class == seq_label, count > 0)
    ptitle <- paste(seq_label, "sequences")
  }
  
  labels <- paste(levels(gg_table$method), "   ")
  
  seq_bars <- ggplot(gg_data, aes_string(x = x_axis)) +
    geom_bar(aes(fill = method), 
             width = 0.8, position = position_dodge(preserve = "single")) +
    scale_fill_manual(name = "Method   ", values = bar_colors, labels = labels) +
    scale_x_discrete(labels = x_labels) +
    labs(title = ptitle, 
         subtitle = psub, x = "dilution", y = "sequences") +
    theme(axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
          plot.subtitle = element_text(size = 16, hjust = 0.5),
          axis.title = element_text(face = "bold", size = 16),
          legend.position = "bottom",
          legend.justification = 0.5,
          legend.title = element_text(face = "bold", size = 16),
          legend.text = element_text(size = 16, hjust = 1),
          legend.key.size = unit(1, "lines")) +
    guides(fill = guide_legend(nrow = 1))
  return(seq_bars)
}


# create dodged column plot of taxa counts for a given label
dodged_col_plot <- function(gg_table, x_axis, seq_label = NULL, bar_colors, x_labels, psub = NULL){
  if (is.null(seq_label)){
    ptitle <- "Total inferred taxa"
    gg_data <- gg_table
  }
  else {
    gg_data <- gg_table %>% filter(class == seq_label)
    ptitle <- paste(seq_label, "taxa")
  }
  
  labels <- paste(levels(gg_table$method), "   ")
  
  tax_cols <- ggplot(gg_data, aes_string(x = x_axis)) +
    geom_col(aes(y = taxa, fill = method), 
             width = 0.8, position = "dodge") +
    scale_fill_manual(name = "Method   ", values = bar_colors, labels = labels) +
    # scale_x_discrete(labels = x_labels) +
    labs(title = ptitle, 
         subtitle = psub, x = "dilution", y = "taxa") +
    guides(fill = guide_legend(nrow = 1))
  return(tax_cols)
}


# stacked bar plots of read counts in each class for each method
class_comp_plot <- function(gg_table, bar_colors){
  class_comp_bars <- ggplot(data = gg_table) +
    geom_col(aes(x = method, y = count, 
                 fill = factor(class, levels = c("Other", "Contam Noisy", "Ref Noisy", "Contaminant", "Reference"))), 
             width = 0.5, position = position_fill()) +
    scale_fill_manual(name = "Classification", values = bar_colors) +
    labs(title = "", x = "method", y = "relative abundance") +
    theme(axis.title = element_text(face = "bold"))
  return(class_comp_bars)
}


# plot sequence lines across dilution series
dilution_line_plot <- function(gg_table, group_var, filter_var = NULL, stat, size = 0.5){
  filter_var <- enquo(filter_var)

  ggplot(gg_table %>% filter(!! filter_var, count > 0), aes(x = sample)) +
    geom_line(aes_string(group = group_var, color = group_var), stat = stat, size = size) +
    geom_point(aes_string(color = group_var), stat = stat, size = 3 * size) +
    scale_x_discrete(labels = dilution_labels) +
    xlab("dilution")
}


# plot taxon lines across dilution series
dilution_tax_line_plot <- function(gg_table, group_var, class_var = NULL, stat, size = 0.5){
  if (!is.null(class_var)){
    gg_table <- gg_table %>% filter(class %in% class_var)
  }

  ggplot(gg_table, aes(x = sample)) +
    geom_line(aes_string(group = group_var, color = group_var, y = "taxa"), stat = stat, size = size) +
    geom_point(aes_string(color = group_var, y = "taxa"), stat = stat, size = 3 * size) +
    scale_x_discrete(labels = dilution_labels) +
    xlab("dilution")
}