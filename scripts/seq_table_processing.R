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
load_uclust <- function(otu_file_path, seqs_file_path, sample_names = NULL, rm_singletons = TRUE, min_abund = 2){
  uclust_table <- as.tibble(read.table(file = otu_file_path, header = TRUE, sep = "\t", 
                                       skip = 1, comment.char = ""))
  colnames(uclust_table)[1] <- "id"
  uclust_table$id <- str_replace(uclust_table$id, "denovo", "denovo_")  # make ids more readable
  
  if (is.null(sample_names)){
    sample_names <- colnames(uclust_table)[-1]  
    sample_names <- sample_names[order(sample_names)]
  }
  uclust_table <- uclust_table %>% select(id, sample_names)
  
  # if (rm_singletons){
  #   uclust_table <- uclust_table[apply(uclust_table[, sample_names], 1, max) > 1, ]  # remove singletons
  # }
  uclust_table <- uclust_table[rowSums(uclust_table[, sample_names]) >= min_abund, ]  # remove sequences with less than min_abund reads across samples
  
  uclust_seqs <- readDNAStringSet(seqs_file_path)
  uclust_seqs <- as.character(uclust_seqs)
  names(uclust_seqs) <- sapply(str_split(names(uclust_seqs), " "), `[`, 1)  # remove annotation from names
  names(uclust_seqs) <- str_replace(names(uclust_seqs), "denovo", "denovo_")  # modify names to match table ids
  uclust_table$sequence <- uclust_seqs[uclust_table$id]  # add sequences to table
  
  return(uclust_table)
}


# function to load UPARSE OTU table
load_uparse <- function(otu_file_path, seqs_file_path, sample_names = NULL, rm_singletons = TRUE, min_abund = 2){
  uparse_table <- as.tibble(read.table(file = otu_file_path, header = TRUE, sep = "\t", 
                                       comment.char = ""))
  colnames(uparse_table)[1] <- "id"
  uparse_table$id <- str_replace(uparse_table$id, "OTU", "OTU_")
  
  if (is.null(sample_names)){
    sample_names <- colnames(uparse_table)[-1]  
    sample_names <- sample_names[order(sample_names)]
  } else colnames(uparse_table)[-1] <- sample_names
  
  if (rm_singletons){
    uparse_table <- uparse_table[apply(uparse_table[, sample_names], 1, max) > 1, ]  # remove singletons
  }
  
  uparse_table <- uparse_table[rowSums(uparse_table[, sample_names]) >= min_abund, ]  # remove sequences with less than min_abund reads across samples

  uparse_seqs <- readDNAStringSet(seqs_file_path)
  uparse_seqs <- as.character(uparse_seqs)
  names(uparse_seqs) <- str_replace(names(uparse_seqs), "OTU", "OTU_")
  uparse_table$sequence <- uparse_seqs[uparse_table$id]  # add sequences to table
  
  return(uparse_table)
}


# function to load UNOISE ZOTU table
load_unoise <- function(otu_file_path, seqs_file_path, sample_names = NULL, rm_singletons = TRUE, min_abund = 2){
  unoise_table <- as.tibble(read.table(file = otu_file_path, header = TRUE, sep = "\t", 
                                       comment.char = ""))
  colnames(unoise_table)[1] <- "id"
  unoise_table$id <- str_replace(unoise_table$id, "OTU", "ZOTU_")
  
  if (is.null(sample_names)){
    sample_names <- colnames(unoise_table)[-1]  
    sample_names <- sample_names[order(sample_names)]
  } else colnames(unoise_table)[-1] <- sample_names
  
  if (rm_singletons){
    unoise_table <- unoise_table[apply(unoise_table[, sample_names], 1, max) > 1, ]  # remove singletons
  }
  
  unoise_table <- unoise_table[rowSums(unoise_table[, sample_names]) >= min_abund, ]  # remove sequences with less than min_abund reads across samples
  
  # load UNOISE sequences
  unoise_seqs <- readDNAStringSet(seqs_file_path)
  unoise_seqs <- as.character(unoise_seqs)
  names(unoise_seqs) <- sapply(names(unoise_seqs), FUN = str_replace, "OTU", "ZOTU_")
  unoise_table$sequence <- unoise_seqs[unoise_table$id]  # add sequences to table
  
  return(unoise_table)
}


# function to load MED Node table
load_med <- function(otu_file_path, seqs_file_path, chimera_file_path, sample_names = NULL, rm_singletons = TRUE){
  med_table <- read.table(file = otu_file_path, header = TRUE, sep = "\t")
  row.names(med_table) <- med_table[, 1]
  med_table <- med_table[, -1]  # remove sample name column
  med_table <- data.frame(t(med_table))   # samples as columns
  med_table <- as.tibble(rownames_to_column(med_table, var = "id"))
  med_table$id <- str_replace(med_table$id, "X", "Node_")
  
  if (is.null(sample_names)){
    sample_names <- colnames(med_table)[-1]  
    sample_names <- sample_names[order(sample_names)]
  } else colnames(med_table)[-1] <- sample_names
  
  # load MED sequences
  med_seqs <- readDNAStringSet(seqs_file_path)
  med_seqs <- as.character(med_seqs)
  med_seqs <- sapply(med_seqs, FUN = str_replace, "-+", "")
  names(med_seqs) <- paste0("Node_", sapply(str_split(names(med_seqs), "\\|"), FUN = `[`, 1))  # remove size annotation
  
  # remove chimeras detected with UCHIME method
  med_seqs.chimeras <- readDNAStringSet(chimera_file_path)
  if (length(med_seqs.chimeras) > 0){
    med_seqs.chimeras <- as.character(med_seqs.chimeras)
    names(med_seqs.chimeras) <- paste0("Node_", sapply(str_split(names(med_seqs.chimeras), ";"), FUN = `[`, 1))
    med_table <- med_table[!med_seqs %in% med_seqs.chimeras, ]
    med_table
    med_seqs <- med_seqs[!med_seqs %in% med_seqs.chimeras]
  }
  
  if (rm_singletons){
    med_table <- med_table[apply(med_table[, sample_names], 1, max) > 1, ]  # remove singletons
  }

  med_table$sequence <- med_seqs[med_table$id]  # add sequences to table
  
  return(med_table)
}


# function to load Deblur sOTU table
load_deblur <- function(otu_file_path, sample_names = NULL, rm_singletons = TRUE){
  deblur_table <- as.tibble(read.table(file = otu_file_path, header = TRUE, sep = "\t", 
                                       skip = 1, comment.char = ""))
  colnames(deblur_table)[1] <- "sequence"
  deblur_table$sequence <- as.character(deblur_table$sequence)
  deblur_table$id <- paste0("sOTU_", 1:nrow(deblur_table))
  
  if (is.null(sample_names)){
    sample_names <- colnames(deblur_table)[c(-1, -length(colnames(deblur_table)))]  
    sample_names <- sample_names[order(sample_names)]
  } else colnames(deblur_table)[c(-1, -length(colnames(deblur_table)))] <- sample_names
  
  deblur_table <- deblur_table %>% select(id, sample_names, sequence)
  
  if (rm_singletons){
    deblur_table <- deblur_table[apply(deblur_table[, sample_names], 1, max) > 1, ]  # remove singletons
  }
  
  return(deblur_table)  
}


# function to load DADA2 ASV table
load_dada2 <- function(otu_file_path, sample_names = NULL, rm_singletons = TRUE, min_abund = 2){
  dada2_table <- read.table(file = otu_file_path, header = TRUE, sep = "\t")
  
  if (length(dada2_table) == 1){
    colnames(dada2_table) <- ifelse(is.null(sample_names), "sample_1", sample_names)
  } else {
    dada2_table <- data.frame(t(dada2_table))
  }
  
  dada2_table <- as.tibble(rownames_to_column(dada2_table, var = "sequence"))
  dada2_table$id <- paste0("ASV_", 1:nrow(dada2_table))
  
  if (is.null(sample_names)){
    sample_names <- colnames(dada2_table)[c(-1, -length(colnames(dada2_table)))]  
    sample_names <- sample_names[order(sample_names)]
  } else colnames(dada2_table)[c(-1, -length(colnames(dada2_table)))] <- sample_names
  
  dada2_table <- dada2_table %>% select(id, sample_names, sequence)
  
  if (rm_singletons){
    dada2_table <- dada2_table[apply(dada2_table[, sample_names], 1, max) > 1, ]  # remove singletons
  }
  
  dada2_table <- dada2_table[rowSums(dada2_table[, sample_names]) >= min_abund, ]  # remove sequences with less than min_abund reads across samples
  
  return(dada2_table)
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
compute_dist_to_ref <- function(seq_table, ref_seqs){
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
isRefNoisy <- function(seq_tab, dist_mat, min_dist = 1, max_dist = 3){
  if (max(colSums(dist_mat == 0)) > 1) {
    print("WARNING: Your distance matrix contains seqeunces that are identical except for length differences.")
    print("This may result in inaccurate counts of noisy seuqences.")
  }
  
  if (sum(!sapply(seq_tab, class) %in% c("integer", "numeric")) > 0){
    stop("The sequence table must only contain column vectors of type 'integer' or 'numeric'.")
  }
  
  mins <- apply(dist_mat, 1, min)   # get minimum distance for each sequence
  noisy_ind <- which(dist_mat >= min_dist & dist_mat <= max_dist, arr.ind = TRUE) 
  if (length(noisy_ind) == 0){
    return(rep(FALSE, nrow(dist_mat)))
  } 
  else {
    noisy_ind <- noisy_ind[dist_mat[noisy_ind] == mins[noisy_ind[, "row"]], , drop = FALSE]   # make sure it's a minimum distance
    has_zero <- apply(dist_mat[, noisy_ind[, "col"], drop = FALSE], 2, min) == 0  # is there a hit to the nearest reference?
    noisy_ind <- noisy_ind[has_zero, , drop = FALSE]
    ref_row <- apply(dist_mat[, noisy_ind[, "col"] , drop = FALSE] == 0, 2, which)   # find the nearest sequence that matches the reference
    noisy_ind <- cbind(noisy_ind, ref = ref_row)
    is_noisy <- ( (seq_tab[noisy_ind[, "row"], ] > 0) & 
                    (seq_tab[noisy_ind[, "row"], ] < seq_tab[noisy_ind[, "ref"], ]) )   # is abundance less than reference?
    is_noisy <- apply(is_noisy, 1, any)   # collapse to vector
    noisy_ind <- noisy_ind[is_noisy, , drop = FALSE]    # remove any with abundance higher than reference
    noisy <- 1:nrow(dist_mat) %in% noisy_ind[, "row"]  
    return(noisy)
  }
}


# function to annotate sequence table with 'Reference' and 'Ref_Noisy'
annotate_ref <- function(seq_table, dist_mat, sample_names, max_dist = 3){
  seq_table$dist_to_ref <- apply(dist_mat, 1, min)
  seq_table$reference <- seq_table$dist_to_ref == 0
  seq_table$ref_noisy <- isRefNoisy(seq_table[, sample_names], dist_mat, 1, max_dist)
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


# function to identify exact matches to BLAST results
isBlastHit <- function(seq_tab, blast_tab){
  perfect_hits <- blast_tab[ (blast_tab$identity == 100) & 
                               (blast_tab$seq_len == blast_tab$aln_len), ]
  return(seq_tab$id %in% perfect_hits$seqID)
}


# function to identify sequences that are similar, but not exact matches, to BLAST results
isBlastNoisy <- function(seq_tab, blast_tab, min_off = 1, max_off = 3){
  # annotate blast table to facilitate computation of differences between query and subject (nt) sequences
  blast_tab$lev_dist <- blast_tab$mismatches + blast_tab$gaps  # compute Levenshtein distances, not counting terminal gaps
  blast_tab$cov_dist <- pmax(blast_tab$seq_len - blast_tab$aln_len, 0)  # compute number of terminal gaps, if any
  blast_tab$tot_dist <- blast_tab$lev_dist + blast_tab$cov_dist  # sum the two to get the total distance between sequences
  
  noisies <- blast_tab %>% group_by(seqID) %>% filter(tot_dist == min(tot_dist)) %>% ungroup() # only consider nearest hit for each sequence
  # noisies <- blast_tab[between(blast_tab$tot_dist, min_off, max_off), ]
  noisies <- noisies[between(noisies$tot_dist, min_off, max_off), ]
  is_noisy <- seq_tab$id %in% noisies$seqID
  
  hits <- isBlastHit(seq_tab, blast_tab)
  is_noisy[hits] <- FALSE
  
  return(is_noisy)
}


# function to annotate sequence table with 'Contaminant' and 'Contam_Noisy'
annotate_contam <- function(seq_table, blast_table, sample_names, max_off = 3){
  seq_table$contaminant <- isBlastHit(seq_table, blast_table) & !seq_table$ref_like
  seq_table$contam_noisy <- isBlastNoisy(seq_table, blast_table, 1, max_off) & !seq_table$ref_like
  seq_table$contam_like <- seq_table$contaminant | seq_table$contam_noisy
  seq_table$other <- !seq_table$ref_like & !seq_table$contam_like
  seq_table$consensus <- apply(seq_table[, sample_names], 1, min) > 0
  return(seq_table)
}


# function to create a summary table from a sequence table
summarize_seqs <- function(seq_table, dist_mat, sample_names, strains, max_dist){
  # summarize counts of various classes
  exp_strains <- rep(length(unique(strains)), length(sample_names))
  total_count <- colSums(seq_table[, sample_names] > 0)
  strain_dist <- collapse_group_dist(dist_mat, strains)
  strain_count <- sapply(seq_table[, sample_names], function(m) sum(colSums(strain_dist[m > 0,]) > 0))
  ref_count <- colSums(seq_table[, sample_names] > 0 & seq_table$reference)
  ref_noisy_count <- countNoisySeqs(seq_table[, sample_names], dist_mat, 1, max_dist)
  ref_like_count <- ref_count + ref_noisy_count[, "ref_noisy"]
  contam_count <- colSums(seq_table[, sample_names] > 0 & seq_table$contaminant)
  contam_noisy_count <- colSums(seq_table[, sample_names] > 0 & seq_table$contam_noisy)
  contam_like_count <- colSums(seq_table[, sample_names] > 0 & seq_table$contam_like)
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
                              contam_like = contam_like_count, other = other_count, 
                              pct_ref = ref_pct, pct_ref_noisy = ref_noisy_pct, 
                              pct_contam = contam_pct, pct_contam_noisy = contam_noisy_pct,
                              pct_other = other_pct,
                              row.names = NULL, check.names = FALSE) #%>% as.tibble
  return(summary_table)
}


# function to add a sanity check to the summary table
sanity_check_summary <- function(sum_table){
  sum_table <- sum_table %>% 
    mutate(check_sum = reference + ref_noisy + contaminant + contam_noisy + other,
           pct_check = pct_ref + pct_ref_noisy + pct_contam+ pct_contam_noisy + pct_other) %>%
    select(-starts_with("pct"), everything())
  return(sum_table)
}


# function to compute precision and recall at the sequence level
compute_pr_seqs <- function(sum_table){
  seq_stats <- sum_table %>% 
    mutate(TP = reference,
           FN = exp_strains - reference,
           FP = ref_noisy + contaminant + contam_noisy + other) %>%
    select(1, TP, FN, FP)
  
  seq_stats <- seq_stats %>%
    mutate(recall = 100 * TP / (TP + FN),
           precision = 100 * TP / (TP + FP)) %>%
    as.tibble()
  
  return(seq_stats)
}


# function to gather sample counts into a single column
gather_samples <- function(seq_table, sample_names, sample_colname){
  gg_table <- gather(seq_table, !!sample_colname, "count", one_of(sample_names)) %>%
    select(id, !!sample_colname, count, everything())
  return(gg_table)
}


# function to annotate normalized sequence counts
annotate_norms <- function(all_seq_table, group){
  all_seq_table <- all_seq_table %>%
    mutate(log10_count = log10(count)) %>%
    group_by_at(group) %>%
    mutate(total_reads = sum(count),
           rel_count = count / total_reads) %>% 
    ungroup()

  if (group == "sample"){
    all_seq_table <- all_seq_table %>%
      mutate(norm_median = round(rel_count * median(total_reads)),
             log10_norm_med = log10(norm_median))
  }
  
  all_seq_table <- all_seq_table %>% select(-sequence, -total_reads, sequence)
  return(all_seq_table)
}


# function to compute precision and recall at the read level
compute_pr_reads <- function(all_seq_tab_gg, group){
  # count the number of true positive, false negative, and false positive reads
  all_seq_tab_gg <- all_seq_tab_gg %>% group_by_at(group)
  TP <- all_seq_tab_gg %>% filter(reference | contaminant) %>% summarize(TP = sum(count))
  FN <- all_seq_tab_gg %>% filter(ref_noisy | contam_noisy) %>% summarize(FN = sum(count))
  FP <- all_seq_tab_gg %>% filter(other) %>% summarize(FP = sum(count))
  TP_ref <- all_seq_tab_gg %>% filter(reference) %>% summarize(TP_ref = sum(count))
  FN_ref <- all_seq_tab_gg %>% filter(ref_noisy) %>% summarize(FN_ref = sum(count))
  
  read_stats <- TP %>% inner_join(FN, by = group) %>% 
    inner_join(FP, by = group) %>%
    inner_join(TP_ref, by = group) %>%
    inner_join(FN_ref, by = group)
  
  read_stats <- read_stats %>% mutate(recall = 100 * TP / (TP + FN),
                                      precision = 100 * TP / (TP + FP),
                                      recall_ref = 100 * TP_ref / (TP_ref + FN_ref),
                                      precision_ref = 100 * TP_ref / (TP_ref + FP)) %>%
                                      # false_neg_rate = 100 * FN / (TP + FN),
                                      # false_disc_rate = 100 * FP / (TP + FP)) %>%
    select(-TP, -FN, -FP, -TP_ref, -FN_ref)
  
  return(read_stats)
}
