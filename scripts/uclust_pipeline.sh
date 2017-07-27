#!/bin/bash

### Author: Vincent Caruso
### Date: 7/13/2017
### This script implements the UCLUST pipeline through the QIIME SOP
### (see www.qiime.org for detailed documentation)

DATA=~/projects/thesis/data

printf "\nReformatting read sequences to QIIME format...\n"

# Reformat sequence files to QIIME format
usearch -fastq_filter $DATA/clean/pooled_filtered.fastq \
        -fastaout $DATA/clean/pooled_filtered.fasta \
        -fastq_ascii 64

### The rest of the pipeline is executed with QIIME scripts

printf "\nIdentifying chimeric sequences...\n"

# Identify chimeric sequences 
identify_chimeric_seqs.py -m usearch61 \
        -i $DATA/clean/pooled_filtered.fasta \
        -r $DATA/reference/gold.fa \
        -o $DATA/uclust/

printf "\nRemoving chimeric sequences from sample reads...\n"

# Filter out identified chimeric sequences
filter_fasta.py -f $DATA/clean/pooled_filtered.fasta \
        -o $DATA/uclust/pooled_nochim.fa \
        -s $DATA/uclust/chimeras.txt \
        -n 

printf "\nPicking OTUs (de novo) and representative seqeunces, assigning"
printf "taxonomy, performing multiple alignments, and building a"
printf "phylogenetic tree..."

# Pick de novo OTUs
# I need to investigate the parameters that this QIIME script passes to UCLUST.
# Are sequences sorted by size?
pick_de_novo_otus.py -i $DATA/uclust/pooled_nochim.fa \
        -o $DATA/uclust/ \
        -f
