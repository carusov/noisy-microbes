#!/bin/bash

### Author: Vincent Caruso
### Date: 7/13/2017
### This script implements the UCLUST pipeline through the QIIME SOP
### (see www.qiime.org for detailed documentation)

DATA=~/projects/thesis/data
RESULTS=~/projects/thesis/results

# Activate the qiime environment in miniconda
source activate qiime1

# Convert .fastq to .fasta for use with QIIME
printf "\nGenerating .fasta file from .fastq file\n"
usearch -fastq_filter $DATA/clean/pooled_filtered.fastq \
        -fastaout $DATA/clean/pooled_filtered.fasta

# Reformat sequence files to QIIME format
printf "\nReformatting read sequences to QIIME format...\n"
sed '/^>/ s/\./_/' $DATA/clean/pooled_filtered.fasta \
    > $DATA/clean/pooled_filtered_qiime.fasta

### The rest of the pipeline is executed with QIIME scripts

if [ ! -d "$RESULTS/uclust" ]; then
    mkdir $RESULTS/uclust
fi

# Identify chimeric sequences 
printf "\nIdentifying chimeric sequences...\n"
identify_chimeric_seqs.py -m usearch61 \
        -i $DATA/clean/pooled_filtered_qiime.fasta \
        -r $DATA/reference/gold.fa \
        -o $RESULTS/uclust/

# Filter out identified chimeric sequences
printf "\nRemoving chimeric sequences from sample reads...\n"
filter_fasta.py -f $DATA/clean/pooled_filtered_qiime.fasta \
        -o $RESULTS/uclust/pooled_nochim.fa \
        -s $RESULTS/uclust/chimeras.txt \
        -n 

# Pick de novo OTUs
# I need to investigate the parameters that this QIIME script passes to UCLUST.
# Are sequences sorted by size?
printf "\nPicking OTUs (de novo) and representative seqeunces, assigning\n"
printf "taxonomy, performing multiple alignments, and building a\n"
printf "phylogenetic tree...\n"
pick_de_novo_otus.py -i $RESULTS/uclust/pooled_nochim.fa \
        -o $RESULTS/uclust/ \
        -f

# Convert the .biom file to a tabbed .txt file
biom convert -i $RESULTS/uclust/otu_table.biom \
     -o $RESULTS/uclust/otu_table.txt \
     --to-tsv
