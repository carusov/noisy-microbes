#!/bin/bash

### Author: Vincent Caruso
### Date: 7/13/2017
### This script implements the UCLUST pipeline through the QIIME SOP
### (see www.qiime.org for detailed documentation)

DATA=~/projects/thesis/data

# First merge forward and reverse reads from each sample, and pool the merged reads
usearch -fastq_mergepairs $DATA/raw/*R1.fastq \
        -fastqout $DATA/clean/pooled_merged.fastq \
	-fastq_maxdiffs 100 \
        -fastq_pctid 50 \
        -fastq_minmergelen 245 \
        -fastq_maxmergelen 255 \
	-relabel @ \
        -report $DATA/reports/merge_report.txt

# Create report of the expected errors of the merged reads
usearch -fastq_eestats2 $DATA/clean/pooled_merged.fastq \
        -output $DATA/reports/merged_eestats.txt \
	-length_cutoffs 200,*,10 \

# Quality filter the reads using a maximum of 2.0 expected errors
usearch -fastq_filter $DATA/clean/pooled_merged.fastq \
        -fastqout $DATA/clean/pooled_filtered.fastq \
	-fastq_maxee 2.0

# Reformat sequence files to QIIME format
usearch -fastq_filter $DATA/clean/pooled_filtered.fastq \
        -fastaout $DATA/clean/pooled_filtered.fasta \
        -fastq_ascii 64

### The rest of the pipeline is executed with QIIME scripts

# Identify chimeric sequences 
identify_chimeric_seqs.py -m usearch61 \
        -i $DATA/clean/pooled_filtered.fasta \
        -r $DATA/references/gold.fa \
        -o $DATA/uchime/

# Filter out identified chimeric sequences
filter_fasta.py -f $DATA/clean/pooled_filtered.fasta \
        -o $DATA/clean/pooled_nochim.fa \
        -s $DATA/uchime/chimeras.txt \
        -n 

# Pick de novo OTUs
# I need to investigate the parameters that this QIIME script passes to UCLUST. Are sequences sorted by size?
pick_de_novo_otus.py -i $DATA/clean/pooled_nochim.fa \
        -o $DATA/uclust/ \
        -f

# Dereplicate 
# usearch -fastx_uniques clean/all_filtered.fastq -fastqout clean/all_uniques.fastq \
#	-sizeout -relabel Uniq

# usearch -cluster_fast clean/s160_MC_Neat_filtered.fastq -id 0.97 \
#	-centroids clustered/centroids.fasta -uc clustered/clusters.uc \
#	-sort size -maxaccepts 1 -maxrejects 80
