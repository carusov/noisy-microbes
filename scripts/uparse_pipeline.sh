#!/bin/bash

### Author: Vincent Caruso
### Date: 7/1/2017
### This script implements the UPARSE pipeline according to its recommended usage
### (see www.drive5.com/usearch for detailed documentation)

DATA_PATH=~/projects/thesis/data

# First truncate low-quality base calls at the beginning and end of the reads
for fq in ls $DATA_PATH/raw/*.fastq
do
    bn=$(basename $fq .fastq)
    usearch -fastq_filter $DATA_PATH/$bn -fastq_trunclen 220 -fastqout $DATA_PATH/clean/$bn".tmp"
    usearch -fastq_filter $DATA_PATH/clean/$bn".tmp" -fastq_stripleft 20 -fastqout $DATA_PATH/clean/$bn"_trunc.fastq"
done
	  
# Then merge forward and reverse reads from each sample, and pool the merged reads
usearch -fastq_mergepairs $DATA_PATH/clean/*R1_trunc.fastq \
	-fastqout $DATA_PATH/clean/pooled_merged.fastq \
	-fastq_maxdiffs 100 -fastq_pctid 50 -fastq_minmergelen 245 -fastq_maxmergelen 255 \
	-relabel @ -report $DATA_PATH/reports/merge_report.txt

# Create a report of the expected errors of the merged reads
usearch -fastq_eestats2 $DATA_PATH/clean/pooled_merged.fastq \
	-output $DATA_PATH/reports/merged_eestats.txt \
	-length_cutoffs 200,*,10 -ee_cutoffs 1.0,2.0,2.5,3.0,3.5,4.0

# Quality filter the reads using a maximum of 2.0 expected errors, after trimming nucleotides from the
# beginning and end of the reads
usearch -fastq_filter $DATA_PATH/clean/pooled_merged.fastq -fastqout $DATA_PATH/clean/pooled_filtered.fastq \
	-fastq_maxee 2.0

# Dereplicate the reads
usearch -fastx_uniques $DATA_PATH/clean/all_filtered.fastq -fastqout $DATA_PATH/clean/pooled_uniques.fastq \
	-sizeout -relabel Uniq

# Now cluster the reads into OTUs using the UPARSE algorithm
usearch -cluster_otus $DATA_PATH/clean/pooled_uniques.fastq -otus $DATA_PATH/uparse_otus/otus.fa \
	-relabel OTU
