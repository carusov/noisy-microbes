#!/bin/bash

### Author: Vincent Caruso
### Date: 7/19/2017
### This script takes truncated forward and reverse reads, merges them, filters them
### using the expected errors criterion, and pools them into a single file. It uses
### the USEARCH software to perform merging and filtering.
### (see www.drive5.com/usearch for detailed documentation)

DATA=~/projects/thesis/data

### I might want to consider using the USEARCH low-complexity read filter, which filters low-complexity
### reads that are potentially the result of seqeuncing past the end of the reverse adapter.
# usearch -filter_lowc reads_R1.fastq \
#         -reverse reads_R2.fastq \
#         -output filtered_R1.fastq \
#         -output2 filtered_R2.fastq \
#         -tabbedout lowc.txt \
#         -hitsout lowc.fastq

### Also, was there any PhiX spike-in used in sequencing that might not have been completely removed?

### And were the primer sequences removed before we received the "raw" data? They don't seem to be
### present in the reads.

# First merge forward and reverse reads from each sample, and pool the merged reads
usearch -fastq_mergepairs $DATA/truncated/*R1.fastq \
	-fastqout $DATA/clean/pooled_merged.fastq \
	-relabel @ \
	-fastq_maxdiffs 100 \
	-fastq_pctid 50 \
	-fastq_minmergelen 220 \
	-fastq_maxmergelen 225 \
        -report $DATA/reports/merge_report.txt

# Create a report with summary stats on the merged reads
usearch -fastx_info $DATA/clean/pooled_merged.fastq \
	-output $DATA/reports/merged_info.txt

# Create a report of the expected errors of the merged reads
usearch -fastq_eestats2 $DATA/clean/pooled_merged.fastq \
	-output $DATA/reports/merged_eestats.txt \
	-length_cutoffs 200,*,10

# Quality filter the reads using a maximum of 2.0 expected errors, after trimming nucleotides from the
# beginning and end of the reads
usearch -fastq_filter $DATA/clean/pooled_merged.fastq \
	-fastqout $DATA/clean/pooled_filtered.fastq \
	-fastq_maxee 2.0

# Create a report with summary stats on the filtered reads
usearch -fastx_info $DATA/clean/pooled_filtered.fastq \
	-output $DATA/reports/filtered_info.txt
