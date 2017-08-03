#!/bin/bash

### Author: Vincent Caruso
### Date: 7/19/2017
### This script takes truncated forward and reverse reads, merges them, filters them
### using the expected errors criterion, and pools them into a single file. It uses
### the USEARCH software to perform merging and filtering.
### (see www.drive5.com/usearch for detailed documentation)

DATA=~/projects/thesis/data

### Merge parameters
maxdiffs=30
pctid=50
minmergelen=230
maxmergelen=235

### Filter parameters
maxee=2.0
maxns=0

### Merge reads from individual samples to get stats for publication
for fq in $(ls $DATA/truncated/*R1.fastq)
do
    bn=$(basename $fq _trunc_R1.fastq)
    nn=$bn"_merged.fastq"
    usearch -fastq_mergepairs $fq \
	    -fastqout $DATA/merged/$nn \
	    -relabel @ \
	    -fastq_maxdiffs $maxdiffs \
	    -fastq_pctid $pctid \
	    -fastq_minmergelen $minmergelen \
	    -fastq_maxmergelen $maxmergelen \
	    -report $DATA/reports/$bn"_merge_report.txt"

    usearch -fastx_info $DATA/merged/$nn \
	    -output $DATA/reports/$bn"_merged_info.txt"
done

### Now filter the individual samples too, in case we need stats on the output of filtering
# for fq in $(ls $DATA/clean/ | grep -E "s[0-9]{3}.+_merged.fastq")
for fq in $(ls $DATA/merged/s1*_merged.fastq)   # I was editing this line when I closed out last time
do
    bn=$(basename $fq _merged.fastq)
    nn=$bn"_filtered.fastq"
    usearch -fastq_filter $fq \
	    -fastqout $DATA/clean/$nn \
	    -fastq_maxee $maxee \
	    -fastq_maxns $maxns
    
    usearch -fastx_info $DATA/clean/$nn \
	    -output $DATA/reports/$bn"_filtered_info.txt"
done

################################################################################
### Now for the actual processing pipeline
################################################################################

### I might want to consider using the USEARCH low-complexity read filter, which filters low-complexity
### reads that are potentially the result of seqeuncing past the end of the reverse adapter.
# usearch -filter_lowc reads_R1.fastq \
#         -reverse reads_R2.fastq \
#         -output filtered_R1.fastq \
#         -output2 filtered_R2.fastq \
#         -tabbedout lowc.txt \
#         -hitsout lowc.fastq

### Also, was there any PhiX spike-in used in sequencing that might not have been completely removed?
### Answer: Any remaining PhiX reads are removed by 'fastqPairedFilter' in DADA2,
### at the same time that the reads are trimmed.

### And were the primer sequences removed before we received the "raw" data? They don't seem to be
### present in the reads.

### First merge forward and reverse reads from each sample, and pool the merged reads
usearch -fastq_mergepairs $DATA/truncated/*R1.fastq \
	-fastqout $DATA/merged/pooled_merged.fastq \
	-relabel @ \
	-fastq_maxdiffs $maxdiffs \
	-fastq_pctid $pctid \
	-fastq_minmergelen $minmergelen \
	-fastq_maxmergelen $maxmergelen \
        -report $DATA/reports/pooled_merge_report.txt

### Create a report with summary stats on the merged reads
usearch -fastx_info $DATA/merged/pooled_merged.fastq \
	-output $DATA/reports/pooled_merged_info.txt

### Create a report of the expected errors of the merged reads
usearch -fastq_eestats2 $DATA/merged/pooled_merged.fastq \
	-output $DATA/reports/pooled_merged_eestats.txt \
	-length_cutoffs 200,*,10

### Quality filter the reads using a maximum of 2.0 expected errors
usearch -fastq_filter $DATA/merged/pooled_merged.fastq \
	-fastqout $DATA/clean/pooled_filtered.fastq \
	-fastaout $DATA/clean/pooled_filtered.fasta \
	-fastq_maxee $maxee \
	-fastq_maxns $maxns

### Create a report with summary stats on the filtered reads
usearch -fastx_info $DATA/clean/pooled_filtered.fastq \
	-output $DATA/reports/pooled_filtered_info.txt
