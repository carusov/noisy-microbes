#!/bin/bash

### Author: Vincent Caruso
### Date: 7/19/2017
### This script takes truncated forward and reverse reads, merges them, filters them
### using the expected errors criterion, and pools them into a single file. It uses
### the USEARCH software to perform merging and filtering.
### (see www.drive5.com/usearch for detailed documentation)

# Define the scripts path
SCRIPTS=~/projects/thesis/noisy-microbes/scripts

# Set the default working directory
WDIR=~/projects/thesis/data/dilution_w_blank

# Set the default truncation parameters
FTRUNC=230
RTRUNC=210

# Set the default merge parameters
MAXDIFFS=30
PCTID=50
MINMERGELEN=220
MAXMERGELEN=225

# Set the default filter parameters
MAXEE=2.0
MAXNS=0

# Parse command-line options
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	-w|--working_dir)
	    WDIR="$2"
	    shift;;
	-f|--ftrunc)
	    FTRUNC="$2"
	    shift;;
	-b|--rtrunc)
	    RTRUNC="$2"
	    shift;;
	-d|--maxdiffs)
	    MAXDIFFS="$2"
	    shift;;
	-p|--pctid)
	    PCTID="$2"
	    shift;;
	-s|--shortest)
	    MINMERGELEN="$2"
	    shift;;
	-l|--longest)
	    MAXMERGELEN="$2"
	    shift;;
	-e|--maxee)
	    MAXEE="$2"
	    shift;;
	-n|--maxns)
	    MAXNS="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: merge_and_filter [-w working_directory]\n"
	    printf "\t\t\t [-f fwd_trunc_pos] [-b rev_trunc_pos]\n"
	    printf "\t\t\t [-d max_merge_differences] [-p min_merge_pct_id]\n"
	    printf "\t\t\t [-s min_merge_length] [-l max_merge_length]\n"
	    printf "\t\t\t [-e max_expected_errors] [-n max_Ns]\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done

printf "\nWORKING DIRECTORY = ""${WDIR}""\n"

# Create the output directory, if necessary
#if [ ! -d $WDIR ]; then
#    mkdir $WDIR
#fi

# Create the 'truncated' directory, if necessary
if [ ! -d $WDIR/truncated ]; then
    mkdir $WDIR/truncated
else
    rm $WDIR/truncated/*
fi

# Create the 'merged' directory, if necessary
if [ ! -d $WDIR/merged ]; then
    mkdir $WDIR/merged
else
    rm $WDIR/merged/*
fi

# Create the 'filtered' directory, if necessary
if [ ! -d $WDIR/filtered ]; then
    mkdir $WDIR/filtered
else
    rm $WDIR/filtered/*
fi

# Create the 'reports' directory, if necessary
if [ ! -d $WDIR/reports ]; then
    mkdir $WDIR/reports
else
    rm $WDIR/reports/*
fi

# First, run the .Rmd script that truncates the .fastq reads
# Make sure we have the latest version of the script
cd $SCRIPTS
rmd2r.R -i fastq_truncate.Rmd
Rscript $SCRIPTS/fastq_truncate.R -d $WDIR -f $FTRUNC -b $RTRUNC

### Merge reads from individual samples to get stats for publication
for fq in $(ls $WDIR/truncated/*_R1.fastq)
do
    bn=$(basename $fq)
    bn=${bn%%_*.fastq}
    nn=$bn"_merged.fastq"
    usearch -fastq_mergepairs $fq \
	    -fastqout $WDIR/merged/$nn \
	    -relabel @ \
	    -fastq_maxdiffs $MAXDIFFS \
	    -fastq_pctid $PCTID \
	    -fastq_minmergelen $MINMERGELEN \
	    -fastq_maxmergelen $MAXMERGELEN \
	    -report $WDIR/reports/$bn"_merge_report.txt"

    usearch -fastx_info $WDIR/merged/$nn \
	    -output $WDIR/reports/$bn"_merged_info.txt"
done

# Now filter the individual samples too, in case we need stats on the output of filtering
for fq in $(ls $WDIR/merged/*_merged.fastq)  
do
    bn=$(basename $fq _merged.fastq)
    nn=$bn"_filtered.fastq"
    usearch -fastq_filter $fq \
	    -fastqout $WDIR/filtered/$nn \
	    -fastq_maxee $MAXEE \
	    -fastq_maxns $MAXNS
    
    usearch -fastx_info $WDIR/filtered/$nn \
	    -output $WDIR/reports/$bn"_filtered_info.txt"
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
usearch -fastq_mergepairs $WDIR/truncated/*R1.fastq \
	-fastqout $WDIR/merged/pooled_merged.fastq \
	-relabel @ \
	-fastq_maxdiffs $MAXDIFFS \
	-fastq_pctid $PCTID \
	-fastq_minmergelen $MINMERGELEN \
	-fastq_maxmergelen $MAXMERGELEN \
        -report $WDIR/reports/pooled_merge_report.txt

### Create a report with summary stats on the merged reads
usearch -fastx_info $WDIR/merged/pooled_merged.fastq \
	-output $WDIR/reports/pooled_merged_info.txt

### Create a report of the expected errors of the merged reads
usearch -fastq_eestats2 $WDIR/merged/pooled_merged.fastq \
	-output $WDIR/reports/pooled_merged_eestats.txt \
	-length_cutoffs 200,*,10

### Quality filter the reads using a maximum of 2.0 expected errors
usearch -fastq_filter $WDIR/merged/pooled_merged.fastq \
	-fastqout $WDIR/filtered/pooled_filtered.fastq \
	-fastaout $WDIR/filtered/pooled_filtered.fasta \
	-fastq_maxee $MAXEE \
	-fastq_maxns $MAXNS

### Create a report with summary stats on the filtered reads
usearch -fastx_info $WDIR/filtered/pooled_filtered.fastq \
	-output $WDIR/reports/pooled_filtered_info.txt

# Reformat sequence files to QIIME format
printf "\nReformatting read sequences to QIIME format...\n"
sed '/^>/ s/\./_/' $WDIR/filtered/pooled_filtered.fasta \
    > $WDIR/filtered/pooled_filtered_qiime.fasta

