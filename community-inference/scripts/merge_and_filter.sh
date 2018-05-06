#!/bin/bash

### Author: Vincent Caruso
### Date: 7/19/2017
### This script takes truncated forward and reverse reads, merges them, filters them
### using the expected errors criterion, and pools them into a single file. It uses
### the USEARCH software to perform merging and filtering.
### (see www.drive5.com/usearch for detailed documentation)

# Define the scripts path
SCRIPTS=~/thesis/noisy-microbes/community-inference/scripts

# Set the default working directory
WDIR=$PWD

# Set the default mothur groups file
GROUPFILE="filtered/mothur.groups"

# Set the default truncation parameters
FTRUNC=230
RTRUNC=210

# Set the default merge parameters
MAXDIFFS=10
#PCTID=80
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
	-g|--groups)
	    GROUPFILE="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: merge_and_filter [-w working_directory]\n"
	    printf "\t\t\t [-f fwd_trunc_pos] [-b rev_trunc_pos]\n"
	    printf "\t\t\t [-d max_merge_differences] [-p min_merge_pct_id]\n"
	    printf "\t\t\t [-s min_merge_length] [-l max_merge_length]\n"
	    printf "\t\t\t [-e max_expected_errors] [-n max_Ns]\n"
	    printf "\t\t\t [-g mothur_groups_file]\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done

# Get full path in case relative path is given by user
WDIR=$(readlink -f $WDIR)
#GROUPFILE=$(readlink -f $GROUPFILE)

printf "\nWORKING DIRECTORY: %s" "$WDIR"
printf "\n\nTRUNCATION parameters:"
printf "\nForward read truncate position: %d" $FTRUNC
printf "\nReverse read truncate position: %d" $RTRUNC
printf "\n\nMERGE parameters:"
printf "\nMaximum differences: %d" $MAXDIFFS
printf "\nMinimum merge length: %d" $MINMERGELEN
printf "\nMaximum merge length: %d" $MAXMERGELEN
printf "\n\nFILTER parameters:"
printf "\nMaximum expected errors: %0.2f" $MAXEE
printf "\nMaximum Ns: %d" $MAXNS
printf "\nMothur groups file: %s/%s\n\n" "$PWD" "$GROUPFILE"

# Create the 'truncated' directory, if necessary
if [ ! -d $WDIR/truncated ]; then
    mkdir $WDIR/truncated
fi

# Create the 'merged' directory, if necessary
if [ ! -d $WDIR/merged ]; then
    mkdir -p $WDIR/merged/pooled
else
    rm $WDIR/merged/*
    rm $WDIR/merged/pooled/*
fi

# Create the 'filtered' directory, if necessary
if [ ! -d $WDIR/filtered ]; then
    mkdir -p $WDIR/filtered/pooled
    mkdir -p $WDIR/filtered/separate
else
    rm $WDIR/filtered/pooled/*
    rm $WDIR/filtered/separate/*
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
if [ ! "$(ls -A "$WDIR"/truncated)" ];then
    pushd $SCRIPTS
    rmd2r.R -i fastq_truncate.Rmd
    Rscript $SCRIPTS/fastq_truncate.R -d $WDIR -f $FTRUNC -b $RTRUNC
    popd
fi

# Merge reads from individual samples
for fq in $(ls $WDIR/truncated/*_R1.fastq)
do
    bn=$(basename $fq)
    bn=${bn%%_*.fastq}
    nn=$bn"_merged.fastq"
    usearch -fastq_mergepairs $fq \
	    -fastqout $WDIR/merged/$nn \
	    -fastq_maxdiffs $MAXDIFFS \
	    -fastq_minmergelen $MINMERGELEN \
	    -fastq_maxmergelen $MAXMERGELEN \
	    -report $WDIR/reports/$bn"_merge_report.txt" \
	    -relabel @

    # Reformat sequence labels for QIIME compatibility
    sed -i '/^@/ s/\(.*\)\./\1_/' $WDIR/merged/$nn

    # Generate merge report
    usearch -fastx_info $WDIR/merged/$nn \
	    -output $WDIR/reports/$bn"_merged_info.txt"

    # Add sequences to pooled file
    cat $WDIR/merged/$nn >> $WDIR/merged/pooled/pooled_merged.fastq

done


# Now filter the individual samples
for fq in $(ls $WDIR/merged/*_merged.fastq)  
do
    bn=$(basename $fq _merged.fastq)
    fasta=$bn"_filtered.fasta"
    fastq=$bn"_filtered.fastq"
    usearch -fastq_filter $fq \
	    -fastaout $WDIR/filtered/separate/$fasta \
	    -fastqout $WDIR/filtered/separate/$fastq \
	    -fastq_maxee $MAXEE \
	    -fastq_maxns $MAXNS 

    # Generate filter report
    usearch -fastx_info $WDIR/filtered/separate/$fastq \
	    -output $WDIR/reports/$bn"_filtered_info.txt"

    # Add sequences to pooled files
    cat $WDIR/filtered/separate/$fasta >> $WDIR/filtered/pooled/pooled_filtered.fasta
    cat $WDIR/filtered/separate/$fastq >> $WDIR/filtered/pooled/pooled_filtered.fastq
    
    # Add sequence labels to mothur groups file
    make_mothur_groups.awk $WDIR/filtered/separate/$fasta >> $GROUPFILE

    # Change sequence labels to sample name and a number
#    printf "\nRelabelling sequences...\n"
#    usearch -fastx_relabel $WDIR/filtered/$nn -prefix $bn"_" -fastaout $WDIR/filtered/$nn".temp"
#    sed -i '/^>/ s/:/_/g' $WDIR/filtered/separate/$nn
#    mv $WDIR/filtered/$nn".temp" $WDIR/filtered/$nn

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
#usearch -fastq_mergepairs $WDIR/truncated/*R1.fastq \
#	-fastqout $WDIR/merged/pooled_merged.fastq \
#	-relabel @ \
#	-fastq_maxdiffs $MAXDIFFS \
#	-fastq_minmergelen $MINMERGELEN \
#	-fastq_maxmergelen $MAXMERGELEN \
 #       -report $WDIR/reports/pooled_merge_report.txt

# Reformat sequence labels for QIIME compatibility
#sed -i '/^@/ s/\(.*\)\./\1_/' $WDIR/merged/pooled_merged.fastq

### Create a report with summary stats on the merged reads
usearch -fastx_info $WDIR/merged/pooled/pooled_merged.fastq \
	-output $WDIR/reports/pooled_merged_info.txt

### Create a report of the expected errors of the merged reads
usearch -fastq_eestats2 $WDIR/merged/pooled/pooled_merged.fastq \
	-output $WDIR/reports/pooled_merged_eestats.txt \
	-length_cutoffs 200,*,10

### Quality filter the reads using a maximum of 2.0 expected errors
#usearch -fastq_filter $WDIR/merged/pooled_merged.fastq \
#	-fastqout $WDIR/filtered/pooled/pooled_filtered.fastq \
#	-fastaout $WDIR/filtered/pooled/pooled_filtered.fasta \
#	-fastq_maxee $MAXEE \
#	-fastq_maxns $MAXNS

### Create a report with summary stats on the filtered reads
usearch -fastx_info $WDIR/filtered/pooled/pooled_filtered.fastq \
	-output $WDIR/reports/pooled_filtered_info.txt

# Reformat sequence files to QIIME format
#printf "\nReformatting read sequences to QIIME format...\n"
#sed '/^>/ s/\(.*\)\./\1_/' $WDIR/filtered/pooled/pooled_filtered.fasta \
#    > $WDIR/filtered/pooled/pooled_filtered_qiime.fasta

