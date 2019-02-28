#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

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
PCTID=80
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
#printf "\nMothur groups file: %s/%s\n\n" "$PWD" "$GROUPFILE"

# Create the 'truncated' directory, if necessary
if [ ! -d $WDIR/truncated ]; then
    mkdir $WDIR/truncated
fi

# Create the 'merged' directory, if necessary
if [ ! -d $WDIR/merged/pooled ]; then
    mkdir -p $WDIR/merged/pooled
else
    # remove fastx files from previous run
    rm -f $WDIR/merged/*.fast*
    rm -f $WDIR/merged/pooled/*.fast*
fi

# Create the 'filtered' directory, if necessary
if [ ! -d $WDIR/filtered/pooled ]; then
    mkdir -p $WDIR/filtered/pooled
else
    # remove fastx files from previous run
    rm -f $WDIR/filtered/pooled/*.fast*
    rm -f $WDIR/filtered/*.fast*
fi

# Create the 'reports' directory, if necessary
if [ ! -d $WDIR/reports ]; then
    mkdir $WDIR/reports
else
    # remove reports from previous run
    rm -f $WDIR/reports/*.txt
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
for fq in $WDIR/truncated/*_R1.fastq
do
    bn=$(basename $fq)
    sn=${bn%%_*.fastq}
    nn=$sn"_merged.fastq"
    usearch -fastq_mergepairs $fq \
	    -fastqout $WDIR/merged/$nn \
	    -fastq_maxdiffs $MAXDIFFS \
	    -fastq_minmergelen $MINMERGELEN \
	    -fastq_maxmergelen $MAXMERGELEN \
	    -report $WDIR/reports/$sn"_merge_report.txt" \
	    -relabel @

    # Generate merge report
    usearch -fastx_info $WDIR/merged/$nn \
	    -output $WDIR/reports/$sn"_merged_info.txt"

    # Add sequences to pooled file
    cat $WDIR/merged/$nn >> $WDIR/merged/pooled/pooled_merged.fastq

done


# Now filter the individual samples
for fq in $WDIR/merged/*_merged.fastq
do
    sn=$(basename $fq _merged.fastq)
    fasta=$sn"_filtered.fasta"
    fastq=$sn"_filtered.fastq"
    usearch -fastq_filter $fq \
	    -fastaout $WDIR/filtered/$fasta \
	    -fastqout $WDIR/filtered/$fastq \
	    -fastq_maxee $MAXEE \
	    -fastq_maxns $MAXNS 

    # Generate filter report
    usearch -fastx_info $WDIR/filtered/$fastq \
	    -output $WDIR/reports/$sn"_filtered_info.txt"

    # Reformat sequence labels for QIIME compatibility
    sed -i '/^>/ s/\(.*\)\./\1_/' $WDIR/filtered/$fasta

    # Add sequences to pooled files
    cat $WDIR/filtered/$fasta >> $WDIR/filtered/pooled/pooled_filtered.fasta
    cat $WDIR/filtered/$fastq >> $WDIR/filtered/pooled/pooled_filtered.fastq
done


### Create a report with summary stats on the pooled merged reads
usearch -fastx_info $WDIR/merged/pooled/pooled_merged.fastq \
	-output $WDIR/reports/pooled_merged_info.txt

### Create a report of the expected errors of the pooled merged reads
usearch -fastq_eestats2 $WDIR/merged/pooled/pooled_merged.fastq \
	-output $WDIR/reports/pooled_merged_eestats.txt \
	-length_cutoffs 200,*,10

### Create a report with summary stats on the pooled filtered reads
usearch -fastx_info $WDIR/filtered/pooled/pooled_filtered.fastq \
	-output $WDIR/reports/pooled_filtered_info.txt


