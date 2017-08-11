#!/bin/bash

### Author: Vincent Caruso
### Date: 7/19/2017
### This script takes truncated forward and reverse reads, merges them, filters them
### using the expected errors criterion, and pools them into a single file. It uses
### the USEARCH software to perform merging and filtering.
### (see www.drive5.com/usearch for detailed documentation)

# Set the default input and output directories
#DATA=~/projects/thesis/data
INDIR=./
OUTDIR=./

# Set the default merge parameters
MAXDIFFS=30
PCTID=50
MINMERGELEN=230
MAXMERGELEN=250

# Set the default filter parameters
MAXEE=2.0
MAXNS=0

# Parse command-line options
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	-i|--input)
	    INDIR="$2"
	    shift;;
	-o|--output)
	    OUTDIR="$2"
	    shift;;
	-m|--maxdiffs)
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
	    printf "USAGE: merge_and_filter -i input_directory -o output_directory [OPTIONS]\n"
	    exit;;
	*)

	;;
    esac
    shift
done

echo INPUT DIRECTORY = "${INDIR}"
echo OUTPUT DIRECTORY = "${OUTDIR}"

# Create the 'merged' directory, if necessary
if [ ! -e $OUTDIR/merged ]; then
    mkdir $OUTDIR/merged
fi

# Create the 'filtered' directory, if necessary
if [ ! -e $OUTDIR/filtered ]; then
    mkdir $OUTDIR/filtered
fi

# Create the 'reports' directory, if necessary
if [ ! -e $OUTDIR/reports ]; then
    mkdir $OUTDIR/reports
fi

### Merge reads from individual samples to get stats for publication
for fq in $(ls $INDIR/*R1.fastq)
do
    bn=$(basename $fq)
    printf "********************\n"
    printf $bn"\n"
    bn=${bn%%_*.fastq}
    printf $bn"\n"
    nn=$bn"_merged.fastq"
    usearch -fastq_mergepairs $fq \
	    -fastqout $OUTDIR/merged/$nn \
	    -relabel @ \
	    -fastq_maxdiffs $MAXDIFFS \
	    -fastq_pctid $PCTID \
	    -fastq_minmergelen $MINMERGELEN \
	    -fastq_maxmergelen $MAXMERGELEN \
	    -report $OUTDIR/reports/$bn"_merge_report.txt"

    usearch -fastx_info $OUTDIR/merged/$nn \
	    -output $OUTDIR/reports/$bn"_merged_info.txt"
done

# Now filter the individual samples too, in case we need stats on the output of filtering
for fq in $(ls $OUTDIR/merged/s1*_merged.fastq)  
do
    bn=$(basename $fq _merged.fastq)
    nn=$bn"_filtered.fastq"
    usearch -fastq_filter $fq \
	    -fastqout $OUTDIR/filtered/$nn \
	    -fastq_maxee $MAXEE \
	    -fastq_maxns $MAXNS
    
    usearch -fastx_info $OUTDIR/filtered/$nn \
	    -output $OUTDIR/reports/$bn"_filtered_info.txt"
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
usearch -fastq_mergepairs $INDIR/*R1.fastq \
	-fastqout $OUTDIR/merged/pooled_merged.fastq \
	-relabel @ \
	-fastq_maxdiffs $MAXDIFFS \
	-fastq_pctid $PCTID \
	-fastq_minmergelen $MINMERGELEN \
	-fastq_maxmergelen $MAXMERGELEN \
        -report $OUTDIR/reports/pooled_merge_report.txt

### Create a report with summary stats on the merged reads
usearch -fastx_info $OUTDIR/merged/pooled_merged.fastq \
	-output $OUTDIR/reports/pooled_merged_info.txt

### Create a report of the expected errors of the merged reads
usearch -fastq_eestats2 $OUTDIR/merged/pooled_merged.fastq \
	-output $OUTDIR/reports/pooled_merged_eestats.txt \
	-length_cutoffs 200,*,10

### Quality filter the reads using a maximum of 2.0 expected errors
usearch -fastq_filter $OUTDIR/merged/pooled_merged.fastq \
	-fastqout $OUTDIR/filtered/pooled_filtered.fastq \
	-fastaout $OUTDIR/filtered/pooled_filtered.fasta \
	-fastq_maxee $MAXEE \
	-fastq_maxns $MAXNS

### Create a report with summary stats on the filtered reads
usearch -fastx_info $OUTDIR/filtered/pooled_filtered.fastq \
	-output $OUTDIR/reports/pooled_filtered_info.txt
