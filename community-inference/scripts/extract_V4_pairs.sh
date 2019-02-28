#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

### Author: Vincent Caruso
### Date: 9/16/2017
### Purpose: This script extracts the paired-end reads that come from the V4
### region of the 16S rRNA gene. The input is a pair of fasta files for the
### same sample from an Illumina paired-end sequencing run. The script works by
### first merging paired reads, assuming that reads from the V4 region will have
### a merge length between 250-255 base pairs. Only read pairs that merge within
### this range are retained.
### Usage: extract_V4_pairs.sh -f sample_R1.fastq -r sample_R2.fastq

# Define default merge parameters
MIN_LEN=250
MAX_LEN=255
MAX_DIFFS=30

# Get command-line options
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	-f|--fwd)
	    FWD="$2"
	    shift;;
	-r|--rev)
	    REV="$2"
	    shift;;
	-s|--min_len)
	    MIN_LEN="$2"
	    shift;;
	-l|--max_len)
	    MAX_LEN="$2"
	    shift;;
	-d|--max_diffs)
	    MAX_DIFFS="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: run_all_pipelines.sh -f fwd_read.fastq [-r rev_read.fastq]\n"
	    printf "\t\t\t [-s min_merge_length] [-l max_merge_length] [-d max_diffs]\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done

#fn=$(basename "$FWD")
sn=${FWD%_R1*.fastq}
ext=${FWD#*R1}
bn=$(basename $sn)

# Check to see if reverse read was specified
if [ -z "$REV" ];then
    REV="$sn"_R2"$ext"
    usearch -fastq_mergepairs "$FWD" \
	    -fastaout "$bn"_merged.fasta \
	    -fastq_minmergelen $MIN_LEN \
	    -fastq_maxmergelen $MAX_LEN \
	    -fastq_maxdiffs $MAX_DIFFS \
	    -tabbedout merge_tabout.txt \
	    -report merge_info.txt
else
    usearch -fastq_mergepairs "$FWD" \
	    -reverse "$REV" \
	    -fastaout "$bn"_merged.fasta \
	    -fastq_minmergelen $MIN_LEN \
	    -fastq_maxmergelen $MAX_LEN \
	    -fastq_maxdiffs $MAXDIFFS \
	    -tabbedout merge_tabout.txt \
	    -report merge_info.txt
fi

# Get only labels of successfully merged reads
grep result=merged merge_tabout.txt | cut -f1 > merged.labels
usearch -fastx_getseqs "$FWD" -labels merged.labels -trunclabels -fastqout "$bn"_V4_R1.fastq
usearch -fastx_getseqs "$REV" -labels merged.labels -trunclabels -fastqout "$bn"_V4_R2.fastq

# Remove intermediate files
rm "$bn"_merged.fasta
rm merge_tabout.txt
rm merged.labels
