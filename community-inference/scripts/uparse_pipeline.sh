#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

### Author: Vincent Caruso
### Date: 7/1/2017
### This script implements the UPARSE pipeline according to its recommended usage
### (see www.drive5.com/usearch for detailed documentation)
### The input to the script is a .fastq file of merged, pooled, and quality-
### filtered sample reads

# Set the default input file, output directory, and raw merged read file
INFILE=""
OUTDIR=uparse
RAW_MERGED_FILE=""

# Parse command-line options
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	-i|--input)
	    INFILE="$2"
	    shift;;
	-o|--output)
	    OUTDIR="$2"
	    shift;;
	-r|--raw-merged)
	    RAW_MERGED_FILE="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: uparse_pipeline.sh -i input_file -o output_directory -r raw_merged_reads_file\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done

printf "\nINPUT FILE = ""${INFILE}""\n"
printf "OUTPUT DIRECTORY = ""${OUTDIR}""\n\n"

# Create the output directory, if necessary
if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR
fi

# Dereplicate the quality-filtered reads
usearch -fastx_uniques $INFILE \
	-fastqout $OUTDIR/pooled_uniques.fastq \
	-sizeout \
	-relabel Uniq

# Now cluster the reads into OTUs using the UPARSE algorithm
usearch -cluster_otus $OUTDIR/pooled_uniques.fastq \
	-otus $OUTDIR/otus.fa \
	-relabel OTU

# Use the generated OTUs to create an OTU table
usearch -otutab $RAW_MERGED_FILE \
	-otus $OUTDIR/otus.fa \
	-otutabout $OUTDIR/otu_table.txt \
	-biomout $OUTDIR/otu_table.biom \
	-mapout $OUTDIR/map.txt \
	-dbmatched $OUTDIR/otus_with_sizes.fa \
	-notmatchedfq $OUTDIR/unmapped_reads.fastq \
	-sizeout
