#!/bin/bash

### Author: Vincent Caruso
### Date: 8/14/2017
### Purpose: This is a wrapper shell script that runs each of the
### 16S processing method scripts in sequence. This script specifies
### paths for the inputs and outputs to each script, and is thus
### specific to my configuration, so it is not convenient for use
### on someone else's machine at the moment.
### Usage: run_all_pipelines.sh

# Set the default input, output, and reference directories
INDIR=~/projects/thesis/data/dilution_w_blank
OUTDIR=~/projects/thesis/results/dilution_w_blank
REFDIR=~/projects/thesis/references

# Set default 
MIN_LEN=230
MAX_LEN=235

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
	-r|--ref)
	    REFDIR="$2"
	    shift;;
	-s|--min_len)
	    MIN_LEN="$2"
	    shift;;
	-l|--max_len)
	    MAX_LEN="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: run_all_pipelines.sh -i input_directory -o output_directory -r reference_directory\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done

printf "\nINPUT DIRECTORY = ""${INDIR}""\n"
printf "OUTPUT DIRECTORY = ""${OUTDIR}""\n"
printf "REFERENCE DIRECTORY = ""${REFDIR}""\n\n"

# Create the output directory, if necessary
if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR
fi

# Run the UCLUST pipeline
printf "\nRunning the UCLUST pipeline...\n\n"
uclust_pipeline.sh -i $INDIR/filtered/pooled_filtered_qiime.fasta \
		   -o $OUTDIR/uclust \
		   -r $REFDIR/gold.fa

# Run the UPARSE pipeline
printf "\nRunning the UPARSE pipeline...\n\n"
uparse_pipeline.sh -i $INDIR/filtered/pooled_filtered.fastq \
		   -o $OUTDIR/uparse \
		   -r $INDIR/merged/pooled_merged.fastq

# Run the UNOISE pipeline
printf "\nRunning the UNOISE pipeline...\n\n"
unoise_pipeline.sh -i $INDIR/filtered/pooled_filtered.fastq \
		   -o $OUTDIR/unoise \
		   -r $INDIR/merged/pooled_merged.fastq

# Run the MED pipeline
printf "\nRunning the MED pipeline...\n\n"
med_pipeline.sh -i $INDIR/filtered/pooled_filtered_qiime.fasta \
		-o $OUTDIR/med

# Run the Deblur pipeline
printf "\nRunning the Deblur pipeline...\n\n"
deblur_pipeline.sh -i $INDIR/filtered/pooled_filtered_qiime.fasta \
		   -o $OUTDIR/deblur

# Run the DADA2 pipeline with defaults
printf "\nRunning the DADA2 pipeline...\n\n"
SCRIPTS=~/projects/thesis/noisy-microbes/scripts
#$SCRIPTS/rmd2r.R -i $SCRIPTS/dada2_pipeline.Rmd -d $SCRIPTS
Rscript $SCRIPTS/dada2_pipeline.R -i $INDIR -o $OUTDIR \
	-s $MIN_LEN -l $MAX_LEN
