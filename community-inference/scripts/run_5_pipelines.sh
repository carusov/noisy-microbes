#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

### Author: Vincent Caruso
### Date: 8/14/2017
### Purpose: This is a wrapper shell script that runs each of the
### 16S processing method scripts in sequence. This script specifies
### paths for the inputs and outputs to each script, and is thus
### specific to my configuration, so it is not convenient for use
### on someone else's machine at the moment.
### Usage: run_all_pipelines.sh

# Set the default reference directory
FASTQ=""
FASTA=""
OUTDIR=""
RAWFILE=""
REFDIR=~/thesis/references
TRIM_LEN=""


# Parse command-line options
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	-q|--fastq_in_file)
	    FASTQ="$2"
	    shift;;
	-a|--fasta_in_file)
	    FASTA="$2"
	    shift;;
	-o|--out_dir)
	    OUTDIR="$2"
	    shift;;
	-r|--raw_file)
	    RAWFILE="$2"
	    shift;;
	-R|--ref)
	    REFDIR="$2"
	    shift;;
	-t|--trim_len)
	    TRIM_LEN="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: run_5_pipelines.sh [-q fastq_input_file]\n"
	    printf "\t\t\t [-a fastq_input_file]\n"
	    printf "\t\t\t [-o output_directory] [-r raw_file]\n"
	    printf "\t\t\t [-R reference_directory] [-t trim_length]\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done

# Convert any relative paths to absolute paths
FASTQ=$(readlink -f "$FASTQ")
FASTA=$(readlink -f "$FASTA")
OUTDIR=$(readlink -f "$OUTDIR")
RAWFILE=$(readlink -f "$RAWFILE")
REFDIR=$(readlink -f "$REFDIR")

printf "\nFASTQ FILE = "${FASTQ}""
printf "\nFASTA FILE = "${FASTA}""
printf "\nRAW FILE = "${RAWFILE}""
printf "\nOUTPUT DIRECTORY = "${OUTDIR}""
printf "\nREFERENCE DIRECTORY = "${REFDIR}""
printf "\nDeblur trim length: %d\n" $TRIM_LEN


# Create the output directory, if necessary
if [ ! -d "$OUTDIR" ]; then
    mkdir "$OUTDIR"
fi

# Run the UCLUST pipeline
printf "\n######################################################################\n"
printf "\nRunning the UCLUST pipeline...\n"
printf "\n######################################################################\n"
uclust_pipeline.sh -i "$FASTA" \
		   -o "$OUTDIR"/uclust \
		   -r "$REFDIR"/gold.fa

# Run the UPARSE pipeline
printf "\n######################################################################\n"
printf "\nRunning the UPARSE pipeline...\n"
printf "\n######################################################################\n"
uparse_pipeline.sh -i "$FASTQ" \
		   -o "$OUTDIR"/uparse \
		   -r "$RAWFILE"

# Run the UNOISE pipeline
printf "\n######################################################################\n"
printf "\nRunning the UNOISE pipeline...\n"
printf "\n######################################################################\n"
unoise_pipeline.sh -i "$FASTQ" \
		   -o "$OUTDIR"/unoise \
		   -r "$RAWFILE"

# Run the MED pipeline
printf "\n######################################################################\n"
printf "\nRunning the MED pipeline...\n"
printf "\n######################################################################\n"
med_pipeline.sh -i "$FASTA" \
		-o "$OUTDIR"/med

# Run the Deblur pipeline
printf "\n######################################################################\n"
printf "\nRunning the Deblur pipeline...\n"
printf "\n######################################################################\n"
deblur_pipeline.sh -i "$FASTA" \
		   -o "$OUTDIR"/deblur \
		   -t $TRIM_LEN
