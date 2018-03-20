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
#INDIR=~/thesis/data/dilution
#OUTDIR=~/thesis/results/dilution
REFDIR=~/thesis/references

# Set default sequence length and filtering parameters
MIN_LEN=220
MAX_LEN=225
FTRUNC=230
RTRUNC=210

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
	-f|--ftrunc)
	    FTRUNC="$2"
	    shift;;
	-b|--rtrunc)
	    RTRUNC="$2"
	    shift;;
	-s|--min_len)
	    MIN_LEN="$2"
	    shift;;
	-l|--max_len)
	    MAX_LEN="$2"
	    shift;;
	-F|--maxee_F)
	    MAXEE_F="$2"
	    shift;;
	-R|--maxee_R)
	    MAXEE_R="$2"
	    shift;;
	-m|--mode)
	    MODE="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: run_all_pipelines.sh [-i input_directory]\n"
	    printf "\t\t\t [-o output_directory] [-r reference_directory]\n"
	    printf "\t\t\t [-f fwd_trunc_pos] [-b rev_trunc_pos]\n"
	    printf "\t\t\t [-s min_merge_length] [-l max_merge_length]\n"
	    printf "\t\t\t [-m mode]\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done

# Convert any relative paths to absolute paths
INDIR=$(readlink -f "$INDIR")
OUTDIR=$(readlink -f "$OUTDIR")
REFDIR=$(readlink -f "$REFDIR")

printf "\nINPUT DIRECTORY = "${INDIR}""
printf "\nOUTPUT DIRECTORY = "${OUTDIR}""
printf "\nREFERENCE DIRECTORY = "${REFDIR}""
printf "\nMinimum merge length: %d" $MIN_LEN
printf "\nMaximum merge length: %d" $MAX_LEN
printf "\nMaximum expected errors (DADA2 forward): %0.2f" $MAXEE_F
printf "\nMaximum expected errors (DADA2 reverse): %0.2f" $MAXEE_R
printf "\nProcessing mode: "${MODE}""


# Create the output directory, if necessary
if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR
fi

# If pooled mode is selected, run five pipelines on the pooled fasta file
if [[ "$MODE" == pooled ]]; then
    run_5_pipelines.sh -q "$INDIR"/filtered/pooled_filtered.fastq \
		       -a "$INDIR"/filtered/pooled_filtered.fasta \
		       -o "$OUTDIR" \
		       -r "$INDIR"/merged/pooled_merged.fastq \
		       -R "$REFDIR" \
		       -t $MIN_LEN

# If separate mode is selected, run five pipelines on each sample fasta separately
elif [[ "$MODE" == separate ]]; then
    for s in $(ls "$INDIR"/filtered | egrep '^s[0-9]{3}.*\.fastq');
    do
	fname="${s%.fastq}"
	sname="${s%_*.fastq}"
	fastq="$s"
	fasta="$fname".fasta
	raw="$sname"_merged.fastq

	run_5_pipelines.sh -q "$INDIR"/filtered/"$fastq" \
			   -a "$INDIR"/filtered/"$fasta" \
			   -o "$OUTDIR"/"$sname" \
			   -r "$INDIR"/merged/"$raw" \
			   -R "$REFDIR" \
			   -t $MIN_LEN

	echo
	echo "Done with sample "${sname}""
	echo
    done
else
    echo "Invalid MODE parameter. Must be one of 'pooled' or 'separate'"
    echo
    exit 1
fi


# Run the DADA2 pipeline with defaults (DADA2 always processes samples separately)
# First, make sure we have the latest version of the script
SCRIPTS=~/thesis/noisy-microbes/scripts
pushd $SCRIPTS
rmd2r.R -i dada2_pipeline.Rmd
popd

# Now run it
printf "\n######################################################################\n"
printf "\nRunning the DADA2 pipeline...\n"
printf "\n######################################################################\n"
Rscript $SCRIPTS/dada2_pipeline.R -i $INDIR -o $OUTDIR/dada2 \
	-s $MIN_LEN -l $MAX_LEN \
	-F $MAXEE_F -R $MAXEE_R
