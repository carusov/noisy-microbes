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
OUTDIR=results
REFDIR=~/thesis/references

# Set default sequence length and filtering parameters
MIN_LEN=221
MAX_LEN=225
FTRUNC=230
RTRUNC=210
MAXEE_F=2.5
MAXEE_R=2.5

# Set the default processing mode
MODE="pooled"

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
	    printf "\nUSAGE: run_all_pipelines.sh [-i --input input_directory] [options]\n"
	    printf "\nOptions: \t[default]"
	    printf "\n-o --output \t[./"$OUTDIR"] \t\t\t\toutput_directory"
	    printf "\n-r --ref \t["$REFDIR"] \treference_directory"
	    printf "\n-f --ftrunc \t["$FTRUNC"] \t\t\t\t\tforward_trunc_position"
	    printf "\n-b --rtrunc \t["$RTRUNC"] \t\t\t\t\treverse_trunc_position"
	    printf "\n-s --min_len \t["$MIN_LEN"] \t\t\t\t\tmin_merge_length"
	    printf "\n-l --max_len \t["$MAX_LEN"] \t\t\t\t\tmax_merge_length"
	    printf "\n-F --maxee_F \t["$MAXEE_F"] \t\t\t\t\tmax_ee_forward"
	    printf "\n-R --maxee_R \t["$MAXEE_R"] \t\t\t\t\tmax_ee_reverse"
	    printf "\n-m --mode \t["$MODE"] \t\t\t\tprocessing mode"
	    printf "\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done


# Create the output directory, if necessary
if [ ! -d "$OUTDIR" ]; then
    mkdir -p $OUTDIR
fi

# Convert any relative paths to absolute paths
INDIR=$(readlink -f "$INDIR")
OUTDIR=$(readlink -f "$OUTDIR")
REFDIR=$(readlink -f "$REFDIR")

printf "\nINPUT DIRECTORY: %s" "$INDIR"
printf "\nOUTPUT DIRECTORY: %s" "$OUTDIR"
printf "\nREFERENCE DIRECTORY: %s" "$REFDIR"
printf "\nMinimum merge length: %d" $MIN_LEN
printf "\nMaximum merge length: %d" $MAX_LEN
printf "\nMaximum expected errors (DADA2 forward): %0.2f" $MAXEE_F
printf "\nMaximum expected errors (DADA2 reverse): %0.2f" $MAXEE_R
printf "\nProcessing mode: %s" "$MODE"
printf "\n\n"


# If pooled mode is selected, run five pipelines on the pooled fasta file
if [[ $MODE == pooled ]]
then
    run_5_pipelines.sh -q "$INDIR"/filtered/pooled/pooled_filtered.fastq \
		       -a "$INDIR"/filtered/pooled/pooled_filtered.fasta \
		       -o "$OUTDIR" \
		       -g "$INDIR"/filtered/mothur.groups \
		       -r "$INDIR"/merged/pooled/pooled_merged.fastq \
		       -R "$REFDIR" \
		       -t $MIN_LEN

# If separate mode is selected, run five pipelines on each sample fasta separately
elif [[ $MODE == separate ]]
then
    for s in $(ls "$INDIR"/filtered/separate | egrep '^s[0-9]{3}.*\.fastq');
    do
	fname="${s%.fastq}"
	sname="${s%_*.fastq}"
	fastq="$s"
	fasta="$fname".fasta
	raw="$sname"_merged.fastq

	run_5_pipelines.sh -q "$INDIR"/filtered/separate/"$fastq" \
			   -a "$INDIR"/filtered/separate/"$fasta" \
			   -o "$OUTDIR"/"$sname" \
			   -g "$INDIR"/filtered/mothur.groups \
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
SCRIPTS=~/thesis/noisy-microbes/community-inference/scripts
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
