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

# Set the default input, output, and reference directories
INDIR=""
OUTDIR=results
REFDIR=~/thesis/references

# Set default sequence length and filtering parameters
FTRUNC=230
RTRUNC=210
MIN_LEN=221
MAX_LEN=225
MAXEE_F=2.5
MAXEE_R=2.5
POOLED=false

# Set the default processing mode
#MODE="pooled"

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
	-F|--maxee_F)
	    MAXEE_F="$2"
	    shift;;
	-R|--maxee_R)
	    MAXEE_R="$2"
	    shift;;
	-p|--pooled)
	    POOLED=true;;
	-h|--help)
	    printf "\nUSAGE: run_all_pipelines.sh [-i --input input_directory] [options]\n"
	    printf "\nOptions: \t[default]"
	    printf "\n-o --output \t[./"$OUTDIR"] \t\t\t\toutput_directory"
	    printf "\n-r --ref \t["$REFDIR"] \treference_directory"
	    printf "\n-s --min_len \t["$MIN_LEN"] \t\t\t\t\tmin_merge_length"
	    printf "\n-l --max_len \t["$MAX_LEN"] \t\t\t\t\tmax_merge_length"
	    printf "\n-F --maxee_F \t["$MAXEE_F"] \t\t\t\t\tmax_ee_forward"
	    printf "\n-R --maxee_R \t["$MAXEE_R"] \t\t\t\t\tmax_ee_reverse"
	    printf "\n-p --pooled \t["$POOLED"] \t\t\t\tpooled_samples"
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
printf "\nPooled samples: %s" "$POOLED"
printf "\n\n"

# First, make sure we have the latest version of the DADA2 script
# Any updates are made to the Rmd file, which then must be knit to an R file
SCRIPTS=~/thesis/noisy-microbes/community-inference/scripts
pushd $SCRIPTS
rmd2r.R -i dada2_pipeline.Rmd
popd

# If pooled mode is selected, run 5 pipelines on the pooled fasta file
if $POOLED
then
    run_5_pipelines.sh -q "$INDIR"/filtered/pooled/pooled_filtered.fastq \
		       -a "$INDIR"/filtered/pooled/pooled_filtered.fasta \
		       -o "$OUTDIR"/ \
		       -r "$INDIR"/merged/pooled/pooled_merged.fastq \
		       -R "$REFDIR" \
		       -t $MIN_LEN

    # Run the DADA2 pipeline with pooled samples
    printf "\n######################################################################\n"
    printf "\nRunning the DADA2 pipeline...\n"
    printf "\n######################################################################\n"
    Rscript "$SCRIPTS"/dada2_pipeline.R -i "$INDIR" -o "$OUTDIR"/dada2 \
	    -s $MIN_LEN -l $MAX_LEN \
	    -F $MAXEE_F -R $MAXEE_R \
	    -p

# If pooled mode is not selected, run five pipelines on each sample fasta separately
else
    if [[ $(ls "$INDIR"/filtered/*.fastq | wc -l) > 1 ]]
    then
	multi=true
    else
	multi=false
    fi
    
    for s in "$INDIR"/filtered/*.fastq
    do
	fastq=$(basename "$s")
	fasta="${fastq/.fastq/.fasta}"
    	sname="${fastq%_*}"
	raw="$sname"_merged.fastq

	if $multi
	then
	    sample_outdir="$OUTDIR"/"$sname"
	else
	    sample_outdir="$OUTDIR"
	fi

	run_5_pipelines.sh -q "$INDIR"/filtered/"$fastq" \
			   -a "$INDIR"/filtered/"$fasta" \
			   -o "$sample_outdir" \
			   -r "$INDIR"/merged/"$raw" \
			   -R "$REFDIR" \
			   -t $MIN_LEN

	echo
	echo "Done with sample "${sname}""
	echo
    done

    # Run the DADA2 pipeline with defaults (samples processed separately)
    printf "\n######################################################################\n"
    printf "\nRunning the DADA2 pipeline...\n"
    printf "\n######################################################################\n"
    Rscript "$SCRIPTS"/dada2_pipeline.R -i "$INDIR" -o "$OUTDIR"/dada2 \
	    -s $MIN_LEN -l $MAX_LEN \
	    -F $MAXEE_F -R $MAXEE_R \
	    
fi
