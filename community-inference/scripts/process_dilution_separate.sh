#!/bin/bash

### This script is mainly meant to serve as a record for how each
### dataset was processed, including the final parameters. It can
### also be run to re-process all datasets if desired.

# Define some constants and default parameters
DATA=~/thesis/data
RESULTS=~/thesis/results
REFDIR=~/thesis/references
MAXDIFFS=10
#PCTID=90
MAXEE_F=2.0
MAXEE_R=2.0
MAXN=0

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
	-t|--trim_len)
	    TRIM_LEN="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: process_dilution_separate.sh [-i input_directory]\n"
	    printf "\t\t\t [-o output_directory] [-r reference_directory]\n"
	    printf "\t\t\t [-t trim_length]\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done

# Create the output directory if necessary
if [ ! -d "$OUTDIR" ];
then
    mkdir "$OUTDIR"
fi

# Convert any relative paths to absolute paths
INDIR=$(readlink -f "$INDIR")
OUTDIR=$(readlink -f "$OUTDIR")
REFDIR=$(readlink -f "$REFDIR")


printf "\nINPUT DIRECTORY = "${INDIR}""
printf "\nOUTPUT DIRECTORY = "${OUTDIR}""
printf "\nREFERENCE DIRECTORY = "${REFDIR}""
printf "\nTrim length: %d" $TRIM_LEN


# Process dilution series samples individually
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing Zymo dilution series samples individually\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"

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
		       -t $TRIM_LEN

    echo
    echo "Done with sample "${sname}""
    echo
done
