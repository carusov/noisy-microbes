#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

### Author: Vincent Caruso
### Date: 7/27/2017
### Purpose: This script implements the MED pipeline
### (see www.merenlab.org for documentation)


# Set the default input file and output directory
INFILE=""
OUTDIR=med

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
	-h|--help)
	    printf "\nUSAGE: med_pipeline.sh -i input_file -o output_directory\n\n"
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

# Open the MED environment using miniconda
set +u
source activate med
set -u

# Pad shorter reads with gaps so all reads have the same length.
# (This is an MED requirement.)
o-pad-with-gaps $INFILE \
		-o $OUTDIR/pooled_filtered_med.fasta

# Finally, run the MED algorithm on the formatted fasta file   
decompose $OUTDIR/pooled_filtered_med.fasta -o $OUTDIR \
	  --skip-check-input-file

# Now check the MED nodes for chimeras using UCHIME
if [ ! -d "$OUTDIR"/uchime ];then
    mkdir "$OUTDIR"/uchime
fi

med2uchime.sh -i "$OUTDIR"/NODE-REPRESENTATIVES.fasta \
	      -o "$OUTDIR"/uchime/NODE-REPS-uchime.fasta

usearch -sortbysize "$OUTDIR"/uchime/NODE-REPS-uchime.fasta \
	-fastaout "$OUTDIR"/uchime/NODE-REPS-uchime.fasta

usearch9 -uchime2_denovo "$OUTDIR"/uchime/NODE-REPS-uchime.fasta \
	 -uchimeout "$OUTDIR"/uchime/uchime_tab.txt \
	 -chimeras "$OUTDIR"/uchime/NODE-chimeras.fasta \
	 -nonchimeras "$OUTDIR"/uchime/NODE-nonchimeras.fasta \
	 -abskew 16

# Deactivate the med environment
set +u
source deactivate
set -u
