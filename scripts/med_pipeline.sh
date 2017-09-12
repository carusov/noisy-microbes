#!/bin/bash

### Author: Vincent Caruso
### Date: 7/27/2017
### Purpose: This script implements the MED pipeline
### (see www.merenlab.org for documentation)


# Set the default input file and output directory
INFILE=~/projects/thesis/data/dilution_w_blank/filtered/pooled_filtered_qiime.fasta
OUTDIR=~/projects/thesis/results/dilution_w_blank/med

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
source activate med

# Pad shorter reads with gaps so all reads have the same length.
# (This is an MED requirement.)
o-pad-with-gaps $INFILE \
		-o $OUTDIR/pooled_filtered_med.fasta

# Finally, run the MED algorithm on the formatted fasta file   
printf "Running the MED pipeline..."
decompose $OUTDIR/pooled_filtered_med.fasta -o $OUTDIR

# Deactivate the med environment
source deactivate
