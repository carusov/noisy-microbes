#!/bin/bash

### Author: Vincent Caruso
### Date: 8/3/2017
### Purpose: This script implement the Deblur pipeline according to its
### recommended usage (see github.com/biocore/deblur).

# Set the default input file and output directory
INFILE=~/projects/thesis/data/dilution_w_blank/filtered/pooled_filtered_qiime.fasta
OUTDIR=~/projects/thesis/results/dilution_w_blank/deblur
REF_FILE=~/projects/thesis/references/silva_nr_v128_prokaryotes.fa

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
	-t|--trim)
	    TRIM="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: deblur_pipeline.sh -i input_file -o output_directory\n\n"
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

# First activate the deblur environment in miniconda
source activate deblur

# Now run the deblur workflow, using the QIIME-formatted file as input
printf "\nRunning the deblur workflow...\n"
deblur workflow --seqs-fp $INFILE \
       --output-dir $OUTDIR \
       -t $TRIM \
       --log-file $OUTDIR/deblur.log \
       --overwrite
#       --pos-ref-fp $REF_FILE

# Convert the .biom files to .txt files
if [ $(stat -c %s $OUTDIR/reference-hit.seqs.fa) != 0 ]; then
    biom convert -i $OUTDIR/reference-hit.biom \
	 -o $OUTDIR/reference-hit.txt \
	 --to-tsv
fi

if [ $(stat -c %s $OUTDIR/reference-non-hit.seqs.fa) != 0 ]; then
    biom convert -i $OUTDIR/reference-non-hit.biom \
	 -o $OUTDIR/reference-non-hit.txt \
	 --to-tsv
fi

if [ $(stat -c %s $OUTDIR/all.seqs.fa) != 0 ]; then
    biom convert -i $OUTDIR/all.biom \
	 -o $OUTDIR/all.txt \
	 --to-tsv
fi

#Deactivate the deblur environmnet
source deactivate
