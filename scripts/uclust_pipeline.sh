#!/bin/bash

### Author: Vincent Caruso
### Date: 7/13/2017
### This script implements the UCLUST pipeline through the QIIME SOP
### (see www.qiime.org for detailed documentation)

# Set the default input file, output directory, and 16S reference db
INFILE=~/projects/thesis/data/dilution_w_blank/filtered/pooled_filtered_qiime.fasta
OUTDIR=~/projects/thesis/results/dilution_w_blank/uclust
REF_FILE=~/projects/thesis/references/gold.fa

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
	-r|--reference)
	    REF_FILE="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: uclust_pipeline.sh -i input_file -o output_directory -r 16S_reference_db.fasta\n\n"
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

# Activate the qiime environment in miniconda
source activate qiime1

### The rest of the pipeline is executed with QIIME scripts

# Identify chimeric sequences 
printf "\nIdentifying chimeric sequences...\n"
identify_chimeric_seqs.py -m usearch61 \
        -i $INFILE \
        -r $REF_FILE \
        -o $OUTDIR/

# Filter out identified chimeric sequences
printf "\nRemoving chimeric sequences from sample reads...\n"
filter_fasta.py -f $INFILE \
        -o $OUTDIR/pooled_nochim.fa \
        -s $OUTDIR/chimeras.txt \
        -n 

# Pick de novo OTUs
# I need to investigate the parameters that this QIIME script passes to UCLUST.
# Are sequences sorted by size?
printf "\nPicking OTUs (de novo) and representative seqeunces, assigning\n"
printf "taxonomy, performing multiple alignments, and building a\n"
printf "phylogenetic tree...\n"
pick_de_novo_otus.py -i $OUTDIR/pooled_nochim.fa \
        -o $OUTDIR/ \
        -f

# Convert the .biom file to a tabbed .txt file
biom convert -i $OUTDIR/otu_table.biom \
     -o $OUTDIR/otu_table.txt \
     --to-tsv

# Deactivate miniconda3 environment
source deactivate
