#!/bin/bash

### Author: Vincent Caruso
### Date: 8/3/2017
### Purpose: This script implement the Deblur pipeline according to its
### recommended usage (see github.com/biocore/deblur).

DATA=~/projects/thesis/data
RESULTS=~/projects/thesis/results
REF=/media/sf_Shared_Ubuntu/references

# First activate the deblur environment in miniconda
printf "\nActivating the deblur environment...\n"
source activate deblur

# Now run the deblur workflow, using the QIIME-formatted file as input
printf "\nRunning the deblur workflow...\n"
deblur workflow --seqs-fp $DATA/clean/pooled_filtered_qiime.fasta \
       --output-dir $RESULTS/deblur \
       -t 230 \
       --pos-ref-fp $DATA/references/silva_nr_v128_train_set.fa \
       --log-file $RESULTS/deblur/deblur.log \
       --overwrite

# Convert the .biom files to .txt files
biom convert -i $RESULTS/deblur/reference-hit.biom \
     -o $RESULTS/deblur/reference-hit.txt \
     --to-tsv
biom convert -i $RESULTS/deblur/reference-non-hit.biom \
     -o $RESULTS/deblur/reference-non-hit.txt \
     --to-tsv
biom convert -i $RESULTS/deblur/all.biom \
     -o $RESULTS/deblur/all.txt \
     --to-tsv
