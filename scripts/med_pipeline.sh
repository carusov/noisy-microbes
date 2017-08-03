#!/bin/bash

DATA=~/projects/thesis/data
RESULTS=~/projects/thesis/results

### First open the MED environment using miniconda
source activate med

### Now, use 'sed' to convert a fasta file from the USEARCH label format
### (using the '-relabel @' option) to the MED format.
printf "Converting filtered fasta file to MED-friendly format...\n"
sed '/^>/ s/\./_/' $DATA/clean/pooled_filtered.fasta \
    > $DATA/clean/pooled_filtered_med.fasta

### Pad shorter reads with gaps so all reads have the same length.
### (This is an MED requirement.)
o-pad-with-gaps $DATA/clean/pooled_filtered_med.fasta \
		-o $DATA/clean/pooled_filtered_med_padded.fasta
mv $DATA/clean/pooled_filtered_med_padded.fasta $DATA/clean/pooled_filtered_med.fasta

### Make sure the output directory exists
if [ ! -d "$RESULTS/med" ]; then
   mkdir $RESULTS/med
fi
   
### Finally, run the MED algorithm on the formatted fasta file   
printf "Running the MED pipeline..."
decompose $DATA/clean/pooled_filtered_med.fasta -o $RESULTS/med
