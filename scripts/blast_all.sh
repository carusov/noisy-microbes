#!/bin/bash

### This is a script to blast sequences from all six 16S processing methods at
### once.

for f in *_seqs.fasta;
do
    name=${f%_seqs.fasta}
    blast_seqs.sh -i $f -o $name"_blast.txt" -t 10
    echo "Done with "$name" BLAST search"
done


#blast_seqs.sh -i uclust_seqs.fasta -o uclust_blast.txt -t 10
#echo "Done with UCLUST BLAST search"
#blast_seqs.sh -i uparse_seqs.fasta -o uparse_blast.txt -t 10
#echo "Done with UPARSE BLAST search"
#blast_seqs.sh -i unoise_seqs.fasta -o unoise_blast.txt -t 10
#echo "Done with UNOISE BLAST search"
#blast_seqs.sh -i med_seqs.fasta -o med_blast.txt -t 10
#echo "Done with MED BLAST search"
#blast_seqs.sh -i deblur_seqs.fasta -o deblur_blast.txt -t 10
#echo "Done with Deblur BLAST search"
#blast_seqs.sh -i dada2_seqs.fasta -o dada2_blast.txt -t 10
#echo "Done with DADA2 BLAST search"
