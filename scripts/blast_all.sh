#!/bin/bash

### This is a script to blast sequences from all six 16S processing methods at
### once.

blast_seqs.sh -i uclust_seqs.fasta -o uclust_blast.txt
blast_seqs.sh -i uparse_seqs.fasta -o uparse_blast.txt
blast_seqs.sh -i unoise_seqs.fasta -o unoise_blast.txt
blast_seqs.sh -i med_seqs.fasta -o med_blast.txt
blast_seqs.sh -i deblur_seqs.fasta -o deblur_blast.txt
blast_seqs.sh -i dada2_seqs.fasta -o dada2_blast.txt
