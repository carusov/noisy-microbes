#!/bin/bash

### Author: Vincent Caruso
### Date: 9/20/2017
### Purpose: This script automates a typical BLAST search of a nucleotide
### sequence against the nt database. It is being written to query 16S
### amplicon data, so it uses the megablast (longer word) strategy. The output
### is in tabular format with comment lines (switch "7"), and the search is
### run remotely, on NCBI computers.
### Usage: blast_seqs.sh -i input_file [-o output_file]

# Parse command-line options
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	-i|--input)
	    INFILE="$2"
	    shift;;
	-o|--output)
	    OUTFILE="$2"
	    shift;;
	-t|--max_targets)
	    TARGETS="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: blast_seqs.sh -i in_file [-o out_file]\n"
	    exit;;
	*)

	;;
    esac
    shift
done

if [ -z "$OUTFILE" ]; then
    OUTFILE=${INFILE%.*}.txt
fi

blastn -query "$INFILE" -db nt -out "$OUTFILE" \
       -task megablast -max_target_seqs 10 \
       -outfmt "7 qseqid sseqid sskingdoms pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore" \
       -remote