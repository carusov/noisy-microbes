#!/bin/bash

echo "This script generates summary stats on all .fastq files in the current working directory."

if [ ! -d fastq_info ]; then
    mkdir fastq_info
fi

for fq in *.fastq
do
    bn=$(basename $fq .fastq)
    usearch -fastx_info $fq -output fastq_info/$bn"_info.txt"
done
