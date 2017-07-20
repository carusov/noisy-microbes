#!/bin/bash

echo "This script generates summary stats on all .fastq files in the current working directory."

mkdir -p ../fastq_info

for fq in *.fastq
do
  usearch -fastx_info $fq -output ../fastq_info/$fq
done
