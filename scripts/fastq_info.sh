#!/bin/bash

echo "This script generates summary stats on all .fastq files in the current working directory."

#alias usearch='usearch10.0.240_i86linux32'
mkdir -p ../fastq_info

for fq in *.fastq
do
#  usearch -fastx_info $fq -output ../fastq_info/$fq
  usearch10.0.240_i86linux32 -fastx_info $fq -output ../fastq_info/$fq
done
