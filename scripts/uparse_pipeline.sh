#!/bin/bash

usearch -fastq_mergepairs raw/*R1.fastq -fastqout clean/all_merged.fastq \
	-fastq_maxdiffs 100 -fastq_pctid 50 -fastq_minmergelen 245 -fastq_maxmergelen 255 \
	-relabel @ -report reports/merge_report.txt

usearch -fastq_eestats2 clean/all_merged.fastq -output reports/merged_eestats.txt \
	-length_cutoffs 200,*,10 -ee_cutoffs 1.0,2.0,2.5,3.0,3.5,4.0

usearch -fastq_filter clean/all_merged.fastq -fastqout clean/all_filtered.fastq \
	-fastq_maxee 3.0

usearch -fastx_uniques clean/all_filtered.fastq -fastqout clean/all_uniques.fastq \
	-sizeout -relabel Uniq

usearch -cluster_otus clean/all_uniques.fastq -otus uparse_otus/otus.fa \
	-relabel OTU
