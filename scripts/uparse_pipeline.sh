#!/bin/bash

### Author: Vincent Caruso
### Date: 7/1/2017
### This script implements the UPARSE pipeline according to its recommended usage
### (see www.drive5.com/usearch for detailed documentation)
### The input to the script is a .fastq file of merged, pooled, and quality-
### filtered sample reads

DATA=~/projects/thesis/data

# Dereplicate the reads
usearch -fastx_uniques $DATA/clean/pooled_filtered.fastq \
	-fastqout $DATA/clean/pooled_uniques.fastq \
	-sizeout \
	-relabel Uniq

# Now cluster the reads into OTUs using the UPARSE algorithm
usearch -cluster_otus $DATA/clean/pooled_uniques.fastq \
	-otus $DATA/uparse/otus.fa \
	-relabel OTU

# Use the generated OTUs to create an OTU table
usearch -otutab $DATA/clean/pooled_merged.fastq \
	-otus $DATA/uparse/otus.fa \
	-otutabout $DATA/uparse/otu_table.txt \
	-biomout $DATA/uparse/otu_table.biom \
	-mapout $DATA/uparse/map.txt \
	-dbmatched $DATA/uparse/otus_with_sizes.fa \
	-notmatchedfq $DATA/uparse/unmapped_reads.fastq \
	-sizeout

# Find out how many reads didn't map to OTUs
usearch -fastx_info $DATA/uparse/unmapped_reads.fastq \
	-output $DATA/uparse/unmapped_reads_info.txt

# How many "hiqh quality" reads don't map to OTUS?
usearch -fastq_filter $DATA/uparse/unmapped_reads.fastq \
	-fastq_maxee 2.0 \
	-fastaout $DATA/uparse/unmapped_hiqual.fa \
	-fastaout_discarded $DATA/uparse/unmapped_loqual.fa

# Get sequences that are predicted to be chimeras
usearch -unoise3 $DATA/clean/pooled_uniques.fastq \
	-zotus $DATA/uparse/zotus.fa \
	-ampout $DATA/uparse/amplicons.fa \
	-tabbedout $DATA/uparse/unoise3.txt

# Filter out the ZOTUS from the predicted amplicons, leaving just the predicted chimeras
usearch -search_exact $DATA/uparse/amplicons.fa \
	-db $DATA/uparse/zotus.fa \
	-strand plus \
	-notmatched $DATA/uparse/chimeras.fa 
	
# Combine predicted OTUs and chimeras into a single database
cat $DATA/uparse/otus.fa $DATA/uparse/chimeras.fa \
    > $DATA/uparse/otus_chimeras.fa

# Now, see how many high-quality, unmapped sequences are within 5% of an OTU or chimeric sequence
usearch -usearch_global $DATA/uparse/unmapped_hiqual.fa \
	-db otus_chimeras.fa \
	-strand plus \
	-id 0.95 \
	-matched $DATA/uparse/unmatched_noisy.fa \
	-notmatched $DATA/uparse/unmatched_hiqual_other.fa

# Final coverage check: see if any leftover high-quality reads map to a large database
#usearch -usearch_global $DATA/uparse/unmatched_hiqual_other.fa \
#	-db silva.udb \
#	-strand both \
#	-idd 0.99 \
#	-alnout unmapped_silva.aln

### I might want to run `uchime2_ref` reference-based chimera checking with `-mode high_confidence`
### to see if it identifies any more chimeras

# usearch -uchime_ref reads.fastq \
#         -db silva.udb \
#         -uchimeout out.txt \
#         -strand plus \
#         -mode high_confidence

### Also, do I want to check for "tight" OTUs (OTUs with >97% identity)? If so I could use the UCLUST
### algorithm:
# usearch -cluster_fast otus.fa \
#         -id 0.97 \
#         -maxaccepts 4 \
#         -maxrejects 128 \
#         -top_hit_only
#         -uc hits.uc
# grep"^H" hits.uc | cut -f4 | sort -g

# Finally, assign taxonomy to the OTUS (I probably want to use the same assignment function for all
# methods being benchmarked, but I'll try this one here for fun).
#usearch -sintax otus.fa \
#	-db silva.udb \
#	-tabbedout sintax_taxonomy.txt \
#	-strand both \
#	-sintax_cutoff 0.8
