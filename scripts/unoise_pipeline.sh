#!/bin/bash

### Author: Vincent Caruso
### Date: 7/1/2017
### This script implements the UNOISE pipeline according to its recommended usage
### (see www.drive5.com/usearch for detailed documentation)
### The input to the script is a .fastq file of merged, pooled, and quality-
### filtered sample reads

DATA=~/projects/thesis/data
RESULTS=~/projects/thesis/results

if [ ! -d "$RESULTS/unoise" ]; then
    mkdir $RESULTS/unoise
fi

# Dereplicate the reads
usearch -fastx_uniques $DATA/clean/pooled_filtered.fastq \
	-fastqout $DATA/clean/pooled_uniques.fastq \
	-sizeout \
	-relabel Uniq

# Now denoise the reads into ZOTUs using the UNOISE2 algorithm
usearch -unoise3 $DATA/clean/pooled_uniques.fastq \
	-zotus $RESULTS/unoise/zotus.fa \
	-ampout $RESULTS/unoise/amplicons.fa \
	-tabbedout $RESULTS/unoise/unoise3.txt

# Filter out the ZOTUS from the predicted amplicons, leaving just the predicted chimeras
usearch -search_exact $RESULTS/unoise/amplicons.fa \
	-db $RESULTS/unoise/zotus.fa \
	-strand plus \
	-notmatched $RESULTS/unoise/chimeras.fa

# Use the generated OTUs to create an OTU table
usearch -otutab $DATA/clean/pooled_merged.fastq \
	-otus $RESULTS/unoise/zotus.fa \
	-otutabout $RESULTS/unoise/zotu_table.txt \
	-biomout $RESULTS/unoise/zotu_table.biom \
	-mapout $RESULTS/unoise/map.txt \
	-dbmatched $RESULTS/unoise/zotus_with_sizes.fa \
	-notmatchedfq $RESULTS/unoise/unmapped_reads.fastq \
	-sizeout

# Find out how many reads didn't map to OTUs
usearch -fastx_info $RESULTS/unoise/unmapped_reads.fastq \
	-output $RESULTS/unoise/unmapped_reads_info.txt

# How many "hiqh quality" reads don't map to OTUS?
usearch -fastq_filter $RESULTS/unoise/unmapped_reads.fastq \
	-fastq_maxee 2.0 \
	-fastaout $RESULTS/unoise/unmapped_hiqual.fa \
	-fastaout_discarded $RESULTS/unoise/unmapped_loqual.fa

# Combine predicted OTUs and chimeras into a single database
cat $RESULTS/unoise/zotus.fa $RESULTS/unoise/chimeras.fa \
    > $RESULTS/unoise/otus_chimeras.fa

# Now, see how many high-quality, unmapped sequences are within 5% of an OTU or chimeric sequence
usearch -usearch_global $RESULTS/unoise/unmapped_hiqual.fa \
	-db $RESULTS/unoise/otus_chimeras.fa \
	-strand plus \
	-id 0.95 \
	-matched $RESULTS/unoise/unmatched_noisy.fa \
	-notmatched $RESULTS/unoise/unmatched_hiqual_other.fa

# Final coverage check: see if any leftover high-quality reads map to a large database
#usearch -usearch_global $RESULTS/unoise/unmatched_hiqual_other.fa \
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
