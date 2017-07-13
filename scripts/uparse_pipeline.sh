#!/bin/bash

### Author: Vincent Caruso
### Date: 7/1/2017
### This script implements the UPARSE pipeline according to its recommended usage
### (see www.drive5.com/usearch for detailed documentation)

DATA_PATH=~/projects/thesis/data

### For now, merging happens on the raw reads. I plan to update this to use reads trimmed by the
### DADA2 command `fastqPairedFilter`. This can't be done effectively using USEARCH, because
### `fastq_filter` doesn't process forward and reverse reads in tandem, so paired order is not
### preserved.

### I might want to consider using the USEARCH low-complexity read filter, which filters low-complexity
### reads that are potentially the result of seqeuncing past the end of the reverse adapter.
# usearch -filter_lowc reads_R1.fastq \
#         -reverse reads_R2.fastq \
#         -output filtered_R1.fastq \
#         -output2 filtered_R2.fastq \
#         -tabbedout lowc.txt \
#         -hitsout lowc.fastq

### Also, was there any PhiX spike-in used in sequencing that might not have been completely removed?

### And were the primer sequences removed before we received the "raw" data? They don't seem to be
### present in the reads.

# First merge forward and reverse reads from each sample, and pool the merged reads
usearch -fastq_mergepairs $DATA_PATH/raw/*R1.fastq \
	-fastqout $DATA_PATH/clean/pooled_merged.fastq \
	-fastq_maxdiffs 100 \
	-fastq_pctid 50 \
	-fastq_minmergelen 245 \
	-fastq_maxmergelen 255 \
	-relabel @ -report $DATA_PATH/reports/merge_report.txt

# Create a report of the expected errors of the merged reads
usearch -fastq_eestats2 $DATA_PATH/clean/pooled_merged.fastq \
	-output $DATA_PATH/reports/merged_eestats.txt \
	-length_cutoffs 200,*,10

# Quality filter the reads using a maximum of 2.0 expected errors, after trimming nucleotides from the
# beginning and end of the reads
usearch -fastq_filter $DATA_PATH/clean/pooled_merged.fastq \
	-fastqout $DATA_PATH/clean/pooled_filtered.fastq \
	-fastq_maxee 2.0

# Dereplicate the reads
usearch -fastx_uniques $DATA_PATH/clean/pooled_filtered.fastq \
	-fastqout $DATA_PATH/clean/pooled_uniques.fastq \
	-sizeout \
	-relabel Uniq

# Now cluster the reads into OTUs using the UPARSE algorithm
usearch -cluster_otus $DATA_PATH/clean/pooled_uniques.fastq \
	-otus $DATA_PATH/uparse_otus/otus.fa \
	-relabel OTU

# Use the generated OTUs to create an OTU table
usearch -otutab $DATA_PATH/clean/pooled_merged.fastq \
	-otus $DATA_PATH/uparse_otus/otus.fa \
	-otutabout $DATA_PATH/uparse_otus/otu_table.txt \
	-biomout $DATA_PATH/uparse_otus/otu_table.biom \
	-mapout $DATA_PATH/uparse_otus/map.txt \
	-dbmatched $DATA_PATH/uparse_otus/otus_with_sizes.fa \
	-notmatchedfq $DATA_PATH/uparse_otus/unmapped_reads.fastq \
	-sizeout

# Find out how many reads didn't map to OTUs
usearch -fastx_info $DATA_PATH/uparse_otus/unmapped_reads.fastq \
	-output $DATA_PATH/uparse_otus/unmapped_reads_info.txt

# How many "hiqh quality" reads don't map to OTUS?
usearch -fastq_filter $DATA_PATH/uparse_otus/unmapped_reads.fastq \
	-fastq_maxee 2.0 \
	-fastaout $DATA_PATH/uparse_otus/unmapped_hiqual.fa \
	-fastaout_discarded $DATA_PATH/uparse_otus/unmapped_loqual.fa

# Get sequences that are predicted to be chimeras
usearch -unoise3 $DATA_PATH/clean/pooled_uniques.fastq \
	-zotus $DATA_PATH/uparse_otus/zotus.fa \
	-tabbedout $DATA_PATH/uparse_otus/unoise3.txt \
	-chimeras $DATA_PATH/uparse_otus/chimeras.fa

# Combine predicted OTUs and chimeras into a single database
cat $DATA_PATH/uparse_otus/otus.fa $DATA_PATH/uparse_otus/chimeras.fa \
    > $DATA_PATH/uparse_otus/otus_chimeras.fa

# Now, see how many high-quality, unmapped sequences are within 5% of an OTU or chimeric sequence
usearch -usearch_global $DATA_PATH/uparse_otus/unmapped_hiqual.fa \
	-db otus_chimeras.fa \
	-strand plus \
	-id 0.95 \
	-matched $DATA_PATH/uparse_otus/unmatched_noisy.fa \
	-notmatched $DATA_PATH/uparse_otus/unmatched_hiqual_other.fa

# Final coverage check: see if any leftover high-quality reads map to a large database
#usearch -usearch_global $DATA_PATH/uparse_otus/unmatched_hiqual_other.fa \
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
