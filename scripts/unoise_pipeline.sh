#!/bin/bash

### Author: Vincent Caruso
### Date: 7/1/2017
### This script implements the UNOISE pipeline according to its recommended usage
### (see www.drive5.com/usearch for detailed documentation)
### The input to the script is a .fastq file of merged, pooled, and quality-
### filtered sample reads


# Set the default input file, output directory, and raw merged read file
INFILE=~/thesis/data/dilution/filtered/pooled_filtered.fastq
OUTDIR=~/thesis/results/dilution/unoise
RAW_MERGED_FILE=~/thesis/data/dilution/merged/pooled_merged.fastq

# Parse command-line options
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	-i|--input)
	    INFILE="$2"
	    shift;;
	-o|--output)
	    OUTDIR="$2"
	    shift;;
	-r|--raw-merged)
	    RAW_MERGED_FILE="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: unoise_pipeline.sh -i input_file -o output_directory -r raw_merged_reads_file\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done

printf "\nINPUT FILE = ""${INFILE}""\n"
printf "OUTPUT DIRECTORY = ""${OUTDIR}""\n\n"

# Create the output directory, if necessary
if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR
fi

# Dereplicate the reads
usearch -fastx_uniques $INFILE \
	-fastqout $OUTDIR/pooled_uniques.fastq \
	-sizeout \
	-relabel Uniq

# Now denoise the reads into ZOTUs using the UNOISE2 algorithm
usearch -unoise3 $OUTDIR/pooled_uniques.fastq \
	-zotus $OUTDIR/zotus.fa \
	-ampout $OUTDIR/amplicons.fa \
	-tabbedout $OUTDIR/unoise3.txt

# USEARCH bug workaround: change 'Zotu' to 'OTU' in 'zotus.fa'
sed -i '/^>/ s/Zotu/OTU/' $OUTDIR/zotus.fa

# Use the generated OTUs to create an OTU table
usearch -otutab $RAW_MERGED_FILE \
	-zotus $OUTDIR/zotus.fa \
	-otutabout $OUTDIR/zotu_table.txt \
	-biomout $OUTDIR/zotu_table.biom \
	-mapout $OUTDIR/map.txt \
	-dbmatched $OUTDIR/zotus_with_sizes.fa \
	-notmatchedfq $OUTDIR/unmapped_reads.fastq \
	-sizeout

# Find out how many reads didn't map to OTUs
usearch -fastx_info $OUTDIR/unmapped_reads.fastq \
	-output $OUTDIR/unmapped_reads_info.txt

# How many "hiqh quality" reads don't map to OTUS?
usearch -fastq_filter $OUTDIR/unmapped_reads.fastq \
	-fastq_maxee 2.0 \
	-fastaout $OUTDIR/unmapped_hiqual.fa \
	-fastaout_discarded $OUTDIR/unmapped_loqual.fa

# Filter out the ZOTUS from the predicted amplicons, leaving just the predicted chimeras
usearch -search_exact $OUTDIR/amplicons.fa \
	-db $OUTDIR/zotus.fa \
	-strand plus \
	-notmatched $OUTDIR/chimeras.fa

# Combine predicted OTUs and chimeras into a single database
cat $OUTDIR/zotus.fa $OUTDIR/chimeras.fa \
    > $OUTDIR/otus_chimeras.fa

# Now, see how many high-quality, unmapped sequences are within 5% of an OTU or chimeric sequence
usearch -usearch_global $OUTDIR/unmapped_hiqual.fa \
	-db $OUTDIR/otus_chimeras.fa \
	-strand plus \
	-id 0.95 \
	-matched $OUTDIR/unmatched_noisy.fa \
	-notmatched $OUTDIR/unmatched_hiqual_other.fa

# Final coverage check: see if any leftover high-quality reads map to a large database
#usearch -usearch_global $OUTDIR/unmatched_hiqual_other.fa \
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
