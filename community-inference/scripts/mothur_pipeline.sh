#!/bin/bash

### Author: Vincent Caruso
### Date: 4/25/18
### Last modified: 4/30/18
### This script implements the mothur pipeline according to the mothur SOP
### (see www.mothur.org/wiki/MiSeq_SOP for detailed documentation)

# Set the default input file, output directory, and 16S reference db
OUTDIR=$PWD/mothur
GROUP_FILE=mothur.groups
REF_ALIGN=~/thesis/references/silva.bacteria/v4.fasta
TAX_DB=~/thesis/references/trainset9_032012.pds.fasta
TAX_LABEL=~/thesis/references/trainset9_032012.pds.tax


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
	-g|--groups)
	    GROUP_FILE="$2"
	    shift;;
	-r|--ref_align)
	    REF_ALIGN="$2"
	    shift;;
	-t|--tax_db)
	    TAX_DB="$2"
	    shift;;
	-l|--tax_labels)
	    TAX_LABEL="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: mothur_pipeline.sh -i input_file -o output_directory [options]"
	    printf "\nOptions \t\t[default]"
	    printf "\n-g --groups \t\t[./%s] \t\t\t\t\t\tmothur groups file" "$GROUP_FILE"
	    printf "\n-r --ref_align \t\t[%s] \t16S alignment file" "$REF_ALIGN"
	    printf "\n-t --tax_db \t\t[%s] \ttaxonomy database" "$TAX_DB"
	    printf "\n-l --tax_labels \t[%s] \ttaxonomy labels" "$TAX_LABEL"
	    printf "\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done

# Create the output directory, if necessary
if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
else
    rm "$OUTDIR"/*
fi

# Check for groups file
if [[ ! -f "$GROUP_FILE" ]]
then
    printf "\nERROR: The mothur groups file was not found. Please double-check the file path."
    printf "\nThe mothur groups file entered was: %s\n\n" "$GROUP_FILE"
    exit 1
else
    GROUP_FILE=$(readlink -f "$GROUP_FILE")
fi


INFILE=$(readlink -f "$INFILE")
OUTDIR=$(readlink -f "$OUTDIR")
REF_ALIGN=$(readlink -f "$REF_ALIGN")
TAX_DB=$(readlink -f "$TAX_DB")
TAX_LABEL=$(readlink -f "$TAX_LABEL")

printf "\nINPUT FILE: %s" "$INFILE"
printf "\nOUTPUT DIRECTORY: %s" "$OUTDIR"
printf "\nGROUP FILE: %s" "$GROUP_FILE"
printf "\nREFERENCE ALIGNMENT FILE: %s" "$REF_ALIGN"
printf "\nTAXONOMY DATABASE: %s" "$TAX_DB"
printf "\nTAXONOMY LABEL FILE: %s" "$TAX_LABEL"
printf "\n\n"

ln -s "$INFILE" "$OUTDIR"
ln -s "$GROUP_FILE" "$OUTDIR"/mothur.groups
name=$(basename "$INFILE" .fasta)
group=mothur.groups

# Mothur puts all outputs in the current directory, so change to the output directory
pushd "$OUTDIR"

### The rest of the pipeline is executed with mothur

# Get unique sequences
mothur "#unique.seqs(fasta = "$name".fasta)"

# Generate a count table
mothur "#count.seqs(name="$name".names, group="$group")"

# Extract just the V4 region of the reference alignment file
#mothur "#pcr.seqs(fasta="$REF_ALIGN", start=11894, end=25319, keepdots=F, processors=4)"
#v4=${REF_ALIGN//.fasta/.pcr.fasta}

# Align sequences to the V4 reference alignment
mothur "#align.seqs(fasta="$name".unique.fasta, reference="$REF_ALIGN")"

# Summarize the alignment
coord=($(mothur "#summary.seqs(fasta="$name".unique.align, count="$name".count_table)" \
		| awk '/Median/ {print $2, $3}'))

echo "$name" ${coord[@]} >> "$OUTDIR"/mothur.coords

# Filter out poorly aligned sequences
mothur "#screen.seqs(fasta="$name".unique.align, count="$name".count_table, summary="$name".unique.summary, start=${coord[0]}, end=${coord[1]}, maxhomop=8)"

# Summarize the filtered alignment
mothur "#summary.seqs(fasta="$name".unique.good.align, count="$name".good.count_table)"

# Filter out terminal gaps and gap-only columns from the alignment
mothur "#filter.seqs(fasta="$name".unique.good.align, vertical=T, trump=.)"

# Unique the sequences again, just in case anything changed
mothur "#unique.seqs(fasta="$name".unique.good.filter.fasta, count="$name".good.count_table)"

# Pre-cluster the sequences, allowing a maximum of 2 differences
mothur "#pre.cluster(fasta="$name".unique.good.filter.unique.fasta, count="$name".unique.good.filter.count_table, diffs=2)"

# Predict chimeras using VSEARCH
mothur "#chimera.vsearch(fasta="$name".unique.good.filter.unique.precluster.fasta, count="$name".unique.good.filter.unique.precluster.count_table, dereplicate=T)"

# Remove chimeras
mothur "#remove.seqs(fasta="$name".unique.good.filter.unique.precluster.fasta, accnos="$name".unique.good.filter.unique.precluster.denovo.vsearch.accnos)"

# Summarize the chimera-free sequences
mothur "#summary.seqs(fasta="$name".unique.good.filter.unique.precluster.pick.fasta, count="$name".unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)"

# Assign taxonomy to remaining sequences to find those that just don't belong
mothur "#classify.seqs(fasta="$name".unique.good.filter.unique.precluster.pick.fasta, count="$name".unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference="$TAX_DB", taxonomy="$TAX_LABEL", cutoff=80)"

# And now remove those undesirables
mothur "#remove.lineage(fasta="$name".unique.good.filter.unique.precluster.pick.fasta, count="$name".unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy="$name".unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"

# Update the taxonomy summary
mothur "#summary.tax(taxonomy="$name".unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, count="$name".unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"

# Compute distances between sequences
mothur "#dist.seqs(fasta="$name".unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.03)"

# Cluster sequences
mothur "#cluster(column="$name".unique.good.filter.unique.precluster.pick.pick.dist, count="$name".unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"

# Make us some shared, I mean, OTU table
mothur "#make.shared(list="$name".unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count="$name".unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)"

# Find representative sequences
mothur "#get.oturep(column="$name".unique.good.filter.unique.precluster.pick.pick.dist, list="$name".unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count="$name".unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta="$name".unique.fasta)"

# Edit representative sequence names to just OTU IDs
awk -F "\t" 'BEGIN {OFS = ""} /^>/ {split($2, otu, "|"); print ">", otu[1]}; /^[^>]/ {print $0}' "$name".unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.fasta > rep_seqs.fasta

# Create links with short names to key output files
ln -s "$name".unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared otu_table.txt
#ln -s "$name".unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.fasta rep_seqs.fasta

# Finally, change back to the directory from which the script was called
popd
