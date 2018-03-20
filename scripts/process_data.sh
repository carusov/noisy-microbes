#!/bin/bash

### This script is mainly meant to serve as a record for how each
### dataset was processed, including the final parameters. It can
### also be run to re-process all datasets if desired.

DATA=~/thesis/data
RESULTS=~/thesis/results
MAXDIFFS=10
#PCTID=90
MAXEE=2.0
MAXN=0


printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing Zymo neat sample...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
merge_and_filter.sh -w $DATA/zymo_neat \
		    -f 230 -b 210 \
		    -s 221 -l 225 \
		    -d $MAXDIFFS  \
		    -e $MAXEE -n $MAXN

run_all_pipelines.sh -i $DATA/zymo_neat -o $RESULTS/zymo_neat \
		     -f 230 -b 210 \
		     -s 221 -l 225 \
		     -F 2.5 -R 2.5
# DADA2 pipeline filter run with maxEE = c(2.5, 2.5)


# Process the "balanced" dataset from Schirmer, et al.
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing Schirmer balanced sample...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
merge_and_filter.sh -w $DATA/schirmer_balanced \
		    -f 240 -b 220 \
		    -s 259 -l 263 \
		    -d $MAXDIFFS  \
		    -e $MAXEE -n $MAXN

run_all_pipelines.sh -i $DATA/schirmer_balanced -o $RESULTS/schirmer_balanced \
		     -f 240 -b 220 \
		     -s 259 -l 263 \
		     -F 2.5 -R 2.5
# DADA2 pipeline filter run with maxEE = c(2.5, 2.5)
		     

# Process dataset 130403 from Kozich, et al.
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing Kozich Mock 1 sample...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
merge_and_filter.sh -w $DATA/kozich_130403 \
		    -f 240 -b 220 \
		    -s 221 -l 225 \
		    -d $MAXDIFFS  \
		    -e $MAXEE -n $MAXN

run_all_pipelines.sh -i $DATA/kozich_130403 -o $RESULTS/kozich_130403 \
		     -f 240 -b 220 \
		     -s 221 -l 225 \
		     -F 2.5 -R 3.0		     
# DADA2 pipeline filter run with maxEE = c(2.5, 3.0)


# Process dataset metaID-90 from D'Amore, et al.
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing D'Amore uneven sample...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
merge_and_filter.sh -w $DATA/damore_uneven \
		    -f 250 -b 240 \
		    -s 260 -l 264 \
		    -d $MAXDIFFS  \
		    -e $MAXEE -n $MAXN

run_all_pipelines.sh -i $DATA/damore_uneven -o $RESULTS/damore_uneven \
		     -f 250 -b 240 \
		     -s 260 -l 264 \
		     -F 2.5 -R 2.5		     
# DADA2 pipeline filter run with maxEE = c(2.5, 2.5)


# Process dataset metaID-88 from D'Amore, et al.
#printf "\n**********************************************************************"
#printf "\n**********************************************************************\n"
#printf "\nProcessing D'Amore balanced sample...\n"
#printf "\n**********************************************************************"
#printf "\n**********************************************************************\n"
#merge_and_filter.sh -w $DATA/damore_balanced \
#		    -f 250 -b 240 \
#		    -s 261 -l 261 \
#		    -d $MAXDIFFS  \
#		    -e $MAXEE -n $MAXN

#run_all_pipelines.sh -i $DATA/damore_balanced -o $RESULTS/damore_balanced \
#		     -f 250 -b 240 \
#		     -s 261 -l 261


# Process dilution series without extraction blank, with pooled samples
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing Zymo dilution series (without blank), with pooled samples...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
merge_and_filter.sh -w $DATA/dilution \
		    -f 230 -b 210 \
		    -s 221 -l 225 \
		    -d $MAXDIFFS  \
		    -e $MAXEE -n $MAXN

run_all_pipelines.sh -i $DATA/dilution -o $RESULTS/dilution_pooled \
		     -f 230 -b 210 \
		     -s 221 -l 225 \
		     -F 2.5 -R 2.5 \
		     -m pooled
# DADA2 pipeline filter run with maxEE = c(2.5, 2.5)

# Process dilution series without extraction blank, with pooled samples
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing Zymo dilution series (without blank), with separated samples...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
run_all_pipelines.sh -i $DATA/dilution -o $RESULTS/dilution_separate \
		     -f 230 -b 210 \
		     -s 221 -l 225 \
		     -F 2.5 -R 2.5 \
		     -m separate
		     

# Process dilution series with extraction blank
#printf "\n**********************************************************************"
#printf "\n**********************************************************************\n"
#printf "\nProcessing Zymo dilution series with blank sample...\n"
#printf "\n**********************************************************************"
#printf "\n**********************************************************************\n"
#merge_and_filter.sh -w $DATA/dilution_w_blank \
#		    -f 230 -b 210 \
#		    -s 223 -l 223 \
#		    -d $MAXDIFFS  \
#		    -e $MAXEE -n $MAXN

#run_all_pipelines.sh -i $DATA/dilution_w_blank -o $RESULTS/dilution_w_blank \
#		     -f 230 -b 210 \
#		     -s 223 -l 223


