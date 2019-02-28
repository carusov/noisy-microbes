#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

### This script is mainly meant to serve as a record for how each
### dataset was processed, including the final parameters. It can
### also be run to re-process all datasets if desired.

DATA=~/repo_test/data
RESULTS=~/repo_test/results
MAXDIFFS=10
MAXEE=2.0
MAXN=0


printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing Zymo neat sample...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
merge_and_filter.sh -w $DATA/high_biomass/zymo \
		    -f 230 -b 210 \
		    -s 221 -l 225 \
		    -d $MAXDIFFS  \
		    -e $MAXEE -n $MAXN

run_all_pipelines.sh -i $DATA/high_biomass/zymo \
		     -o $RESULTS/high_biomass/zymo \
		     -s 221 -l 225 \
		     -F 2.5 -R 2.5 
# DADA2 pipeline filter run with maxEE = c(2.5, 2.5)


# Process dataset 130403 from Kozich, et al.
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing Kozich Mock 1 sample...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
merge_and_filter.sh -w $DATA/high_biomass/kozich \
		    -f 240 -b 220 \
		    -s 221 -l 225 \
		    -d $MAXDIFFS  \
		    -e $MAXEE -n $MAXN

run_all_pipelines.sh -i $DATA/high_biomass/kozich \
		     -o $RESULTS/high_biomass/kozich \
		     -s 221 -l 225 \
		     -F 2.5 -R 3.0		     
# DADA2 pipeline filter run with maxEE = c(2.5, 3.0)


# Process the "balanced" dataset from Schirmer, et al.
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing Schirmer balanced sample...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
merge_and_filter.sh -w $DATA/high_biomass/schirmer \
		    -f 240 -b 220 \
		    -s 259 -l 263 \
		    -d $MAXDIFFS  \
		    -e $MAXEE -n $MAXN

run_all_pipelines.sh -i $DATA/high_biomass/schirmer \
		     -o $RESULTS/high_biomass/schirmer \
		     -s 259 -l 263 \
		     -F 2.5 -R 2.5
# DADA2 pipeline filter run with maxEE = c(2.5, 2.5)
		     

# Process dataset metaID-90 from D'Amore, et al.
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing D'Amore uneven sample...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
merge_and_filter.sh -w $DATA/high_biomass/damore \
		    -f 250 -b 240 \
		    -s 260 -l 264 \
		    -d $MAXDIFFS  \
		    -e $MAXEE -n $MAXN

run_all_pipelines.sh -i $DATA/high_biomass/damore \
		     -o $RESULTS/high_biomass/damore \
		     -s 260 -l 264 \
		     -F 2.5 -R 2.5		     
# DADA2 pipeline filter run with maxEE = c(2.5, 2.5)
