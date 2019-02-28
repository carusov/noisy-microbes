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


# Process dilution series with separate samples
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing Zymo dilution series with separated samples...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
merge_and_filter.sh -w $DATA/dilution/zymo \
		    -f 230 -b 210 \
		    -s 221 -l 225 \
		    -d $MAXDIFFS  \
		    -e $MAXEE -n $MAXN

run_all_pipelines.sh -i $DATA/dilution/zymo \
		     -o $RESULTS/dilution/dilution_separate \
		     -s 221 -l 225 \
		     -F 2.5 -R 2.5 \

		     
# Process dilution series with pooled samples
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing Zymo dilution series with pooled samples...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
run_all_pipelines.sh -i $DATA/dilution/zymo \
		     -o $RESULTS/dilution/dilution_pooled \
		     -s 221 -l 225 \
		     -F 2.5 -R 2.5 \
		     --pooled
# DADA2 pipeline filter run with maxEE = c(2.5, 2.5)

