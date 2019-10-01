#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

### This script is mainly meant to serve as a record for how each
### dataset was processed, including the final parameters. It can
### also be run to re-process all datasets if desired.


WDIR=""
MAXDIFFS=10
MAXEE=2.0
MAXN=0


# Parse command-line options
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	-w|--working_dir)
	    WDIR="$2"
	    shift;;
	-d|--maxdiffs)
	    MAXDIFFS="$2"
	    shift;;
	-e|--maxee)
	    MAXEE="$2"
	    shift;;
	-n|--maxns)
	    MAXNS="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: process_dilution.sh [-w --working working_directory] [options]\n"
	    printf "\nOptions: \t[default]"
	    printf "\n-e --maxee \t["$MAXEE"] \t\t\t\t\tmax_expected_errors"
	    printf "\n-d --maxdiffs \t["$MAXDIFFS"] \t\t\t\t\tmax_differences"
	    printf "\n-n --maxns \t["$MAXN"] \t\t\t\t\tmax_Ns"
	    printf "\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done


WDIR=$(readlink -f "$WDIR")
DATA="$WDIR"/data
RESULTS="$WDIR"/results


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

