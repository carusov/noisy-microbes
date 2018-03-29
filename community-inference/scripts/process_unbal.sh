#!/bin/bash

DATA=~/thesis/data
RESULTS=~/thesis/results

# Process dataset metaID-49 from D'Amore, et al.
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing D'Amore unbalanced large sample...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
merge_and_filter.sh -w $DATA/damore_unbal_lg \
		    -f 240 -b 230 \
		    -s 263 -l 268 \
		    -d 30 -p 80 \
		    -e 2.0 -n 0

run_all_pipelines.sh -i $DATA/damore_unbal_lg -o $RESULTS/damore_unbal_lg \
		     -f 240 -b 230 \
		     -s 263 -l 268


# Process dataset metaID-90 from D'Amore, et al.
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
printf "\nProcessing D'Amore unbalanced small sample...\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n"
merge_and_filter.sh -w $DATA/damore_unbal_sm \
		    -f 240 -b 230 \
		    -s 258 -l 263 \
		    -d 30 -p 80 \
		    -e 2.0 -n 0

run_all_pipelines.sh -i $DATA/damore_unbal_sm -o $RESULTS/damore_unbal_sm \
		     -f 240 -b 230 \
		     -s 258 -l 263
