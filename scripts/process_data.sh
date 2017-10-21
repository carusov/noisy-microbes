#!/usr/bin

### This script is mainly meant to serve as a record for how each
### dataset was processed, including the final parameters. It can
### also be run to re-process all datasets if desired.

DATA=~/projects/thesis/data
RESULTS=~/projects/thesis/results

# Process dilution series with extraction blank
merge_and_filter.sh -w $DATA/dilution_w_blank \
		    -f 230 -b 210 \
		    -s 220 -l 225 \
		    -d 30 -p 50 \
		    -e 2.0 -n 0

run_all_pipelines.sh -i $DATA/dilution_w_blank -o $RESULTS/dilution_w_blank \
		     -f 230 -b 210 \
		     -s 220 -l 225


# Process dilution series without extraction blank
merge_and_filter.sh -w $DATA/dilution \
		    -f 230 -b 210 \
		    -s 220 -l 225 \
		    -d 30 -p 50 \
		    -e 2.0 -n 0

run_all_pipelines.sh -i $DATA/dilution -o $RESULTS/dilution \
		     -f 230 -b 210 \
		     -s 220 -l 225


# Process only the "neat" sample from the dilution series
merge_and_filter.sh -w $DATA/neat \
		    -f 230 -b 210 \
		    -s 220 -l 225 \
		    -d 30 -p 50 \
		    -e 2.0 -n 0

run_all_pipelines.sh -i $DATA/neat -o $RESULTS/neat \
		     -f 230 -b 210 \
		     -s 220 -l 225


# Process the "balanced" dataset from Schirmer, et al.
merge_and_filter.sh -w $DATA/schirmer_balanced \
		    -f 240 -b 220 \
		    -s 258 -l 263 \
		    -d 30 -p 80 \
		    -e 2.0 -n 0

run_all_pipelines.sh -i $DATA/schirmer_balanced -o $RESULTS/schirmer_balanced \
		     -f 240 -b 220 \
		     -s 258 -l 263


# Process dataset 130403 from Kozich, et al.
merge_and_filter.sh -w $DATA/kozich_130403 \
		    -f 240 -b 220 \
		    -s 220 -l 225 \
		    -d 30 -p 80 \
		    -e 2.0 -n 0

run_all_pipelines.sh -i $DATA/kozich_130403 -o $RESULTS/kozich_130403 \
		     -f 240 -b 220 \
		     -s 220 -l 225

# Process dataset metaID-46 from D'Amore, et al.
merge_and_filter.sh -w $DATA/damore_balanced \
		    -f 250 -b 240 \
		    -s 258 -l 263 \
		    -d 30 -p 80 \
		    -e 2.0 -n 0

run_all_pipelines.sh -i $DATA/damore_balanced -o $RESULTS/damore_balanced \
		     -f 250 -b 240 \
		     -s 258 -l 263
