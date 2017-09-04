#!/bin/bash

### Author: Vincent Caruso
### Date: 8/14/2017
### Purpose: This is a wrapper shell script that runs each of the
### 16S processing method scripts in sequence. This script specifies
### paths for the inputs and outputs to each script, and is thus
### specific to my configuration, so it is not convenient for use
### on someone else's machine at the moment.
### Usage: run_all_pipelines.sh

# Run the UCLUST pipeline with defaults
printf "\nRunning the UCLUST pipeline...\n\n"
uclust_pipeline.sh

# Run the UPARSE pipeline with defaults
printf "\nRunning the UPARSE pipeline...\n\n"
uparse_pipeline.sh

# Run the UNOISE pipeline with defaults
printf "\nRunning the UNOISE pipeline...\n\n"
unoise_pipeline.sh

# Run the MED pipeline with defaults
printf "\nRunning the MED pipeline...\n\n"
med_pipeline.sh

# Run the Deblur pipeline with defaults
printf "\nRunning the Deblur pipeline...\n\n"
deblur_pipeline.sh

# Run the DADA2 pipeline with defaults
printf "\nRunning the DADA2 pipeline...\n\n"
SCRIPTS=~/projects/thesis/noisy-microbes/scripts
Rscript $SCRIPTS/rmd2r.R
Rscript $SCRIPTS/dada2_pipeline.R
