#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

WDIR=""

# Parse command-line options
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	-w|--working_dir)
	    WDIR="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: process_data.sh [-w --working working_directory] [options]\n"
	    printf "\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done

if [[ -z "$WDIR" ]]
then
    echo "ERROR: You must specify the top-level working directory."
    exit 1
fi


WDIR=$(readlink -f "$WDIR")


printf "\n**********************************************************************"
printf "\n**********************************************************************\n\n"
printf "\nProcessing high biomass samples...\n\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n\n"
process_high_biomass.sh -w "$WDIR"

printf "\n**********************************************************************"
printf "\n**********************************************************************\n\n"
printf "\nProcessing dilution series samples...\n\n"
printf "\n**********************************************************************"
printf "\n**********************************************************************\n\n"
process_dilution.sh -w "$WDIR"
