#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# This is a script to convert the MED 'NODE-REPRESENTATIVES.fasta'
# output file to the format required by 'USEARCH -uchime3_denovo'

# Set default parameter values
INFILE=""
OUTFILE=""

while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	-i|--infile)
	    INFILE="$2"
	    shift;;
	-o|--outfile)
	    OUTFILE="$2"
	    shift;;
	-h|--help)
	    printf "\nUSAGE: med2uchime.sh -i in_file -o out_file\n\n"
	    exit;;
	*)

	;;
    esac
    shift
done

sed -r '/^>/ s/\|(size):([0-9]+)/;\1=\2;/' <"$INFILE" >"$OUTFILE"
sed -ri '/^[^>]/ s/-+//' "$OUTFILE"
