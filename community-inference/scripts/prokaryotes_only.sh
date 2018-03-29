#!/bin/bash

### Author: Vincent Caruso
### Date: 8/14/2017
### Purpose: Filter a fasta taxonomy database file to include only
### Bacteria and Archaea sequences
### Usage: prokaryotes_only.sh input_file output_file

# Delete all identifiers and sequences that are not Bacteria or Archaea
sed '/^>Bacteria\|^>Archaea/,+1! d' $1 > $2
