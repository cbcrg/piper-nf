#!/bin/bash 

#
# A generic wrapper a for BLAST formatting 
#
# It takes as input a genome file and format it 
# to BLAST db binary format 
#
# CLI Parameters:
# - $1: Blast strategy to use
# - $2: The genome FASTA file to format 
# - $3: The resulting BLAST DB

set -e
set -u

case "$1" in
'ncbi-blast')
makeblastdb -dbtype nucl -in $2 -out $3/db 
;;

'wu-blast')
xdformat -n -o $3/db - < $2
;;

*) echo "Not a valid BLAST strategy: $1"; exit 1
;;

esac 

 