#!/bin/bash 

#
# A generic BLAST command wrapper
#
# CLI Parameters:
# - $1: Blast strategy to use
# - $2: Blast DB
# - $3: Query file

set -e
set -u

case "$1" in
'ncbi-blast')
blastn -db $2/db -query $3 -outfmt '6 qseqid sseqid evalue score qgi bitscore length nident positive mismatch pident ppos qacc gaps gaopen qaccver qlen qframe qstart qend sframe sstart send'
;;

'wu-blast')
wu-blastn $2/db $3 -mformat=2 -e 0.00001 -cpus 1 -filter=seg -lcfilter
;;

*) echo "Not a valid BLAST strategy: $1"; exit 1
;;

esac 

 
