#!/bin/bash
set -e 

on_exit() {
  exit_status=$?
  [[ $exit_status != 0 ]] && echo "*** VALIDATION FAILED ***"
  exit $exit_status
}

trap on_exit EXIT
        
NXF_CMD=${NXF_CMD:-nextflow}
NXF_RUN="$NXF_CMD run ." 
[[ $WITH_DOCKER ]] && NXF_RUN+=' -with-docker'

rm -rf db work result

#
# run normal mode 
#
$NXF_RUN | tee stdout

[[ 2 == $(grep -c 'Submitted process > formatChr' stdout) ]] || false
[[ 2 == $(grep -c 'Submitted process > formatBlast' stdout) ]] || false
[[ 2 == $(grep -c 'Submitted process > blast' stdout) ]] || false
[[ 2 == $(grep -c 'Submitted process > exonerate' stdout) ]] || false
[[ 2 == $(grep -c 'Submitted process > normExonerate' stdout) ]] || false
[[ 5 == $(grep -c 'Submitted process > similarity' stdout) ]] || false
[[ 1 == $(grep -c 'Submitted process > matrix' stdout) ]] || false


[[ 1 == $(grep -c 'F35C12.2b,86.67,100.00' stdout) ]] || false
[[ 1 == $(grep -c 'F35C1eee,83.28,100.00' stdout) ]] || false
[[ 1 == $(grep -c 'T24D1.1a.1,0,100.00' stdout) ]] || false
[[ 1 == $(grep -c 'T24D1.1b,0,100.00' stdout) ]] || false
[[ 1 == $(grep -c 'W04A8.1a,0,100.00' stdout) ]] || false


#
# run resume mode 
#
$NXF_RUN -resume | tee stdout

[[ 2 == $(grep -c 'Stored process > formatChr' stdout) ]] || false
[[ 2 == $(grep -c 'Stored process > formatBlast' stdout) ]] || false
[[ 2 == $(grep -c 'Cached process > blast' stdout) ]] || false
[[ 2 == $(grep -c 'Cached process > exonerate' stdout) ]] || false
[[ 2 == $(grep -c 'Cached process > normExonerate' stdout) ]] || false
[[ 5 == $(grep -c 'Cached process > similarity' stdout) ]] || false
[[ 1 == $(grep -c 'Cached process > matrix' stdout) ]] || false


[[ 1 == $(grep -c 'F35C12.2b,86.67,100.00' stdout) ]] || false
[[ 1 == $(grep -c 'F35C1eee,83.28,100.00' stdout) ]] || false
[[ 1 == $(grep -c 'T24D1.1a.1,0,100.00' stdout) ]] || false
[[ 1 == $(grep -c 'T24D1.1b,0,100.00' stdout) ]] || false
[[ 1 == $(grep -c 'W04A8.1a,0,100.00' stdout) ]] || false

