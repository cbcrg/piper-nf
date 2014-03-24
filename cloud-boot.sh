#!/usr/bin/env bash
set -e
REPO='/glusterfs/users/tcoffee/piper-nf/'

# Load the PIPE-R image
docker pull cbcrg/piper-nf

# Install NEXTFLOW and launch it 
wget -q -O - http://get.nextflow.io > nextflow && chmod +x nextflow
./nextflow -d -daemon.interface eth0  

# Save this node IP address to a shared file 

ifconfig eth0 | grep "inet addr" | awk -F: '{print $2}' | awk '{print $1}' >> $REPO/nodes.txt
