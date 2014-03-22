#!/usr/bin/env bash
REPO=/glusterfs/users/tcoffee/piper-nf/

# Load the PIPE-R image
docker pull cbcrg/piper-nf

# Install NEXTFLOW and launch it 
export NXF_PACK='hz'
wget -q -O -  get.nextflow.io | bash -x
./nextflow -d -daemon.interface eth0  

# Save this node IP address to a shared file 

ifconfig eth0 | grep "inet addr" | awk -F: '{print $2}' | awk '{print $1}' >> $REPO/nodes.txt
