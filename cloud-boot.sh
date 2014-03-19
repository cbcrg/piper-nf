#!/usr/bin/env bash
REPO=/glusterfs/users/tcoffee/piper-nf/

# Load the PIPE-R image
docker pull cbcrg/piper-nf

# Copy Nextflow and launch the daemon
cp $REPO/nextflow $HOME
chmod +x $HOME/nextflow
nohup $HOME/nextflow -daemon.interface eth0 &> $HOME/log & 

# Save this node IP address to a shared file 

ifconfig eth0 | grep "inet addr" | awk -F: '{print $2}' | awk '{print $1}' >> REPO/nodes.txt