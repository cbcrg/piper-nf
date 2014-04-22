#!/usr/bin/env bash
set -e
set -x 
REPO='/glusterfs/users/tcoffee/piper-nf/'
MAX=10
count=0
while [ $count -lt $MAX ]; do
  count=$(( $count + 1 ))
  ping -c1 index.docker.io >/dev/null  

  [[ $? -eq 0 || $count -eq $MAX ]] && break
  echo "No connection: sleeping for a while ... ($count)"
  sleep $(( $count * 5 ))
done

# Load the PIPE-R image
docker pull cbcrg/piper-nf

# Install NEXTFLOW and launch it 
export NXF_VER=0.7.1
wget -q -O - http://get.nextflow.io > nextflow && chmod +x nextflow
./nextflow -bg \
  -daemon.interface eth0 \
  -daemon.tcp.socketTimeout 120s \
  -daemon.tcp.heartbeatFrequency 5s \
  -daemon.tcp.maxMissedHeartbeats 5 \
  -daemon.tcp.reconnectCount 15 \
  -daemon.tcp.networkTimeout 30s \
  -daemon.tcp.ackTimeout 30s \
  -daemon.join file:/glusterfs/users/tcoffee/piper-nf/cluster/
  

# save the environment for debugging 
env | sort > .boot.env

# Save this node IP address to a shared file 
ifconfig eth0 | grep "inet addr" | awk -F: '{print $2}' | awk '{print $1}' >> $REPO/nodes.txt
