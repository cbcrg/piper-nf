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
wget -q -O - http://get.nextflow.io > nextflow && chmod +x nextflow
./nextflow -bg \
  -daemon.interface eth0 \
  -daemon.tcp.socketTimeout 35s \
  -daemon.tcp.heartbeatFrequency 3s \
  -daemon.tcp.maxMissedHeartbeats 3 \
  -daemon.tcp.reconnectCount 15 \
  -daemon.tcp.networkTimeout 10s 
  

# save the environment for debugging 
env | sort > .boot.env

# Save this node IP address to a shared file 
ifconfig eth0 | grep "inet addr" | awk -F: '{print $2}' | awk '{print $1}' >> $REPO/nodes.txt
