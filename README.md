PIPER-R NF
==========

A pipeline for the detection and mapping of long non-coding RNAs


Quick start
-----------

Make sure you have all the required dependencies listed in the last section. 

Install the Nextflow runtime by running the following command: 

    $ curl -fsSL get.nextflow.io | bash


When done, you can launch the pipeline execution by entering the command shown below:

    $ ./nextflow run cbcrg/piper-nf


By default the pipeline is executed against the provided example dataset. Check the *Pipeline parameters* section below
to see how specify your input data on the program command line.

Pipeline parameters
-------------------

#####query

  * The query transcripts file in multi-fasta format
  * Example: `nextflow run piper-nf --query=/some/path/query.fa`

#####genomes-file

  * The file listing the full paths to the genomes files
  * Example: `nextflow run piper-nf --genomes-file=my-genomes.txt`


#####genomes-db

  * The location where the BLAST formatted *DB* are stored
  * Example: `nextflow run piper-nf --genomes-db=/my/db/path`


#####query-chunk-size

  * Number of sequences in each chunck in which is sliced the query file
  * Example: `nextflow run piper-nf --query-chunk-size=50`


#####result-dir

  * The location where the result files are stored.
  * Please note: if the folder exists, the all existing content will be deleted without further notice
  * Example: `nextflow run piper-nf --result-dir=./my-result/`


#####blast-strategy

  * Which BLAST program to be used, `ncbi-blast` (default) or `wu-blast`
  * Example: `nextflow run piper-nf --blast-strategy=wu-blast`


Run with Docker 
---------------- 

Piper-nf dependecies are also distributed by using a [Docker](http://www.docker.com) container which frees you from 
the installation and configuration of all the pieces of software required by Piper-nf. 

The Piper-nf Docker image is published at this address https://registry.hub.docker.com/u/cbcrg/piper-nf/

If you have Docker in installed in your computer pull this image by entering the following command: 

    $ docker pull cbcrg/piper-nf
  
  
After that you will be able to run Piper-nf using the following command line: 

    $ ./nextflow run cbcrg/piper-nf -with-docker




Cluster support
---------------

Piper-NF execution relies on [Nextflow](http://nextflow.io) framework which provides an abstraction between
the pipeline functional logic and the underlying processing system.

Thus it is possible to execute it on your computer or any cluster resource
manager without modifying it.

Currently the following clusters are supported:

  + Oracle/Univa/Open Grid Engine (SGE)
  + Platform LSF
  + SLURM
  + PBS/Torque


By default the pipeline is parallelized by spanning multiple threads in the machine where the script is launched.

To submit the execution to a SGE cluster create a file named `nextflow.config`, in the directory
where the pipeline is going to be launched, with the following content:

    task {
      processor='sge'
      queue='<your queue name>'
    }

In doing that, tasks will be executed through the `qsub` SGE command, and so your pipeline will behave like any
other SGE job script, with the benefit that *Nextflow* will automatically and transparently manage the tasks
synchronisation, file(s) staging/un-staging, etc.

Alternatively the same declaration can be defined in the file `$HOME/.nextflow/config`.

To lean more about the avaible settings and the configuration file read the Nextflow documentation 
 http://www.nextflow.io/docs/latest/config.html


Dependencies
------------
 * Java VM 7 (or higher) - http://www.oracle.com/technetwork/java/javase/downloads/
 * NCBI BLAST+ - http://blast.ncbi.nlm.nih.gov/
 * T-Coffee - http://tcoffee.org
 * Exonerate 2.2 - http://www.ebi.ac.uk/~guy/exonerate/ 
 * chr_subseq - unknown author 
 * gnu-sed (Mac OSX only)
 * gnu-coreutils (Mac OSX only)
