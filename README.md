PIPER-R NF
==========

A pipeline for the detection and mapping of long non-coding RNAs


Quick start
-----------

Make sure you have installed the required dependencies listed below. 
Clone the git repository on your computer with the following command: 

    $ git clone git@github.com:cbcrg/piper-nf.git


When done, move in the project root folder just created `piper-nf`, 
it contains an example dataset in the `tutorial` folder. 

Launch the pipeline execution by entering the following command 
on your shell terminal:

    $ ./nextflow run piper.nf

By default the pipeline is executed against the provided example dataset. Check the *Pipeline parameters* section below
to see how enter your data on the program command line.

Pipeline parameters
-------------------

#####query

  * The query transcripts file in multi-fasta format
  * Example: `nextflow piper.nf --query=/some/path/query.fa`

#####genomes-file

  * The file listing the full paths to the genomes files
  * Example: `nextflow piper.nf --genomes-file=my-genomes.txt`


#####genomes-db

  * The location where the BLAST formatted *DB* are stored
  * Example: `nextflow piper.nf --genomes-db=/my/db/path`


#####query-chunk-size

  * Number of sequences in each chunck in which is sliced the query file
  * Example: `nextflow piper.nf --query-chunk-size=50`


#####result-dir

  * The location where the result files are stored.
  * Please note: if the folder exists, the all existing content will be deleted without further notice
  * Example: `nextflow piper.nf --result-dir=./my-result/`


#####blast-strategy

  * Which BLAST program to be used, `ncbi-blast` (default) or `wu-blast`
  * Example: `nextflow piper.nf --blast-strategy=wu-blast`


Run using Vagrant
-----------------

To avoid having to install all the pipeline dependencies, you may test the pipeline using 
the provided Vagrant VM, which downloads and configures all the required pieces 
of software for you. See http://www.vagrantup.com for more details about Vagrant.

The Vagrant environment uses the Ubuntu Precise 64 virtual machine, if you don't have it 
in your Vagrant boxes list, it will be downloaded automatically. 

To launch the VM move to the pipeline root folder `grape-nf` and enter the following command:
  
    $ vagrant up 


When it boots up and the configuration steps are terminated, login into the VM instance 
and move to the Grape pipeline folder 

    $ vagrant ssh 
    $ cd /vagrant
    
Now you can launch the pipeline as shown: 

	 $ ./nextflow piper.nf


When finished, stop the VM using the command `vagrant halt` or `vagrant destroy`, depending if you
want to temporary stop the execution or delete permanently the VM with all its files. 


Cluster support
---------------

Piper-NF execution relies on [Nextflow](http://nextflow-project.org) framework which provides an abstraction between
the pipeline functional logic and the underlying processing system.

Thus it is possible to execute it on your computer or any cluster resource
manager without modifying it.

Currently the following clusters are supported:

  + Oracle Grid Engine (SGE)
  + SLURM (beta)
  + Platform LSF (beta)


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


Run in the OSDC platform
-------------------------------

Piper can be executed in the  [Open Science Data Cloud](https://www.opensciencedatacloud.org) platform.

Launch instances of type `cbcrg_v.1.2`. One of them will used as *master* node i.e. where the 
application will be started, the remaining nodes will run the pipeline tasks. 

The *master* and *executor* nodes have to be instantiated separately. 

When launching the *executor* nodes specify the following code in the *Customization Script* field 
in order to configure the cluster:
    
    #!/bin/bash
    curl -fsSL https://raw.github.com/cbcrg/piper-nf/master/cloud-boot.sh | bash


Dependencies
------------
 
 * NCBI BLAST+ - http://blast.ncbi.nlm.nih.gov/
 * T-Coffee - http://tcoffee.org
 * Exonerate 2.2 - http://www.ebi.ac.uk/~guy/exonerate/ 
 * chr_subseq - unknown author 
 * gnu-sed (Mac OSX only)
 * gnu-coreutils (Mac OSX only)
