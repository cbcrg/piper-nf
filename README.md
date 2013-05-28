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
on your shell temrinal:

    $ ./nextflow piper.nf


Dependencies 
------------
 
 * NCBI BLAST+ - http://blast.ncbi.nlm.nih.gov/
 * T-Coffee - http://tcoffee.org
 * Exonerate 2.2 - http://www.ebi.ac.uk/~guy/exonerate/ 
 * chr_subseq - unknown author 
 * gnu-sed (Mac OSX only)
 * gnu-csplit (Mac OSX only)