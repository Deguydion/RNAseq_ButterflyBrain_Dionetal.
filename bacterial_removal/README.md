## Purpose

Remove contaminant sequences from bacterial origines.
We use bbsplit from https://github.com/BioInfoTools/BBMap/blob/master/sh/bbsplit.sh
We use a concatenation of common bacterial genomes that are available on NCBI.
The file is called bacterialgenome.fa

## Usage 

Install BBmap and bb split.
Edit the path to the software in the script.
Open the Terminal from the directory containing the sequence files.
Run the script by typing: bash trim.sh
The script is a loop that takes all the samples from the folder. 
Script can also be done for samples one by one.


