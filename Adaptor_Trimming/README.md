**Purpose**

The script is used to do adaptor sequenc trimming from the raw RNAseq data.
We use a concatenation of common bacterial genomes that are available on NCBI.
#The file is called bacterialgenome.fa

**Usage**

Install Trimmomatics from https://github.com/timflutre/trimmomatic
Edit the path to the software in the script
Adjust the script depending on the type of data (pair ended etc)
Open the Terminal from the directory containing the sequence files
Run the script by typing: bash bbsplit.sh
The script is a loop that takes all the samples from the folder 
Script can also be done for samples one by one.
