## Purpose
Update annotations (from the genome annotation) with insect-only gene functions, add annotation to transcripts not annotated in the genome. 
We use diamond (https://github.com/bbuchfink/diamond) to perform a local blast on the insect only database, that we created using blastdb.makedb (https://www.ncbi.nlm.nih.gov/books/NBK569841/).

## Process
1. Make a fasta file containing all sequences of interest (here sequences.fa)
2. Create a local insect database using blastdb by downloading the nr database and extracting the insects using the insect db id number (50557)
3. Blastx the sequences on the insect db
4. Open blast output into Blast2go for further analyzes
