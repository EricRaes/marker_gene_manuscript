#! /bin/bash

#make a fasta file with size annotation from seqs in index col of a zotu table

#make a tab delimited file of the sequence and size annotations from the zotu table

echo 'creating tab file'

/your/path/to/file/otu.py ASV_abundance.txt > AASV_abundance.tab


#turn this table into fasta format

echo 'creating fasta'

awk '{printf(">%s\n%s\n",$2,$1);}' ASV_abundance.tab > ASV_abundance.fasta
