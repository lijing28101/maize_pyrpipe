#!/bin/bash

mkdir -p maize_out

wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz -O maize_out/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -O "+workingDir+"/uniprot_sprot.fasta.gz


cd maize_out
gunzip -f *.gz

#run hisat2 index
mkdir -p maizeIndex
hisat2-build -p 16 -a -q Zm-B73-REFERENCE-NAM-5.0.fa maizeIndex/maizeInd
