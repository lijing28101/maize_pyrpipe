#!/bin/bash

module load blast-plus

#buid index for blast, protein from 7 closded species
makeblastdb -in combined.faa -title combined.db -dbtype prot -out combined.db -parse_seqids

#transcript against protein database by blast
blastx -query transcripts.fa -db combined.db -num_threads 16 -outfmt '6 qseqid sseqid qstart qend sstart send evalue score' >> blast.txt

#extract orphan and non.orphan list
grep ">" transcripts.fa | cut -f1 -d" " | sed 's/>//' > all.gene
awk '$7<0.00001' blast.txt | cut -f1 | sort -u > non.orphan.blast
cat all.gene non.orphan.blast | sort | uniq -u > orphan.blast

#extract orphan and non-orphan sequence as fasta file
module load seqtk
seqtk subseq transcripts.fa orphan.blast > orphan.fa
seqtk subseq transcripts.fa non.orphan.blast > non.orphan.fa

#length and GC content summary
seqtk comp orphan.fa | cut -f1,2,3,4,5,6 | awk -F'\t' 'BEGIN {OFS=FS}  $7=($4+$5)/$2*100' > orphan.summary
seqtk comp non.orphan.fa | cut -f1,2,3,4,5,6 | awk -F'\t' 'BEGIN {OFS=FS}  $7=($4+$5)/$2*100' > non.orphan.summary
