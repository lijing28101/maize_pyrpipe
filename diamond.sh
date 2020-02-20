#!/bin/bash

module load diamond

#buid index for diamond
diamond makedb --in uniprot_sprot.fasta -d diamond_index/diamondDB -p 16

#transcript against protein database by diamond 
diamond blastx -d diamond_index/diamondDB -q transcripts.fa -p 16 -o diamond.txt -f 6

#extract orphan and non.orphan list
grep ">" transcripts.fa | cut -f1 -d" " | sed 's/>//' > all.gene
awk '$11<0.00001' diamond.txt | cut -f1 | sort -u > non.orphan.diamond
cat all.gene non.orphan.diamond | sort | uniq -u > orphan.diamond

