import yaml
import sys
import os
from pyrpipe import sra,qc,mapping,tools,assembly,quant
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
from pyrpipe.runnable import Runnable

####Read config#####
configfile: "config.yaml"
DIR = config['DIR']
THREADS=config['THREADS']

##check required files
GENOME=config['genome']

if not pu.check_files_exist(GENOME):
        pu.mkdir("maize_out")
        pe.execute_command('wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz -O maize_out/Zm-B73-REFERENCE-NAM-5.0.fa.gz'.split(),logs=False)
        pe.execute_command('gunzip maize_out/Zm-B73-REFERENCE-NAM-5.0.fa.gz'.split(),logs=False)

if pu.check_files_exist(GENOME):
        pu.print_green("Genome files...OK")
else:
        pu.print_boldred("Genome files...NOT FOUND")
        sys.exit(1)


#####Read SRR ids######
with open ("SRRids.txt") as f:
        SRR=f.read().splitlines()

#######Read Params###########
with open('params.yaml') as file:
       documents = yaml.full_load(file)

##bbduk
bbduk_params=documents['bbduk']

##hisat2
hisat2_params=documents['hisat2']

##StringTie##
stringtie_params=documents['stringtie']


#####create objects######

##BBDUK##
BBDUK=qc.BBmap(**bbduk_params)
##HISAT2##
HISAT=mapping.Hisat2(index="maize_out/maizeIndex/maizeInd",**hisat2_params)
##StringTie1##
ST=assembly.Stringtie(**stringtie_params)


rule all:
        input:
                expand("{wd}/{sample}/{sample}_hisat2_sorted_stringtie_merged.gtf",sample=SRR,wd=DIR),
                expand("{wd}/transcripts.fa",wd=DIR),
                expand("{wd}/blastout.txt",wd=DIR),
                expand("{wd}/orphan.id",wd=DIR),
                expand("{wd}/orphan_transcripts_final.fa",wd=DIR),
                expand("{wd}/salmon_index",wd=DIR),
                expand("{wd}/{sample}/salmon_out/quant.sf",wd=DIR,sample=SRR)

                
rule process:
        output:
                gtf="{wd}/{sample}/{sample}_hisat2_sorted_stringtie.gtf"

        threads: THREADS

        run:
                gtffile=str(output.gtf)
                srrid=gtffile.split("/")[1]
                sra.SRA(srrid,directory=DIR).trim(BBDUK).align(HISAT).assemble(ST)


rule create_GTFlist:
        input:
                g=expand("{wd}/{sample}/{sample}_hisat2_sorted_stringtie.gtf",sample=SRR,wd=DIR,num=[1,2])
        output:
                "{wd}/gtflist.txt"
        run:
                filepath=DIR+"/gtflist.txt"
                f=open(filepath,"w")
                f.write("\n".join(input.g))
                f.close()

rule stringtieMerge:
        input:
                gtflist="{wd}/gtflist.txt"
        output:
                mergedgtf="{wd}/gtflist_stringtieMerge.gtf"

        run:
                STMERGE=assembly.Stringtie()
                STMERGE.run(input.gtflist,subcommand='--merge',**{'-p':"16",'-o':output.mergedgtf})


rule gffread:
        input:
                mergedgtf="{wd}/gtflist_stringtieMerge.gtf"
        output:
                transcripts="{wd}/transcripts.fa"

        shell:
                "gffread -w {output.transcripts} -g {GENOME} {input.mergedgtf}"


rule blast:
        input:
                transcripts="{wd}/transcripts.fa",
                target="{wd}/uniprot_sprot.fasta"
        output:
                blastout="{wd}/blastout.txt"

        run:
                blastdb=Runnable(command='makeblastdb')
                dbparams={'-in':input.target,'-dbtype':"prot",'-parse_seqids':"",'-out':"{wd}/blastdb"}
                blastdb.run(**dbparams)
                blastx=Runnable(command='blastx')
                blastparams={'-max_target_seqs':"5",'-num_threads':"16",'-query':input.transcripts,'-outfmt':"6",'-db':"{wd}/blastdb",'-out':output.blastout,'-evalue':"0.001"}
                blastx.run(**blastparams)


rule get_orphan:
        input:
                transcripts="{wd}/transcripts.fa",
                blastout="{wd}/blastout.txt"
        output:
                orphanid="{wd}/orphan.id"

        shell:
                """
                grep ">" {input.transcripts} | cut -f1 -d" " | sed 's/>//' > {DIR}/id.all
                cut -f1 {input.blastout} | sort -u > {DIR}/nonorphan.id
                cat {DIR}/id.all {DIR}/nonorphan.id | sort | uniq -u > {output.orphanid}
                """


rule get_orphanFa:
        input:
               transcripts="{wd}/transcripts.fa",
               orphanid="{wd}/orphan.id"
        output:
               orphanFa="{wd}/orphan_transcripts_final.fa"

        shell:
               "seqtk subseq {input.transcripts} {input.orphanid} > {output.orphanFa}"


rule salmon_index:
        input:
                transcripts="{wd}/transcripts.fa"
        output:
                index=directory("{wd}/salmon_index")

        run:
                salmon=quant.Salmon(index=output.index,transcriptome=input.transcripts,threads=16)


rule salmon:
        input:
                index="{wd}/salmon_index"
        output:
                quant_file="{wd}/{sample}/salmon_out/quant.sf"

        run:
                salmon=quant.Salmon(index=input.index,threads=16)
                outfile=str(output.quant_file)
                srrid=outfile.split("/")[1]
                sra.SRA(srrid,directory=DIR).quant(salmon)



              
