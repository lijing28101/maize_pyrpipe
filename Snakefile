import yaml
import sys
import os
from pyrpipe import sra,qc,mapping,tools,assembly
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe


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

##StringTie1##
stringtie1_params=documents['stringtie1']

##StringTie2##
stringtie2_params=documents['stringtie2']


#####create objects######

#SRA
SOBJS={}
for s in SRR:
        SOBJS[s]=sra.SRA(s,DIR)
        #OBLIST.append(sra.SRA(s,DIR))

##BBDUK##
BBDUK=qc.BBmap(**bbduk_params)
##HISAT2##
HISAT=mapping.Hisat2(hisat2_index="hisat_index/maizeInd",**hisat2_params)
if not HISAT.check_index():
        #build hisat2 index
        pu.print_green("Building Hisat2 index")
        HISAT.build_index("hisat_index","maizeInd",GENOME,**{"-p":"16","-a":"","-q":""})
##SAMTOOLS
ST=tools.Samtools(**{"-@":str(THREADS)})

##StringTie1##
STIE1=assembly.Stringtie(**stringtie1_params)


rule all:
        input:
                #expand("{wd}/{sample}/{sample}_1_bbduk.fastq",sample=SRR,wd=DIR,num=[1,2]),
                #expand("{wd}/{sample}/{sample}_2_bbduk.fastq",sample=SRR,wd=DIR,num=[1,2])
                #expand("{wd}/{sample}/{sample}_hisat2.sam",sample=SRR,wd=DIR,num=[1,2])
                #expand("{wd}/{sample}/{sample}_hisat2_sorted_stringtie1.gtf",sample=SRR,wd=DIR,num=[1,2]),
                #expand("{wd}/{sample}/{sample}_hisat2_sorted_stringtie_merged.gtf",sample=SRR,wd=DIR,num=[1,2]),
                #expand("{wd}/transcripts.fa",wd=DIR,num=[1,2]),
                expand("{wd}/maize.vioplot.jpeg",wd=DIR,num=[1,2])


rule download:
        output:
                f1="{wd}/{sample}/{sample}_1_bbduk.fastq",
                f2="{wd}/{sample}/{sample}_2_bbduk.fastq"
                #fid="{sample}"


        threads: THREADS

        run:
                fqfile=str(output.f1)
                srr=pu.get_file_basename(fqfile).split("_")[0]
                SOBJS[srr].download_fastq(procs=int(threads))
                SOBJS[srr].perform_qc(BBDUK,deleteRawFastq=False)

rule hisat:
        input:
                "{wd}/{sample}/{sample}_1_bbduk.fastq",
                "{wd}/{sample}/{sample}_2_bbduk.fastq"

        output:
                sam="{wd}/{sample}/{sample}_hisat2.sam"

        threads: THREADS

        run:
                samfile=str(output.sam)
                srr=pu.get_file_basename(samfile).split("_")[0]
                HISAT.perform_alignment(SOBJS[srr],**{"-p":str(threads)})

rule samtoBam:
        input:
                sam="{wd}/{sample}/{sample}_hisat2.sam"
        output:
                bam="{wd}/{sample}/{sample}_hisat2_sorted.bam"
        threads: THREADS
        run:
                print({input.sam})
                bamf=str({output.bam})
                samf=str({input.sam})
                srr=pu.get_file_basename(bamf).split("_")[0]
                ST.sam_sorted_bam(input.sam,delete_sam=True,delete_bam=True,objectid=srr)


rule stringtie1:
        input:
                bam="{wd}/{sample}/{sample}_hisat2_sorted.bam"
        output:
                gtf="{wd}/{sample}/{sample}_hisat2_sorted_stringtie.gtf"

        threads: THREADS

        run:
                f=str(output.gtf)
                srr=pu.get_file_basename(f).split("_")[0]
                STIE1.perform_assembly(input.bam,objectid=srr)

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
                STMERGE.stringtie_merge(input.gtflist,**{"-p":"16"})

rule stringtie2:
        input:
                bam="{wd}/{sample}/{sample}_hisat2_sorted.bam",
                mergedgtf="{wd}/gtflist_stringtieMerge.gtf"
        output:
                gtf="{wd}/{sample}/{sample}_hisat2_sorted_stringtie_merged.gtf",
                abundance="{wd}/{sample}/{sample}_abundance.tab"

        threads: THREADS

        run:
                f=str(output.gtf)
                srr=pu.get_file_basename(f).split("_")[0]
                STIE2=assembly.Stringtie(reference_gtf=input.mergedgtf,**stringtie2_params)
                STIE2.perform_assembly(input.bam,out_suffix="_stringtie_merged",objectid=srr)


rule gffread:
        input:
                mergedgtf="{wd}/gtflist_stringtieMerge.gtf"
        output:
                transcripts="{wd}/transcripts.fa"

        threads: THREADS

        run:
                cmd="gffread -w {wd}/transcripts.fa -g maize_out/Zm-B73-REFERENCE-NAM-5.0.fa {wd}/gtflist_stringtieMerge.gtf"
                pe.execute_command(cmd.split(" "),verbose=False,quiet=False,logs=True,command_name="gffread")



rule get_orphan:
        input:
                transcripts="{wd}/transcripts.fa",
                target="{wd}/combined.faa"
        output:
                orphan="{wd}/orphan.summary",
                nonorphan="{wd}/non.orphan.summary"

        shell:
                "{wd}/blast.sh"

rule summary:
        input:
                orphan="{wd}/orphan.summary",
                nonorphan="{wd}/non.orphan.summary"
        output:
                "{wd}/maize.vioplot.jpeg"

        shell:
                "Rscript DIR/plot.R"
