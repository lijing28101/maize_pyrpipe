import yaml
import sys
import os
from pyrpipe import sra,qc,mapping,tools,assembly,quant
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
from pyrpipe.runnable import Runnable


with open ("SRRids.txt") as f:
        srrid=f.read().splitlines()
	
workingDir="maize_out"
if not pu.check_paths_exist(workingDir):
    pu.mkdir(workingDir)

#Genome files	
GENOME=workingDir+"/Zm-B73-REFERENCE-NAM-5.0.fa"
if not pu.check_files_exist(GENOME):
    print("Downloading genome fasta file")
    wget="wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz -q -O "+GENOME+".gz"
    pe.execute_command(wget.split(),verbose=True,logs=False)
    pe.execute_command(['gunzip',GENOME+".gz"],verbose=True,logs=False)
    

# sra to fastq and quality trim
pathToAdapters="adapters2.fa"
bbdOpts={"ktrim":"r","k":"23","mink":"11","qtrim":"'rl'","trimq":"10","--":("-Xmx10g",),"ref":pathToAdapters}
BBDUK=qc.BBmap(**bbdOpts)

# run hisat2
hsOpts={"--dta-cufflinks":"","-p":"16"}
hs=mapping.Hisat2(index=workingDir+"/maizeIndex/maizeInd")
hisat2_buildArgs={"-p":"16","-a":"","-q":""}

if hs.check_index():
    print("Index {} exists".format(hs.hisat2_index))
else:
	if hs.build_index(workingDir+"/maizeIndex","maizeInd",GENOME,**hisat2_buildArgs) :
		print("Indexing done.")
	else:
		print("Index failed")
		sys.exit()

#perform assembly
stOpts={"-m":"150","-p":"16"}
st=assembly.Stringtie(**stOpts)

#download sra data, quality trimming, alignment and assemly
gtfList=[]
for x in srrid:
	gtfList.append(sra.SRA(x,directory=workingDir).trim(BBDUK).align(hs).assemble(st).gtf)


#merge stringtie output
filepath=workingDir+"/gtflist.txt"
f=open(filepath,"w")
f.write("\n".join(gtfList))
f.close()
mergedgtf=workingDir+"/gtflist_stringtieMerge.gtf"
STMERGE=assembly.Stringtie()
STMERGE.run(filepath,subcommand='--merge',**{'-p':"16",'-o':mergedgtf})

#get transcript cDNA fasta
gffread_out=workingDir+"/transcripts.fa"
gffread="gffread -w "+gffread_out+" -g "+GENOME+" "+mergedgtf
pe.execute_command(gffread.split(),verbose=True,logs=False)

#start blast
target=workingDir+"/uniprot_sprot.fasta"
blastout=workingDir+"/blastout.txt"
blastdb=Runnable(command='makeblastdb')
dbparams={'-in':input.target,'-dbtype':"prot",'-parse_seqids':"",'-out':"{wd}/blastdb"}
blastdb.run(**dbparams)
blastx=Runnable(command='blastx')
blastparams={'-max_target_seqs':"5",'-num_threads':"16",'-query':input.transcripts,'-outfmt':"6",'-db':"{wd}/blastdb",'-out':output.blastout,'-evalue':"0.001"}
blastx.run(**blastparams)


#find tx which are not in diamond match
matches=set()
with open(blastout) as f:
        lines=f.read().splitlines()

for l in lines:
        matches.add(l.split('\t')[0])

#print(len(matches))

#read all ids in query
with open(gffread_out) as f:
        lines=f.read().splitlines()
allQ=set()
for l in lines:
        if ">" in l:
                thisq=l.split()[0].split(">")[1]
                allQ.add(thisq)
unm=allQ.difference(matches)
print("Total queries: {}. Matched: {}. Orphans: {}".format(len(allQ),len(matches),len(unm)))

#extract all orphans
#use Biopython for more efficient handling of seq data
towrite=""
writeF=False
for l in lines:
        if ">" in l:
                thisq=l.split()[0].split(">")[1]
                if thisq in unm:
                        writeF=True
                else:
                        writeF=False
        if writeF:
                towrite+=l+"\n"

f=open(workingDir+"/orphan_transcripts.fa","w")
f.write(towrite)
f.close()
print("Orphan transcripts written to {}".format(workingDir+"/orphan_transcripts_final.fa"))


#quantify expression with salmon
salmon=quant.Salmon(index=output.index,transcriptome=input.transcripts,threads=16)
sra.SRA(srrid,directory=DIR).quant(salmon)









