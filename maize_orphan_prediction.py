from pyrpipe import sra,mapping,assembly,qc,tools,quant
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
import sys
from gffutils.iterators import DataIterator

#maizeRun=['SRR1573523','SRR999058','SRR520999','SRR1168424','SRR1621015','SRR3084882','SRR1620828','SRR3053545','SRR1620949','SRR1620947']
maizeRun=['SRR1573523','SRR999058']
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
    
sraObjects=[]

#download runs from SRA
for x in maizeRun:
    thisSraOb=sra.SRA(x,workingDir)
    if thisSraOb.download_fastq():
        sraObjects.append(thisSraOb)
    else:
        print("Download failed:"+x)


print("Following runs downloaded:")
for ob in sraObjects:
    print(ob.srr_accession)

#check fastq
for ob in sraObjects:
    print("SRR Accession: {}, fastq files: {}. {}".format(ob.srr_accession,ob.localfastq1Path,ob.localfastq2Path)) 
    if ob.fastqFilesExistsLocally():
          print("Both files exist!!")
    else:
          print("Error")
          raise Exception("Fastq files not found")

"""
# sra to fastq and quality trim
pathToAdapters="adapters2.fa"
bbdOpts={"ktrim":"r","k":"23","mink":"11","qtrim":"'rl'","trimq":"10","--":("-Xmx10g",),"ref":pathToAdapters}
bbdOb=qc.BBmap(**bbdOpts)
for ob in sraObjects:
    ob.perform_qc(bbdOb)


# run hisat2
hsOpts={"--dta-cufflinks":"","-p":"16"}
hs=mapping.Hisat2(hisat2_index=workingDir+"/maizeIndex/maizeInd",**hsOpts)
hisat2_buildArgs={"-p":"16","-a":"","-q":""}

if hs.check_index():
    print("Index {} exists".format(hs.hisat2_index))
else:
	if hs.build_index(workingDir+"/maizeIndex","maizeInd",GENOME,**hisat2_buildArgs) :
		print("Indexing done.")
	else:
		print("Index failed")
		sys.exit()


samList=[]
for ob in sraObjects:
    print("Processing {}...".format(ob.srr_accession))
    thisSam=hs.perform_alignment(ob,**{"-p":"16"}) 
    if thisSam:
        samList.append(thisSam)
print("Alignment done!! Sam files:"+ ",".join(samList))

#sam to sorted bam
samOb=tools.Samtools(**{"-@":"16"})
bamList=[]
i=0
for sam in samList:
    print("Processing:"+sam)
    thisBam=samOb.sam_sorted_bam(sam,delete_sam=True,delete_bam=True,objectid=sraObjects[i].srr_accession) 
    i+=1
    if thisBam:
        bamList.append(thisBam)
print("Sorted bam files:"+",".join(bamList))

#perform assembly
st=assembly.Stringtie()
gtfList=[]
i=0
for bam in bamList:
    print("Processing:"+bam)
    gtfList.append(st.perform_assembly(bam,objectid=sraObjects[i].srr_accession,**{"-m":"50","-p":"28"}))
    i+=1

print("Final GTFs:"+",".join(gtfList))

mergedGTF=st.stringtie_merge(*gtfList,**{"-p":"16"})


myGTF = mergedGTF
myFasta = GENOME

#not the fastest method
toWrite=[]
i=0
for feature in DataIterator(myGTF):
        if feature.featuretype == "transcript":
		print('processing gtf...')
                txName=feature.attributes['transcript_id'][0]
                gene=feature.attributes['gene_id'][0]
                thisSeq=">"+txName+" gene="+gene+" chr="+feature.seqid+" len="+str(len(feature))
                toWrite.append(thisSeq)
                toWrite.append(feature.sequence(myFasta))
                i+=1

f=open(workingDir+"/transcripts.fa","w")
f.write("\n".join(toWrite))
print("{} seq written".format(len(toWrite)/2))
"""


#start diamond
c1="wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -O "+workingDir+"/uniprot_sprot.fasta.gz"
pe.execute_command(c1.split(),logs=False)
c2="gunzip "+workingDir+"/uniprot_sprot.fasta.gz"
pe.execute_command(c2.split(),logs=False)

dindex=workingDir+"/dindex"
dm=tools.Diamond(index=dindex)
if not dm.check_index(dindex):
	print("Building diamond index")
	dm.build_index(workingDir+"/uniprot_sprot.fasta", "dindex", out_dir=workingDir, threads=20)

diamond_out="maize_diamond_out"
dm.run_align(workingDir+"/transcripts.fa", diamond_out, command="blastx", out_fmt=6, fmt_string="qseqid sseqid evalue pident", out_dir=workingDir, threads=25, objectid='diamond', **{"--more-sensitive":""})

#find tx which are not in diamond match
matches=set()
with open(workingDir+"/"+diamond_out) as f:
	lines=f.read().splitlines()

for l in lines:	
	matches.add(l.split('\t')[0])

#print(len(matches))

#read all ids in query
with open(workingDir+"/transcripts.fa") as f:
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
print("Orphan transcripts written to {}".format(workingDir+"/orphan_transcripts.fa"))


#use blastx with orphan transcripts for more sensitive results

#make blastdb
mkdb='makeblastdb -in '+ workingDir+"/uniprot_sprot.fasta"+' -dbtype prot -parse_seqids -out '+ workingDir+'/blastdb'
#pe.execute_command(mkdb.split(),logs=True)
blastOut=workingDir+'/blastout'
blastdb=workingDir+'/blastdb'
query=workingDir+"/orphan_transcripts.fa"
blastcmd='blastx -max_target_seqs 5 -num_threads 28 -query '+query+' -outfmt 6 -db '+blastdb+' -out '+ blastOut+' -evalue 1.0e-5'
#pe.execute_command(blastcmd.split())


#find final orphans from blast results
blastmatches=set()
with open(blastOut) as f:
	blastres=f.read().splitlines()
for l in blastres:
	blastmatches.add(l.split('\t')[0])

unm=matches.difference(blastmatches)
print("Orphans after blastxs: {}.".format(len(unm)))

#extract final orphans
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

f=open(workingDir+"/orphan_transcripts_final.fa","w")
f.write(towrite)
f.close()
print("Orphan transcripts written to {}".format(workingDir+"/orphan_transcripts_final.fa"))


#quantify expression with salmon
sl=quant.Salmon(salmon_index="")
sl.build_index(index_path=workingDir+"/salmonIndex",index_name="salIndex",fasta=workingDir+"/transcripts.fa")

for ob in sraObjects:
	sl.perform_quant(ob)








