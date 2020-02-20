from pyrpipe import sra,mapping,assembly,qc,tools
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe

maizeRun=['SRR1573523','SRR999058','SRR520999','SRR1168424','SRR1621015','SRR3084882','SRR1620828','SRR3053545','SRR1620949','SRR1620947']
workingDir="maize_out"
if not pu.check_paths_exist(workingDir):
    pu.mkdir(workingDir)
	
GENOME=workingDir+"/Zm-B73-REFERENCE-NAM-5.0.fa"
if not pu.check_files_exist(GENOME):
    print("Downloading genome fasta file")
    wget="wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz -q -O "+GENOME+".gz"
    pe.execute_command(wget.split(),verbose=True,logs=False)
    pe.execute_command(['gunzip',GENOME+".gz"],verbose=True,logs=False)
    
sraObjects=[]

for x in maizeRun:
    thisSraOb=sra.SRA(x,workingDir)
    if thisSraOb.download_fastq():
        sraObjects.append(thisSraOb)
    else:
        print("Download failed:"+x)


print("Following runs downloaded:")
for ob in sraObjects:
    print(ob.srr_accession)
	
pathToAdapters="adapters2.fa"
bbdOpts={"ktrim":"r","k":"23","mink":"11","qtrim":"'rl'","trimq":"10","--":("-Xmx2g",),"ref":pathToAdapters}
bbdOb=qc.BBmap(**bbdOpts)
for ob in sraObjects:
    ob.perform_qc(bbdOb)

for ob in sraObjects:
    print("SRR Accession: {}, fastq files: {}. {}".format(ob.srr_accession,ob.localfastq1Path,ob.localfastq2Path))
    
    if ob.fastqFilesExistsLocally():
          print("Both files exist!!")
    else:
          print("Error")
          raise Exception("Fastq files not found")
		  
hsOpts={"--dta-cufflinks":"","-p":"16"}
hs=mapping.Hisat2(hisat2_index="",**hsOpts)
hisat2_buildArgs={"-p":"16","-a":"","-q":""}

if hs.build_index(workingDir+"/maizeIndex","maizeInd",GENOME,**hisat2_buildArgs) :
    print("Indexing done.")

if hs.check_index():
    print("Index {} exists".format(hs.hisat2_index))

samList=[]
for ob in sraObjects:
    print("Processing {}...".format(ob.srr_accession))
    thisSam=hs.perform_alignment(ob,**{"-p":"16"}) 
    if thisSam:
        samList.append(thisSam)
print("Alignment done!! Sam files:"+ ",".join(samList))

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


st=assembly.Stringtie()
gtfList=[]
i=0
for bam in bamList:
    print("Processing:"+bam)
    gtfList.append(st.perform_assembly(bam,objectid=sraObjects[i].srr_accession,**{"-m":"150","-p":"16"}))
    i+=1

print("Final GTFs:"+",".join(gtfList))

filepath=workingDir+"/gtflist.txt"
f=open("gtflist.txt","w")
f.write("\n".join(gtfList))
f.close()

st.stringtie_merge(workingDir+"/gtflist.txt",**{"-p":"16"})

GTF=workingDir+"/gtflist_stringtieMerge.gtf"
st2=assembly.Stringtie(reference_gtf=GTF)
gtfList2=[]
i=0
for bam in bamList:
    abundance=workingDir+"/"+sraObjects[i].srr_accession+"/"+sraObjects[i].srr_accession+"_abundance.tab"
    print("Processing:"+bam)
    gtfList2.append(st2.perform_assembly(bam,objectid=sraObjects[i].srr_accession,out_suffix="_MergedST",**{"-p":"16","-e":"","-A":abundance}))
    i+=1

print("Final GTFs:"+",".join(gtfList2))

tx_file=workingDir+"/transcripts.fa"
cmd="gffread -w "+tx_file+" -g maize_out/Zm-B73-REFERENCE-NAM-5.0.fa maize_out/gtflist_stringtieMerge.gtf"
pe.execute_command(cmd.split(" "),verbose=False,quiet=False,logs=True,command_name="gffread")

