#pairwisesort.py
#take fas, output sorted GN fas for alignment. 

from Bio.Seq import Seq
import sys
import os
from Bio import SeqIO
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pp

importDir = sys.argv[1]

def seqdict(importDir):
  faList=[]
  names=[]
  seqs=[]
  seqDict={}
  for root,dirs,files in os.walk(importDir):
    for F in files:
      print(F)
      if ".fa" in F:
        faList.append(root+"/"+F)
  #print(faList)
  allseqs=[]
  for F in faList:
    header=(F.split("/")[-1].replace(".fa",""))
    fasta_sequences = SeqIO.parse(open(F),'fasta')
    seqsList = [[fasta.id,fasta.id.split("GN=")[1][:4].replace("_",""),str(fasta.seq)] for fasta in fasta_sequences]
    allseqs.extend([seq[:] for seq in seqsList])
  subunits=set([s[1] for s in allseqs])
  for sU in subunits:
    #print(sU)
    suOut=[">"+line[0]+"\n"+line[2]+"\n" for line in allseqs if line[1]==sU ]
    #print(suOut)
    with open(sU+"_Out.fa",'w') as OF:
    	for su in suOut:
       		OF.write(su)

seqDict=seqdict(importDir)
