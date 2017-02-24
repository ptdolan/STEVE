import os
import numpy as np
from matplotlib import pyplot as pp

charKey={"a":1,"g":1,"t":1,"c":1,"y":1,"r":1,"n":1,"d":1,"e":1,"f":1,"g":1,"h":1,"i":1,"k":1,"l":1,"m":1,"n":1,"p":1,"q":1,"r":1,"s":1,"t":1,"w":1,"v":1,"x":2,"*":2,"-":0}


##########################################################

def readDir(eveDir,pattern):
  faList=[]
  for root,dirs,files in os.walk(eveDir):
    for F in files:
      if pattern in F:
        faList.append(root+"/"+F)
  #print(faList)
  return(faList)

def convert(sequences):
  converted=[[charKey[c] for c in sequence] for sequence in sequences]
  converted=np.array(converted)
  print(converted)
  return(converted)
  
def openFile(F):
  names=[]
  print(F)
  with open(F,'r') as IF:
    lines=[line.strip() for line in IF]
  newSeq=1
  sequences=[]
  sequence=[]
  for line in lines:
    if line.startswith(">"):
      if len(sequence)>0:
        #print("writing",line)
        sequences.append(sequence)
      sequence=[]
      names.append(line.replace(">",""))
      newSeq=1
    elif line.startswith(">")!=True and newSeq==1:
      sequence=line.strip().lower()
      newSeq=0
    else:
      sequence="".join([sequence,line.strip().lower()])
      #print(sequence)
  return(sequences)


########################################################## 

eveDir="/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data"
pattern=".aln"

dirList=readDir(eveDir,pattern)

for F in dirList:
  sequences = openFile(F)
  aarray=convert(sequences)
  print(np.shape(aarray)[0])
  if np.shape(aarray)[0]>1:
    pp.matshow(aarray)
    coverage = np.sum(aarray,axis=0)
    cleansequences = [[seq[char] for char in range(len(seq)) if coverage[char]>20] for seq in sequences]
    after = convert(cleansequences)
    pp.matshow(after)
    pp.show()

