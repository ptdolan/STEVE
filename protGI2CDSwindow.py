#protGI2CDSwindow.py
import os
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import matplotlib as mpl
from matplotlib import pyplot as pp

path="/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/EVEphyloAnalysis/"
eveDir=path+"ContigEVEFastas/"
infiles= os.listdir(path+"ResultFiles/") 
infiles = [f for f in infiles if ".result" in f]
#print (infiles)

def evedict(eveDir):
  faList=[]
  names=[]
  seqs=[]
  eveDict={}
  for root,dirs,files in os.walk(eveDir):
    for F in files:
      if ".fasta" in F:
        faList.append(root+"/"+F)
  #print(faList)
  for F in faList:
    print(F)
    fasta_sequences = SeqIO.parse(open(F),'fasta')

    eveDict.update({fasta.id:fasta.seq for fasta in fasta_sequences})
  #print(len(eveDict.keys()))
  print("Created EVE dict with "+str(len(eveDict.keys()))+" keys.")
  return(eveDict)
  
def extractEVEs(infile):
  with open(path+"ResultFiles/"+infile,'r') as IF:
    eveinfo = [line for line in IF]
    #print(eveinfo)
    EVEid = [line.split("\t")[0] for line in eveinfo]
    #print(EVEid)
    EVEcoord =  [[int(line.split("\t")[1]),int(line.split("\t")[2]) ]for line in eveinfo]
    EVEstrand =  [coord[0]<coord[1] for coord in EVEcoord]

    EVEstart =  [min(coord) for coord in EVEcoord]
    EVEend =  [max(coord) for coord in EVEcoord]
    EVEregister =  [line.split("\t")[5] for line in eveinfo]
    #print(EVEregister)


    protIDlist = [line.split("\t")[3].split("|") for line in eveinfo]
    protID=[]
    for R in protIDlist:
    #  print(R)
      if len(R)>1:
        protID.append(R[3])
      elif len(R)==1:
        protID.append(R[0].split(" ")[0])
    viruscoord = [[int(line.split("\t")[8]),int(line.split("\t")[9])] for line in eveinfo] #
    virusstrand = [interval[1]>interval[0] for interval in viruscoord] #
    virusstart = [min(interval) for interval in viruscoord] #
    virusend = [max(interval) for interval in viruscoord]

    return(EVEid,protID,EVEstart,EVEend,virusstart,virusend,virusstrand,EVEstrand,EVEregister)

##
## Use blast hit info to make profile alignments of the EVE hits
##
def idCDS(eveDict,protID,EVEstart,EVEend,EVEid,virusstart,virusend,vstrand,estrand,infile,EVEregister):
  entrycount=0
  ### Entrez search 
  Entrez.email = "ptdolan@stanford.edu"       #Always introduce yourself to NCBI.

  outputFasta=path+infile+"input.fa"       # make new file for each profile
  print("OUTPUT FILE: "+path+infile+"input.fa")  
  open(outputFasta, 'w').close()
  count=0
#  print(protGIs)
  for prot in list(set(protID)):              # For unique template protein
    print ("Searching Prot GI: "+prot)
    nucID=""
    protInfo = Entrez.efetch(db="protein",
                         id=prot,
                         rettype='gb',
                         retmode="xml")
    protrecord = Entrez.read(protInfo)                #pull GenBank record and read xml
    for sub in protrecord[0]["GBSeq_feature-table"]:  #parsing prot GenBank table...
      if sub["GBFeature_key"]=="source":
        sourcefeatures=sub["GBFeature_quals"]
        for annot in sourcefeatures:
          print(annot)
          if annot["GBQualifier_name"]=="organism":
            virusID=annot["GBQualifier_value"]
            print(virusID)
      if sub["GBFeature_key"]=="CDS":                 #find CDS...
        CDSfeatures=sub["GBFeature_quals"]
        #print(CDSfeatures)
        for feature in CDSfeatures:
          if feature["GBQualifier_name"]=="coded_by":
            templateFeature=feature["GBQualifier_value"]
            if(templateFeature.startswith("join")):
              print("Multiple Coding Regions")
              templateString =templateFeature.replace("join(","").replace(")","")
              #print(templateString)
              templateList = templateString.split(",")
              #print(templateList)
              nucID=[entry.split(":")[0] for entry in templateList]
              nucID=nucID[0]
              positions=[entry.split(":")[1] for entry in templateList]
              templateStart = [int(entry.split("..")[0]) for entry in positions]
              templateEnd = [int(entry.split("..")[1] )for entry in positions]
            else:
              print("Single Coding Region")
              splitTemplate=templateFeature.split(":")
              nucID=splitTemplate[0].strip()
              templateStart = int(splitTemplate[1].split("..")[0])
              templateEnd = int(splitTemplate[1].split("..")[1])

            #print(nucID+": "+str(templateStart)+str(templateEnd))

    #Find NT template for EVE
    print ("Matched in NCBI: "+ nucID)
    nuchandle = Entrez.efetch(db="nuccore",
                         id=nucID,
                         rettype='fasta',
                         retmode="xml")

    record=Entrez.read(nuchandle)
    #print(record[0]["TSeq_sequence"])
    sequenceStr=record[0]["TSeq_sequence"]
    if(isinstance( templateStart, int )):
      templateSeq=sequenceStr[(templateStart-1):(templateEnd-1)]
    else:
      templateSeq=[sequenceStr[(int(i)-1):(int(j)-1)] for (i,j) in (templateStart,templateEnd)]
      templateSeq="".join(templateSeq)
    #print(templateSeq)

    with open(outputFasta, 'a') as OF:
        OF.write(">"+virusID.replace(" ","_")+"."+str(entrycount)+"\n"+templateSeq+"\n")
    
    # for EVE in range(len(protID)):  
    #   if protID[EVE] == prot:
    #     count+=1
    #     print("Matching EVEs to "+prot+" from "+EVEid[EVE]+":"+str(EVEstart[EVE])+" to "+str(EVEend[EVE]))
    #     EVELab=">"+infile+"-"+str(count)

        #seq=str(eveDict[EVEid[EVE]])[(EVEstart[EVE]-1):(EVEend[EVE]-1)] #sequence region from contig
        
        # if vstrand[EVE] is False:
        #   print("Rev. Comp.'d 1")
        #   mystrand=Seq(seq, generic_dna)  
        #   seq=str(mystrand.reverse_complement())

        #print(estrand[EVE])
        #print(EVEregister[EVE])

        # if int(EVEregister[EVE])<0:
        #   print(seq[0:60])
        #   print("Rev. Comp.'d 2")
        #   mystrand=Seq(seq, generic_dna)  
        #   seq=str(mystrand.reverse_complement())
        #   print(seq[0:60])

        #with open(outputFasta, 'a') as OF:
        #  OF.write(EVELab+"\n"+seq+"\n")

    print("\n")
    entrycount+=1

eveDict=evedict(eveDir)
for infile in infiles:
  EVEid,protID,EVEstart,EVEend,virusstart,virusend,vstrand,estrand,EVEregister=extractEVEs(infile)
  EVElengths=[EVEend[i]-EVEstart[i] for i in range(len(EVEstart))]
  Virlengths=[virusend[i]-virusstart[i] for i in range(len(virusstart))]
  print(EVElengths)
  pp.scatter(EVElengths,Virlengths)
  #print(protID)
  idCDS(eveDict,protID,EVEstart,EVEend,EVEid,virusstart,virusend,vstrand,estrand,infile.replace(".result",""),EVEregister)

pp.show()

