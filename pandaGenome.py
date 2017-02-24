import pandas as pd
import numpy as np

##########################################

def generateInts(selectS,selectE): #generates the int numpy arrays for comparison.
	n=0
	x=0
	intList=[]
	#print(selectS)
	if len(selectS)>0:
		for i in selectS:
			if x==0:
				intList=np.arange(i,selectE[n]+1,1)
				x=1
			else:
				inter=np.arange(i,selectE[n]+1,1)
				intList=np.concatenate((intList,inter))
			n+=1

	Sites,Counts=np.unique(intList,return_counts=True) #super handy function for counting number of mapped elements.
	return(Sites,Counts)

##########################################
print("Reading Data...")

#READ in each file as a panda DF

conFile= "/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/Aag2_Contigs.fa.fai"
teFile="/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/Aag2_Contigs_TEs.bed"
eveFile="/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/Aag2_Contigs_EVEs_sorted.bed_withTaxonomy.txt"
piFile="/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/Aag2_piRNAs/I234_dsFluc-B_Aag2-PB.map_piRNA.csv"

conDF=pd.read_csv(conFile,sep="\t",header=None)
teDF=pd.read_csv(teFile,sep="\t",header=None)
eveDF=pd.read_csv(eveFile,sep="\t")
piDF=pd.read_csv(piFile,sep=",",header=None)

print("...done.")

##########################################

#Extract Start and End from each data type
piDF['contigPi']=[i[0] for i in piDF[2].str.split("|")]
startsPi=piDF[3].values
#print(startsPi)
endsPi=startsPi #Just a single Nt (5' start position)
contigPi=piDF['contigPi'].values

contigs=teDF[0].values
starts=teDF[1].values
ends=teDF[2].values

startEVE=np.array(eveDF['EVEstart'].values)
endEVE=np.array(eveDF['EVEend'].values)
contigEVE=eveDF['ContigEVE'].values

#Get set up for looping and annotating over contigs.
ucontigs=set(contigs)
outputDF=pd.DataFrame()
for c in ucontigs:
	print(c)
	contigInfo=conDF.loc[conDF[0]==c,1]

	selectSte=starts[contigs==c]
	selectEte=ends[contigs==c]
	
	selectSeve=startEVE[contigEVE==c]
	selectEeve=endEVE[contigEVE==c]

	selectSpi=np.array(startsPi[contigPi==c])
	selectEpi=endsPi[contigPi==c]

	#print("Generating Intervals...")
	eveSites,eveCounts=generateInts(selectSeve,selectEeve)	#generates the integer lists for comparison
	teSites,teCounts=generateInts(selectSte,selectEte)
	piSites,piCounts=generateInts(selectSpi,selectEpi)
	piStats=pd.DataFrame({"piSites":piSites,"piCounts":piCounts})
	
	intersection=np.intersect1d(eveSites,teSites)
	piintersection=np.intersect1d(eveSites,piStats['piSites'])
	
	selCount=piStats.loc[piStats['piSites'].isin(piintersection),'piCounts']
	#print(selCount)
	piReads=np.sum(piCounts)
	piReadsEVE=np.sum(selCount)
	outputDFlist=pd.DataFrame({'contig':c,'length':contigInfo,'teSites':len(teSites),'eveSites':len(eveSites),"eveTEoverlap":len(intersection),"piReads":piReads,"piEVEReads":piReadsEVE,"piSites":len(piSites),"piEVEoverlap":len(piintersection)})
	outputDFlist["teProp"]=outputDFlist['teSites']/outputDFlist['length']
	outputDFlist["eveTEProp"]=outputDFlist['eveTEoverlap']/outputDFlist['eveSites']
	outputDFlist["pi/eveProp"]=outputDFlist['piEVEoverlap']/outputDFlist['eveSites']
	outputDFlist["eve/piProp"]=outputDFlist['piEVEoverlap']/outputDFlist['piSites']
	outputDFlist["piReads/EVEProp"]=outputDFlist['piEVEReads']/outputDFlist['piReads']
	outputDF=pd.concat([outputDF,outputDFlist])
	print(len(outputDF))
	outputDF.to_csv("EVE_TEpi_Analysis_byContig.csv")

#EVE analysis 
allFrames=pd.DataFrame()
allStats=pd.DataFrame()
piTotal=len(piDF)
for i in eveDF.iterrows():
    contig=i[1]['ContigEVE']
    family=i[1]['family']
    species=i[1]['species']
    ID=i[1]['EVEdescription']
    piSel=piDF.loc[piDF.contigPi==contig,3]
    interval=np.arange(i[1]["EVEstart"],i[1]["EVEend"],1)
    piSel,piIndex,piCounts=np.unique(piSel,return_index=True,return_counts=True)
    piFrame=pd.DataFrame({"piSel":piSel,"piIndex":piIndex,"piCounts":piCounts,"piID":ID,"contig":contig})
    allFrames=pd.concat([allFrames,piFrame])
    print(i[0])

    piEVEinter=np.intersect1d(piFrame['piSel'],interval)
    piEVEreads=piFrame.loc[piFrame['piSel'].isin(piEVEinter).values,'piCounts']
#    print(piEVEreads)
    piEVEoverlap=len(piEVEinter)

    piStats=pd.DataFrame({"ID":ID,"piPos":[len(piSel)],"piReads":[sum(piCounts)],"piEVEreads":[np.sum(piEVEreads)],"piEVEPos":[piEVEoverlap],"family":[family],"species":[species],"EVEPos":[len(interval)],"contig":contig,"piTotal":[piTotal]})
    allStats=pd.concat([allStats,piStats])
    
allStats["piHit"]=allStats['piEVEreads']>0
allStats["piReadsPerPos"]=allStats['piEVEreads']/allStats['EVEPos']
allStats["piEVEProp"]=allStats['piEVEreads']/allStats['piTotal']
allStats["EVEpiCover"]=allStats['piEVEPos']/allStats['EVEPos']

allStats['piPerTESite']=allStats['piTotal']/TEtotal
allStats['piPerSite']=allStats['piTotal']/1723952533
allStats['propGenome']=allStats['EVEPos']/1723952533
allStats['enrichGenome']=allStats['piEVEProp']/allStats['propGenome']
allStats['enrichGenome']=allStats['piEVEProp']/allStats['propGenome']

teDF['length']=(teDF[1]-teDF[2]).abs()
allStats.to_csv("allStats_byEVE.csv")