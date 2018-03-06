import pandas as pd
import numpy as np

##########################################
def generateInts(selectS,selectE):
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

	Sites,Counts=np.unique(intList,return_counts=True)
	return(Sites,Counts)

##########################################
print("Reading Data...")
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

piDF['contigPi']=[i[0] for i in piDF[2].str.split("|")]
startsPi=piDF[3].values
#print(startsPi)
endsPi=startsPi
contigPi=piDF['contigPi'].values

contigs=teDF[0].values
starts=teDF[1].values
ends=teDF[2].values

ucontigs=set(contigs)
startEVE=np.array(eveDF['EVEstart'].values)
endEVE=np.array(eveDF['EVEend'].values)
contigEVE=eveDF['ContigEVE'].values

outputDF=pd.DataFrame()
# for c in ucontigs:
# 	print(c)
# 	contigInfo=conDF.loc[conDF[0]==c,1]

# 	selectSte=starts[contigs==c]
# 	selectEte=ends[contigs==c]
	
# 	selectSeve=startEVE[contigEVE==c]
# 	selectEeve=endEVE[contigEVE==c]

# 	selectSpi=np.array(startsPi[contigPi==c])
# 	selectEpi=endsPi[contigPi==c]

# 	#print(selectSpi)
# 	#print("Generating Intervals...")
# 	eveSites,eveCounts=generateInts(selectSeve,selectEeve)
# 	teSites,teCounts=generateInts(selectSte,selectEte)
# 	piSites,piCounts=generateInts(selectSpi,selectEpi)
# 	piStats=pd.DataFrame({"piSites":piSites,"piCounts":piCounts})
	
# 	intersection=np.intersect1d(eveSites,teSites)
# 	piintersection=np.intersect1d(eveSites,piStats['piSites'])
	
# 	selCount=piStats.loc[piStats['piSites'].isin(piintersection),'piCounts']
# 	#print(selCount)
# 	piReads=np.sum(piCounts)
# 	piReadsEVE=np.sum(selCount)
# 	outputDFlist=pd.DataFrame({'contig':c,'length':contigInfo,'teSites':len(teSites),'eveSites':len(eveSites),"eveTEoverlap":len(intersection),"piReads":piReads,"piEVEReads":piReadsEVE,"piSites":len(piSites),"piEVEoverlap":len(piintersection)})
# 	outputDFlist["teProp"]=outputDFlist['teSites']/outputDFlist['length']
# 	outputDFlist["eveTEProp"]=outputDFlist['eveTEoverlap']/outputDFlist['eveSites']
# 	outputDFlist["pi/eveProp"]=outputDFlist['piEVEoverlap']/outputDFlist['eveSites']
# 	outputDFlist["eve/piProp"]=outputDFlist['piEVEoverlap']/outputDFlist['piSites']
# 	outputDFlist["piReads/EVEProp"]=outputDFlist['piEVEReads']/outputDFlist['piReads']
# 	outputDF=pd.concat([outputDF,outputDFlist])
# 	print(len(outputDF))
# outputDF.to_csv("EVE_TEpi_Analysis.txt")

print(eveDF.datadtypes)


