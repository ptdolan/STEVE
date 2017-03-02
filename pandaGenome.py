import pandas as pd
import numpy as np
import sys

"""
# ## ## ## ## # ## ## ## ## ## ## ## ## #
# ## ## ## ##  FUNCTIONS  # ## ## ## ## #
# ## ## ## ## # ## ## ## ## ## ## ## ## #
"""

def readFiles(conFile,teFile,eveFile,piFile,clustFile,clustTable):
    print("\nReading Data...")
    conDF=pd.read_csv(conFile,sep="\t",header=None)

    teDF=pd.read_csv(teFile,sep="\t",header=None)
    eveDF=pd.read_csv(eveFile,sep="\t")
    piDF=pd.read_csv(piFile,sep=",",header=None)

    piDF.columns=['piSeq','strand','contigInfo','piPos','complement']
    piDF['contigPi']=[i[0] for i in piDF.contigInfo.str.split("|")]
    clustDF=pd.read_csv(clustFile,sep="\t",header=None)
    clustInfo=pd.read_csv(clustTable,sep="\t",header=None,skiprows=59) #this is a real shit format.
    print("done.")

    return(conDF,teDF,eveDF,piDF,clustDF,clustInfo)

######
# Progress bar
# Adapted from: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
######

def progBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 40, fill = '|'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    sys.stdout.write('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix))
    sys.stdout.flush()

    # Print New Line on Complete
    if iteration == total:
        print()

################
# generates dataframe of piRNA coordinates, counts, and index.
################

def piFrame(piSel,contig):
    piSel,piCounts=np.unique(piSel,return_counts=True)
    piF=pd.DataFrame({"piSel":piSel,"piCounts":piCounts,"contig":contig})
    return(piF)

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

################
# By Contig Analysis of Eves and piRNAs
################
def contiganalysis(piDF,eveDF,teDF):
    print("\nBy-Contig Analysis...")

#piRNAs
    startsPi=piDF[3].values
    endsPi=startsPi #Just a single Nt (5' start position)
    contigPi=piDF['contigPi'].values

#TEs
    contigs=teDF[0].values
    starts=teDF[1].values
    ends=teDF[2].values

#EVEs
    startEVE=eveDF['EVEstart'].values
    endEVE=eveDF['EVEend'].values
    contigEVE=eveDF['ContigEVE'].values

    ucontigs=np.unique(contigs,return_counts=True)
    outputDF=pd.DataFrame()
    l=len(ucontigs)
    for C in ucontigs.iteritems():
        progBar(C[0],l)
        contigInfo=conDF.loc[conDF[0]==C,1]
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
        return(outputDF)

#######################
# By EVE analysis of piRNAs
#######################
def EVEanalysis(piDF,eveDF,teDF,TEtotal):
    print("\nby-EVE analysis...")
    allFrames=pd.DataFrame()
    allStats=pd.DataFrame()
    piTotal=len(piDF)
    l=len(eveDF)
    for i in eveDF.iterrows():
        progBar(i[0],l)
        contig=i[1]['ContigEVE']
        family=i[1]['family']
        species=i[1]['species']
        ID=i[1]['EVEdescription']
        piSel=piDF.loc[piDF.contigPi==contig,'piPos']
        interval=np.arange(i[1]["EVEstart"],i[1]["EVEend"]+1,1) #EVE interval positions array
        piF=piFrame(piSel,contig)
        allFrames=pd.concat([allFrames,piF])
        piEVEinter=np.intersect1d(piF['piSel'],interval)
        countVector=piF['piCounts'].loc[piF['piSel'].isin(piEVEinter).values]
        piEVEreads=countVector
        piEVEoverlap=len(piEVEinter)
        piStats=pd.DataFrame({"ID":ID,"piPos":[len(piSel)],"piReads":[piF.piCounts.sum()],"piEVEreads":[np.sum(piEVEreads)],"piEVEPos":[piEVEoverlap],"family":[family],"species":[species],"EVEPos":[len(interval)],"contig":contig,"piTotal":[piTotal]})
        allStats=pd.concat([allStats,piStats])

    allStats["piReadsPerPos"]=allStats['piEVEreads']/allStats['EVEPos']
    allStats["piEVEProp"]=allStats['piEVEreads']/allStats['piTotal']
    allStats["EVEpiCover"]=allStats['piEVEPos']/allStats['EVEPos']

    allStats['piPerTESite']=allStats['piTotal']/TEtotal
    allStats['piPerSite']=allStats['piTotal']/1723952533
    allStats['propGenome']=allStats['EVEPos']/1723952533
    allStats['enrichGenome']=allStats['piEVEProp']/allStats['propGenome']
    allStats['enrichGenome']=allStats['piEVEProp']/allStats['propGenome']

    teDF['length']=(teDF[1]-teDF[2]).abs()
    return(allStats)

############################
# by piCluster ANALYSIS OF EVEs and piRNAs
############################
def piClustanalysis(piDF,clustDF,eveDF,TEtotal):
    print("\nby-Cluster analysis...")
    allFrames=pd.DataFrame()
    allStatsClust=pd.DataFrame()
    piTotal=len(piDF)
    l=len(clustDF)
    for i in clustDF.iterrows():
        progBar(i[0],l)
        contig=i[1][0]
        clusterS=i[1][1]
        clusterE=i[1][2]
        clustInterval=np.arange(clusterS,clusterE+1,1)
        piSel=piDF.loc[piDF.contigPi==contig,'piPos'] #piSel: piRNAs mapping to a given contig
#        pSel.columns=['','','','']
        piF=piFrame(piSel,contig)
        allFrames=pd.concat([allFrames,piF])
        piClustinter=np.intersect1d(piF['piSel'],clustInterval)
        piClustreads=piF.loc[piF['piSel'].isin(piClustinter).values,'piCounts']
        eveSel=eveDF.loc[eveDF.ContigEVE==contig,:]
        EVEhit=0
        EVEtotal=eveSel.EVEdescription.count()
        EVEClustOlen=0
        ClustEVEpiOlen=0
        ClustEVEpiOreads=0
        EVEClustoverlapreads=0
        EVEf=[]
        if EVEtotal>0:
            EVEfamily=eveSel.family
            EVEspecies=eveSel.species
            EVEID=eveSel.EVEdescription
            for EVE in eveSel.iterrows():
                EVEinterval=np.arange(EVE[1]["EVEstart"],EVE[1]["EVEend"]+1,1)
                EVEClustoverlap=np.intersect1d(EVEinterval,clustInterval)   #cluster-EVE overlap
                EVEClustoverlapreads+=piF.loc[piF['piSel'].isin(EVEClustoverlap).values,'piCounts'].sum()
                EVEClustOlen+=len(EVEClustoverlap)
                ClustEVEpioverlap=np.intersect1d(EVEClustoverlap,piClustinter)
                ClustEVEpiOlen+=len(ClustEVEpioverlap)
                ClustEVEpiOreads+=piF.loc[piF['piSel'].isin(np.intersect1d(EVEClustoverlap,piClustinter)).values,'piCounts'].sum()
                if len(EVEClustoverlap)>0:
                    F=EVE[1]['family']
                    EVEf.append(str(F))
                    EVEhit+=1
                    #print(EVEhit)
        eveF=len(set(EVEf))
        clustEVEcover=EVEClustOlen/i[1][5]
        index=i[0]
        #print(index)
    allStatsClust=pd.concat([allStatsClust,pd.DataFrame({"piTotal":[piTotal],"EVEhits":[EVEhit],"EVEfamilies":[eveF],"EVE-ClusterOverlap":[EVEClustOlen],"ClusterEVEpiOverlap":[ClustEVEpiOlen],"ClusterEVEpiReads":[ClustEVEpiOreads],"clusterEVEcoverage":[clustEVEcover],"piClust_reads":[piClustreads.sum()],"piClust_piSites":[len(piClustinter)]},index=[i[0]])])
    return(allStatsClust)
