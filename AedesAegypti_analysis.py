import pandaGenome as PG
##########################################
#	Define input file paths.
##########################################
conFile= "/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/Aag2_Contigs.fa.fai"
teFile="/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/Aag2_Contigs_TEs.bed"
eveFile="/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/Aag2_Contigs_EVEs_sorted.bed_withTaxonomy.txt"
piFile="/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/Aag2_piRNAs/I234_dsFluc-B_Aag2-PB.map_piRNA.csv"
clustFile="/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/proTRAC_piRNAs.map_2016y11m1d22h53m3s/Aag2_piClusters.bed"
clustTable="/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/proTRAC_piRNAs.map_2016y11m1d22h53m3s/results.table"

#	READ in each file as a panda DF
conDF,teDF,eveDF,piDF,clustDF,clustInfo=PG.readFiles(conFile,teFile,eveFile,piFile,clustFile,clustTable)
print("...done.")
#   Individual comparison scripts. Generate a pdDF and write to CSV

# contigOut=PG.contiganalysis(piDF,eveDF,teDF,TEtotal)
# contigOut.to_csv("EVE_TE_pi_Analysis_byContig.csv")
# TEtotal = contigOut.teSites.sum()
TEtotal=926471910
print("Total TEs:"+str(TEtotal))
#clustOut=PG.piClustanalysis(piDF,clustDF,eveDF,TEtotal=926471910)
#clustOut.to_csv("EVE_TE_pi_Analysis_byClust.csv")

eveOut=PG.EVEanalysis(piDF,eveDF,teDF,TEtotal=926471910)
eveOut.to_csv("EVE_TE_pi_Analysis_byEVE.csv")
