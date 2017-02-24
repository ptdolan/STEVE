from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pp



AagContigs = "/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/Aag2_Contigs.fa"
AagEVEs = "/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/Flavivirae_EVEs/Aag2_Contigs_EVEs.fasta"
def contigDict(contigfasta):
	with open(contigfasta,'r') as cfa:
		contigs=[line.strip() for line in cfa]
		names=[entry.split(" ")[0][1:] for entry in contigs if entry[0]==">"]
		sequences=[entry.split(" ")[0] for entry in contigs if entry[0]!=">"]
		contigdict=dict(zip(names,sequences))
		return(contigdict)

contigDict=contigDict(AagContigs)

with open(AagEVEs,'r') as EVEs:
	EVEslist=[line for line in EVEs]
	

