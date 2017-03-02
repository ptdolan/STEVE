#
#   piCluster_dotplot.py
#
#   Adapted from K.M. Dalton's dotsplot.py
#   https://github.com/kmdalton/dotsplot
#

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
hfont = {'fontname':'Helvetica'}

AagContigs = "/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/Aag2_assembly/Aag2_Contigs.fa"
LVPContigs = "/Users/ptdolan/Research/EVEsAndpiRNA/Genomes/Aedes_aegypti-LVP-Contigs.fasta"
piClusterFile="/Users/ptdolan/Research/EVEsAndpiRNA/Frozen_Data/proTRAC_piRNAs.map_2016y11m1d22h53m3s/Aag2_piClusters.bed"
emem="/Users/ptdolan/Research/EVEsAndpiRNA/e-mem_1000.out" # from EForm.py

window=5000
def contigDict(contigfasta):
	contigdict = SeqIO.to_dict(SeqIO.parse(contigfasta, "fasta"))
	return(contigdict)

def pullSeqWindow(p):
    pcontig=p[1][0]
    piStart=min(p[1][1],p[1][2])
    piEnd=max(p[1][1],p[1][2])
    print("pi positions"+str(piStart)+" : "+str(piEnd))
    aligns=EMEM.loc[EMEM[0]==pcontig,:]
    sortAlign=aligns.sort(1)
    #print(p)
    if sortAlign.empty!=True:
        #print(pcontig[0])
        upperneeb=sortAlign.loc[sortAlign[1]>(piEnd+window),:].iloc[1]
        lowerneeb=sortAlign.loc[sortAlign[1]<(piStart-window),:].iloc[-1]

        print(lowerneeb[1])
        print(upperneeb[1])

        if(upperneeb[4]==lowerneeb[4]):
            AAGseq=str(cD_A[pcontig].seq)[lowerneeb[1]:upperneeb[1]]
            LVPseq=str(cD_L[upperneeb[4]].seq)[lowerneeb[2]:upperneeb[2]]
            adj=piStart-lowerneeb[1]
            print(adj)
            print("Length of Window: "+"LVP: "+str(len(LVPseq))+" AAG: "+str(len(AAGseq)))
        else:
            LVPseq=""
            AAGseq=""
            adj=0
    return(AAGseq,LVPseq,(piStart,piEnd,adj,lowerneeb[4]))

def dotplot(seqA,seqB,window):
    strarrayA=pd.dataframe([seqA[i:(i+window)] for i in range(len(seqA)-window)]*(len(seqB)-window))
    strarrayB=pd.dataframe([seqB[i:(i+window)] for i in range(len(seqB)-window)]*(len(seqA)-window)).transpose
    arr=strarrayA==strarrayB
    return(pd.dataframe(array,dtype="int"))

def encode_sequence(sequence, k, **kwargs):
    """
    encode_sequence(sequence, k, **kwargs)

    Parameters
    ----------
    sequence : str
        Sequences you wish to encode
    k : int
        Word length
    Returns
    -------
    encoded_seq : numpy.ndarray
        array of integers describing the encoded sequence
    kmer_encoding : dict
        the kmers corresponding to each int in encoded_seq
    """
    kmer_encoding = kwargs.get('kmer_encoding', {})
    n = 0
    encoded = []
    for i in range(len(sequence) - k):
        subsequence = sequence[i:i+k]
        if subsequence not in kmer_encoding:
            kmer_encoding[subsequence] = n
            n += 1
        encoded.append(kmer_encoding[subsequence])

    return np.array(encoded), {v:k for k,v in kmer_encoding.iteritems()}

def revcomp(dna):
    fordna = Seq(dna, generic_dna)
    revdna = fordna.reverse_complement()
    return(str(revdna))

def get_matches(eseqA, eseqB):
    """
    encode_sequence(sequence, k, **kwargs)

    Parameters
    ----------
    encoded_sequence : numpy.ndarray
        array of integers describing the sequence. generate with dotsplot.encode_sequence

    Returns
    -------
    X : numpy.ndarray
        array of matches, first coordinate
    Y : numpy.ndarray
        array of matches, second coordinate
    """
    X,Y = [],[]
    for x,kmer in enumerate(eseqA):
        #print np.where(encoded_seq == kmer)
        y = np.where(eseqB == kmer)[0]
        X = np.concatenate((X, [x]*len(y)))
        Y = np.concatenate((Y, y))
        #X = np.concatenate((X, [x]*len(y)))
    return X,Y

#pp.savefig("piCluster"+str(p[0])+".png", dpi=300, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None,transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)

#######

if __name__=="__main__":
    cD_A=contigDict(AagContigs)
    cD_L=contigDict(LVPContigs)

    EMEM=pd.read_csv(emem,header=None,sep="\t")

    pi=pd.read_csv(piClusterFile,sep="\t",header=None)
    pi.head()
    k = 30
    for p in pi.iterrows():
        print(p)
        LVPseq=""
        AAGseq=""
        try:
            AAGseq,LVPseq,piCoord=pullSeqWindow(p)
        except: pass
        if AAGseq is not "" and LVPseq is not "":
        # seqlen = 1000
        # seqA = ''.join(np.random.choice(['A', 'T', 'C', 'G'], seqlen))
        # seqB = ''.join(np.random.choice(['A', 'T', 'C', 'G'], 2*seqlen))

            encoded_seqA, kmer_encodingA = encode_sequence(AAGseq, k)
            encoded_seqB, kmer_encodingB = encode_sequence(LVPseq, k)

            X,Y = get_matches(encoded_seqA, encoded_seqB)
            plt.rc('text', usetex=True )
            plt.figure(figsize=(5,5))
            plt.rc('font', family='Helvetica')
            plt.scatter(X, Y,s=0.1,c='black',marker=".")
            plt.xlabel('Aag2 sequence',**hfont)
            plt.ylabel('LVP sequence',**hfont)
            #print(piCoord[0]-piCoord[2])
            plt.axvline(piCoord[2],color='grey')
            plt.axvline(piCoord[2]+piCoord[1]-piCoord[0],color='grey')
            plt.suptitle(p[1][0]+" "+str(p[1][1])+":"+str(p[1][2])+"\n"+piCoord[3])
            plt.savefig(p[1][0]+str(p[0])+"_ForWin-"+str(window)+".pdf")
            plt.clf()

            A,B = get_matches(encoded_seqA, encoded_seqA)
            len(encoded_seqA)
            plt.figure(figsize=(5,5))
            plt.rc('font', family='Helvetica')
            plt.scatter(A, B,s=0.1,c='black',marker=".")
            plt.xlabel('Aag2 sequence',**hfont)
            plt.ylabel('Aag2 sequence',**hfont)
            #print(piCoord[0]-piCoord[2])
            plt.axvline(piCoord[2],color='grey')
            plt.axvline(piCoord[2]+piCoord[1]-piCoord[0],color='grey')
            plt.axvspan(piCoord[2],piCoord[2]+piCoord[1]-piCoord[0],color='grey', alpha=0.25)
            plt.suptitle(p[1][0]+" "+str(p[1][1])+":"+str(p[1][2]))
            plt.savefig(p[1][0]+str(p[0])+"Aag2Self_ForWin-"+str(window)+".pdf")
            plt.clf()
