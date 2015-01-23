#  #  #  #  #  #  #  #  #  #  #
# contigQuery.py
# Purpose: Uses BLAST hit tables (in .csv format with header) and genomic DNA sequences (in .fasta format)
# Usage: python contigQuery.py <.csv file> <genomic contig fasta file>
#  #  #  #  #  #  #  #  #  #  #

## IMPORTS ##

import sys #, os

## FUNCTIONS ##

def readfile(filename):             #Reads Contig .fasta file
    print "Reading contig file..."
    with open(filename) as file:
        print "...getting names"
        names=[part[0].split(' ')[0] for part in [entry.partition('\n') for entry in file.read().split('>')[1:]]]
    with open(filename) as file:
        print "...getting sequences"
        seqs=[part[2].replace("\n","") for part in [entry.partition('\n') for entry in file.read().split('>')[1:]]]
    print "...Done.\nTotal Contigs:", len(seqs)
    return zip(names,seqs) #return tuples of name and sequence


def parseQueries(queryFile):      #Pulls out relevant columns from BLAST table.
    with open(queryFile) as qF:
        names = [part[0] for part in [entry.split(',') for entry in qF.read().split('\r')[1:]]]
    with open(queryFile) as qF:
        starts = [part[5] for part in [entry.split(',') for entry in qF.read().split('\r')[1:]]]
    with open(queryFile) as qF:
        ends = [part[6] for part in [entry.split(',') for entry in qF.read().split('\r')[1:]]]
    with open(queryFile) as qF:
        vIDs= [part[12] for part in [entry.split(',') for entry in qF.read().split('\r')[1:]]]
    with open(queryFile) as qF:
        queryFrame = [part[11] for part in [entry.split(',') for entry in qF.read().split('\r')[1:]]]


    print "Total Queries:", len(names)
    queries = zip(names,starts, ends, vIDs, queryFrame)
    queries = [X for X in queries if X != ['','','']]
    return queries


def search(queries,contigs):                #Search for
    sequences=[]
    for query in queries:                   # Pull out sequences -- should be list comprehension
        for contig in contigs:
            if contig[0]==query[0]:
                coord1 = int(query[1])
                coord2 = int(query[2])
                sequence = contig[1][min(coord1,coord2)-1: max(coord1,coord2)]
                if int(query[4]) < 0:
                    sequence = revComp(sequence)
        sequences.append(sequence)
    return zip(queries,sequences)


def revComp(seq):
    compDict = {'A':'T','C':'G','G':'C','T':'A','a':'t','c':'g','g':'c','t':'a','n':'n','N':'N'}
    letters = list(seq)
    letters = [compDict[base] for base in letters]
    forComp = ''.join(letters)
    reverseComp=forComp[::-1]
    return reverseComp


def formatOutput(hits):
    outfile=open(multiFastaOutput,'w')
    for hit in hits:
        outfile.write(">{0}_{1}-{2}_{3}\n{4}\n".format(hit[0][0],hit[0][1],hit[0][2],hit[0][3].replace(" ","_").partition("(")[0],hit[1]))
    outfile.close()
    print "Output written to {0}".format(multiFastaOutput)

## GLOBALS ##

queryFile = sys.argv[1]
contigFile = sys.argv[2]
multiFastaOutput = "{0}{1}".format(queryFile.partition(".")[0],"_fasta.txt")
#other output formats...

## MAIN ##

contigs=readfile(contigFile)
queries=parseQueries(queryFile)
hits=search(queries,contigs)
formatOutput(hits)