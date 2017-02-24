# e-mem reformatting script
# EForm.py <E-MEM file> <.tsv output >
# 
#
import sys,csv
import operator
import matplotlib.pyplot as pp
import matplotlib as mpl

infile=sys.argv[1]
output=[]

with open(infile,'r') as IF: 
	info=[line.strip().split() for line in IF]
#print(info)

outfile=infile+".out"
of = open(outfile,'w') 
contigCount=dict()
for line in info:
	if(line[0]==">"):
#		print("Contig")
		sourceContig=line[1]
		rev="For"
		if(line[2] == "Reverse"):
			#print(sourceContig)
			rev="Rev"
	else:
		line.append(sourceContig)
		line.append(rev)
		for x in range(len(line)):
			of.write(line[x])
			if x <len(line):
				of.write("\t")
		of.write("\n")
		keystr=sourceContig+"_"+line[0]
		if keystr in contigCount.keys():
			contigCount[keystr]=contigCount[keystr]+1
		else: contigCount[keystr]=1
of.close()

sortedPairs = sorted(contigCount.items(), key=operator.itemgetter(1))

rankeddegree=[ pair[1] for pair in sortedPairs ]

pp.plot(range(len(rankeddegree)), rankeddegree)
pp.yscale("log")
pp.savefig(infile+".pdf")


