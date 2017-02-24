library(ggplot2)
library(limma)
library(reshape2)
library(pbapply)
Crange=1:3

inDir<-"~/Research/EVEsAndpiRNA/"
setwd(inDir)

load("Frozen_Data/Contigs.RData")
load("Alignment.RData")

edgeCounts=NULL

allalign$LVPAAG=F
# filter for ambiguity in LVP->Aag2
# must be mapped thrice to a given contig target in Aag2. Will be useful for visualization, I think.
allalign$jointname=as.factor(paste(allalign$Contig1,allalign$Contig2))

totalalignedc1<-by(allalign$start.c1,allalign$jointname,FUN = function(X){abs(range(X)[1]-range(X)[2])})
totalalignedc2<-by(allalign$start.c2,allalign$jointname,FUN = function(X){abs(range(X)[1]-range(X)[2])})
names<-names(totalaligned)
lengthsc1<-matrix(totalalignedc1)
lengthsc2<-matrix(totalalignedc2)
tol<-(lengthsc1/lengthsc2)>0.8&(lengthsc1/lengthsc2)<1.2

lengthDF<-data.frame(names,lengthsc1,tol)

ggplot(lengthDF)+geom_hex(aes(lengths))+scale_x_log10()+scale_y_log10()

filtLengths<-lengthDF[lengthDF$lengthsc1>10000&tol==TRUE,]
ggplot(filtLengths)+geom_point(aes(c2Length,lengths))+scale_x_log10()+scale_y_log10()

filtered<-allalign[allalign$jointname%in%filtLengths$names,]
filtered<-allalign[allalign$jointname%in%names(table(allalign$jointname)[table(allalign$jointname)>15]),]
filtered$window
filtered$length<-filtered$length*2*(as.integer(filtered$dir)-1.5)

for(contig in unique(filtered$Contig1[allContigs$ID=="Aag2_Contigs"])){
  print(contig)
  matches=allalign[allalign$Contig1==contig,]
  edgeCountDF=data.frame(matches)
  edgeCountDF=edgeCountDF[order(edgeCountDF$start.c2),]
  edgeCountDF=edgeCountDF[order(edgeCountDF$start.c1),]
  filteredCounts=table(edgeCountDF$Contig2)[table(edgeCountDF$Contig2)>1]
  edgeCountDF<-edgeCountDF[edgeCountDF$Contig2%in%names(filteredCounts),]
  edgeCountDF=data.frame(matches)
  edgeCountDF=edgeCountDF[order(edgeCountDF$start.c2),]
  edgeCountDF=edgeCountDF[order(edgeCountDF$start.c1),]
  edgeCountDF$adjPosC2=0
  #lengths=data.frame()
  
  G<-ggplot(edgeCountDF)+
    ylab("LVP Contig Position")+
    xlab(paste("Aag2 Contig:",contig))+
    geom_point(aes(x=start.c1,y=start.c2,color=reorder(Contig2,start.c1)),cex=0.2)+ 
    coord_fixed(ratio = 1)+theme_bw()+theme(legend.position="none")
  
  Glin<-ggplot(edgeCountDF)+
    ylab("LVP Contig Position")+
    xlab(paste("Aag2 Contig:",contig))+
    geom_segment(aes(x=start.c1,xend=start.c1+length,y=reorder(Contig2,start.c1,mean),yend=reorder(Contig2,start.c1),color=reorder(Contig2,start.c1)))+ 
    theme_bw()+theme(legend.position="none")
  
  ggsave(width = 5,height=3.5,paste(contig,"alignment.pdf",sep=""),G)
  ggsave(width = 5,height=3.5,paste(contig,"alignmentstack.pdf",sep=""),Glin)
}


out2<-by(allalign$Contig2,allalign$Contig1,FUN=function(X){length(unique(X))})
out1<-by(allalign$Contig1,allalign$Contig2,FUN=function(X){length(unique(X))})


par(mar=c(10,6,8,6))
boxplot(pch=16,alpha=0.2,cex=0.3,list(as.vector(out1),as.vector(out2)),jitter=0.2,log = "y",las=2,ylab="Number of Aligned Contigs",names = c("PacBio -> Sanger","Sanger -> PacBio"),col = cos[2:1])

