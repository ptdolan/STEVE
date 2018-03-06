#TE EVE
library(limma)
library(data.table)
library(ggplot2)
library(reshape2)
library(viridisLite)
library(viridis)
library("beeswarm")

#byEVE
EVEpiInfo<-read.csv("/Users/ptdolan/Research/EVEsAndpiRNA/FrozenData_4-27-17/Scripts/TE-EVE-piRNAexpression/EVE_TE_pi_Analysis_byEVE.csv")
plot(EVEpiInfo$negreads+EVEpiInfo$posreads,EVEpiInfo$piEVEreads)

#bypiCluster
piClusterInfo<-read.csv("/Users/ptdolan/Research/EVEsAndpiRNA/FrozenData_4-27-17/Scripts/TE-EVE-piRNAexpression/EVE_TE_pi_Analysis_byClust.csv")
piClusterTable<-read.delim("/Users/ptdolan/Research/EVEsAndpiRNA/FrozenData_4-27-17/proTRAC_piRNAs.map_2016y11m1d22h53m3s/results.clean.txt",sep = "\t",header = F)
piClusterInfo$Length=piClusterInfo$End-piClusterInfo$Start
piClusterInfo$density<-piClusterInfo$piClust_reads/piClusterInfo$Length

brokenTableNames<-apply(piClusterTable,MARGIN = 2,FUN = function(X){unique(strsplit2(X,": ")[,1])})
brokenTable<-as.data.frame(apply(piClusterTable,MARGIN = 2,FUN = function(X){strsplit2(X,": ")[,-1]}))
TFinfo<-apply(brokenTable[,-c(1:11)],MARGIN = 1,FUN = function(X){paste(X,sep="",collapse = "")})
fixedTable<-data.frame(brokenTable[1:11],TFinfo)
fixedTable$V1<-paste("Cluster",1:nrow(fixedTable))
fixedTable$X<-0:(nrow(fixedTable)-1)
Names<-append(unlist(brokenTableNames[-1])[-13],values = c("ClusterID","X"))
colnames(fixedTable)<-Names
output<-apply(MARGIN = 1,fixedTable,FUN = function(X){data.frame(cbind(TFs=rownames(table(strsplit2(X[12]," \\(.{1,30}\\)[ ]{0,1}" ))),count=table(strsplit2(X[12]," \\(.{1,30}\\)[ ]{0,1}" )),contig=X[1],ClusterID=X[13]))})
out<-plyr::ldply(output)
mout<-melt(out,measure.vars = "count")
castTFcounts<-dcast(mout,formula = ClusterID+contig~TFs,fill = 0)
piClusterTable<-merge(fixedTable,castTFcounts,by = "ClusterID")
str(piClusterTable)
m<-melt(piClusterTable,measure.vars = 16:23,variable.name = "TF",value.name = "count")

m$count<-as.numeric(m$count)
colnames(m)
m$x1<-as.numeric(strsplit2(m$Coordinates,"-")[,1])
m$x2<-as.numeric(strsplit2(m$Coordinates,"-")[,2])
m$`Size [bp]`<-as.integer(as.character(m$`Size [bp]`))
m$`Hits (normalized)`<-as.numeric(as.character(m$`Hits (normalized)`))
m$`Hits (absolute)`<-as.numeric(as.character(m$`Hits (absolute)`))
m$`Hits (normalized) per kb`<-as.numeric(as.character(m$`Hits (normalized) per kb`))
m$`Normalized hits with 1T`<-as.numeric(stringr::str_replace(as.character(m$`Normalized hits with 1T`),"\\%",""))
m$`Normalized hits with 10A`<-as.numeric(stringr::str_replace(as.character(m$`Normalized hits with 10A`),"\\%",""))
m$`Normalized hits 24-32 nt`<-as.numeric(stringr::str_replace(as.character(m$`Normalized hits 24-32 nt`),"\\%",""))
m$`Normalized hits on the main strand(s)`<-as.numeric(stringr::str_replace(as.character(m$`Normalized hits on the main strand(s)`),"\\%",""))

##MERGED DFs
mergeCluster<-merge(x = m ,by = "X",y = piClusterInfo,all = T)

logistic_smooth<-function(X) (geom_smooth(method='glm', method.args = list(family = "binomial")))

#pdf("TFandtranscriptionData.pdf",width = 10,height = 8)
mCgg<-ggplot(mergeCluster[mergeCluster$TF!="NA",])+scale_alpha(0.6)
mCgg+geom_point(aes(count,`Hits (absolute)`,color=EVEhits>0))+geom_smooth(level=.99,method = 'lm',aes(count,`Hits (absolute)`))+facet_wrap(~TF,scales = 'free')
mCgg+geom_point(aes(count,`Hits (normalized) per kb`,color=clusterEVEcoverage>0))+geom_smooth(method = 'lm',aes(count,`Hits (normalized) per kb`))+facet_wrap(~TF,scales = 'free')+coord_cartesian(ylim =c(0,700))


mCgg+geom_violin(aes(clusterEVEcoverage>0,`Hits (normalized) per kb`,fill=clusterEVEcoverage>0))+
  geom_boxplot(aes(clusterEVEcoverage>0,`Hits (normalized) per kb`),outlier.alpha = 1,width = 0.3)+
  coord_cartesian(ylim =c(10,2500))+
  scale_y_log10()

mCgg+geom_boxplot(aes(clusterEVEcoverage>0,`Hits (normalized) per kb`,color=clusterEVEcoverage>0))+coord_cartesian(ylim =c(10,2500))+scale_y_log10()

wilcox.test(mergeCluster$piClust_reads~mergeCluster$EVEhits>0)
wilcox.test(mergeCluster$`Hits (normalized) per kb`~mergeCluster$EVEhits>0)
wilcox.test(mergeCluster$`Normalized hits with 1T`~mergeCluster$EVEhits>0)
wilcox.test(mergeCluster$`Normalized hits with 10A`~mergeCluster$EVEhits>0)
wilcox.test(mergeCluster$density~mergeCluster$EVEhits>0)
wilcox.test(mergeCluster$count~mergeCluster$TF)

# mCgg+geom_point(aes(count,`Normalized hits on the main strand(s)`,color=clusterEVEcoverage>0))+facet_wrap(~TF,scales = 'free')
# mCgg+geom_boxplot(aes(clusterEVEcoverage>0,`Normalized hits on the main strand(s)`,color=clusterEVEcoverage>0))+facet_wrap(~TF,scales = 'free')

# mCgg+geom_point(aes(`Hits (normalized) per kb`,`Normalized hits 24-32 nt`,color=EVE.ClusterOverlap>0))+scale_x_log10()
# mCgg+geom_point(aes(`Hits (normalized) per kb`,`Normalized hits with 1T`,color=EVE.ClusterOverlap>0))+scale_y_log10()+scale_x_log10()
# mCgg+geom_density(aes(`Normalized hits with 10A`,fill=EVE.ClusterOverlap>0))
# mCgg+geom_density(aes(`Normalized hits with 1T`,fill=EVE.ClusterOverlap>0))

ggplot(mergeCluster,aes(EVE.ClusterOverlap>0,count))+
  geom_violin(scale = 'width',aes(fill=EVE.ClusterOverlap>0))+
  geom_boxplot(width = 0.3)+scale_y_sqrt()+
  facet_grid(~TF,margins=TRUE,scales="free")

mCgg+geom_point(aes(rank(ties.method = 'last',density),density))+scale_y_log10()

dev.off()

TEinfo<-read.csv("~/Research/EVEsAndpiRNA/FrozenData_4-27-17/EVE_TE_pi_Analysis_byContig.csv")

piClusterInfo$X<-with(piClusterInfo,reorder(X,piClust_reads,median))

ggplot(piClusterInfo)+theme_classic()+
  geom_point(aes(X,piClust_reads/piTotal,fill=EVE.ClusterOverlap>0,cex=EVEhits),pch = 21,stat="identity")+
  xlab("piCluster")+
  ylab("Proportion of piRNA reads")+
  scale_y_log10()

ggplot(piClusterInfo)+theme_classic()+
  geom_point(aes(X,density,fill=EVE.ClusterOverlap>0,cex=EVEhits),pch = 21,stat="identity")+
  scale_y_log10()

ggplot(piClusterInfo)+theme_classic()+
  geom_violin(aes(X,density,fill=EVE.ClusterOverlap>0,cex=EVEhits),pch = 21,stat="identity")+
  xlab("piCluster")+
  ylab("Reads per piCluster bp")+
  scale_y_log10()

ggplot(piClusterInfo)+
  theme_classic()+
  geom_violin(alpha = 1,aes(EVE.ClusterOverlap>0,piClust_reads/piTotal,fill=`EVE.ClusterOverlap`>0))+
  geom_boxplot(aes(EVE.ClusterOverlap>0,piClust_reads/piTotal),width=0.3)+
  scale_y_log10()

ggplot(piClusterInfo)+
  theme_classic()+
  geom_violin(alpha = 1,aes(EVE.ClusterOverlap>0,density,fill=`EVE.ClusterOverlap`>0))+
  geom_boxplot(aes(EVE.ClusterOverlap>0,density),width=0.3)+
  scale_y_log10()

Wp<-wilcox.test(piClusterInfo$density~piClusterInfo$EVE.ClusterOverlap>0)
#View(EVEpiInfo)
virus=as.factor(strsplit2(strsplit2(EVEpiInfo$ID,"\\[")[,2],"\\]")[,1])
piPropperEVE<-ggplot(data=EVEpiInfo)+geom_point(aes(piEVEProp,EVEpiCover,color=family))+scale_x_log10()+scale_y_log10()
piPosReads<-ggplot(data=EVEpiInfo)+geom_point(aes(piReadsPerPos,piEVEProp,size=piEVEProp,color=family))+scale_x_log10()+scale_y_log10()

pdf("EVE_pi_barcharts.pdf",width=3,height=4)
EVEpiInfo$family<-with(EVEpiInfo,reorder(family,family,length))

B=ggplot(data=EVEpiInfo,lwd=1)+theme_minimal()

plot(chart_EVECount<-B+geom_bar(aes(x=family,fill=family),color='black',lwd=0.2,width=0.5,stat='count',position='dodge')+
  scale_fill_brewer(palette = "Spectral")+xlab("EVE Virus Family")+ylab("Number of EVEs")+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "None")
)

plot(chart_EVEbp<-B+geom_bar(aes(x=family,y=EVEPos,fill=family),color='black',lwd=0.2,width=0.6,stat='identity',position='stack')+
  scale_fill_brewer(palette = "Spectral")+xlab("EVE Virus Family")+ylab("Length of EVEs")+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "None")
)

plot(chart_piRNAs<-B+geom_bar(aes(x=family,y=piEVEProp,fill=family),color='black',lwd=0.2,width=0.6,stat='identity',position='stack')+
  scale_fill_brewer(palette = "Spectral")+xlab("EVE Virus Family")+ylab("Proportion of piRNA Reads")+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "None")
)


plot(chart_piRNAs<-B+geom_bar(aes(x=reorder(species,piEVEProp,sum),y=piEVEProp,fill=family),color='black',lwd=0.2,width=0.6,stat='identity',position='stack')+
       scale_fill_brewer(palette = "Spectral")+xlab("EVE Virus Family")+ylab("Proportion of piRNA Reads")+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "None")
)

plot(chart_piRNAs<-B+geom_bar(aes(x=reorder(species,propGenome,sum),y=propGenome,fill=family),color='black',lwd=0.2,width=0.6,stat='identity',position='stack')+
       scale_fill_brewer(palette = "Spectral")+xlab("EVE Virus Species")+ylab("Proportion of genome")+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "None")
)

plot(chart_piRNAs<-B+geom_bar(aes(x=reorder(species,piEVEProp/propGenome,sum),y=piEVEProp/propGenome,fill=family),color='black',lwd=0.2,width=0.6,stat='identity',position='stack')+
       scale_fill_brewer(palette = "Spectral")+xlab("EVE Virus Species")+ylab("Normalized")+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "None")
)

plot(chart_piRNAs<-B+geom_bar(aes(x=reorder(species,piEVEProp/propGenome,sum),y=piEVEProp/propGenome,fill=family),color='black',lwd=0.2,width=0.6,stat='identity',position='stack')+
       scale_fill_brewer(palette = "Spectral")+xlab("EVE Virus Species")+ylab("Normalized")+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "None")
)

dev.off()

B+scale_fill_brewer(palette = "Spectral")+
  geom_boxplot(aes(family,piEVEProp,fill=family),color="black",scale = 'width')+
  geom_jitter(width = .2,pch=21,aes(family,piEVEProp),fill="black")+
  #geom_text(data =EVEpiInfo[EVEpiInfo$enrichGenome>1500,],nudge_x = -0.02,hjust = 1,aes(EVEPos,piEVEProp,label=species),color="black",alpha=0.7)+
  scale_y_log10()+ylab("Proportion of piRNA Reads")+
  xlab("EVE Length")#+scale_size(trans='sqrt',breaks = c(0,10,100,1000))

bs<-as.data.frame(with(EVEpiInfo,beeswarm(log10(piEVEProp+.00000001))))
bs$att<-EVEpiInfo$family

ggplot(bs)+geom_point(aes(x,y,color=att))+scale_color_brewer(palette = "Spectral")

B+geom_point(aes(EVEPos,piEVEProp,color=piReadsPerPos,size=piReadsPerPos))+scale_x_log10()

chart_meanLength<-
       ggplot(data=EVEpiInfo,aes(x=reorder(family,family,length)))+
       geom_bar(width=0.6,stat='identity',position='dodge')+scale_fill_brewer(palette = "Spectral")+xlab("EVE Virus Family")
#dev.off()

ggplot(piClusterInfo)+scale_fill_viridis(discrete = T)+
  geom_histogram(aes(piClust_reads/piClust_piSites,fill=EVEhits>0),color="black",alpha=0.7)+scale_x_log10()

ggplot(piClusterInfo)+geom_point(aes(EVEhits,piClust_reads))+scale_y_log10()+geom_smooth(aes(EVEhits,piClust_reads),formula = piClust_reads~EVEhits)

plot(table(EVEpiInfo$piPos))
table(EVEpiInfo$EVEpiCover>0)

#  scale_x_log10()+ylab("")
#  xlab("EVE Length")+scale_size(trans='sqrt',breaks = c(10,100,1000))

piPosTotal<-sum(EVEpiInfo$piPos)
#plot(EVEpiInfo)
piReadsTotal<-sum(EVEpiInfo$piReads)
piEVEReadsTotal<-sum(EVEpiInfo$piEVEreads)
piEVEReadsTotal/piReadsTotal

piEVEposTotal<-sum(EVEpiInfo$piEVEPos)
evePosTotal<-sum(EVEpiInfo$EVEPos)
piEVEposTotal/evePosTotal

sum(TEinfo$teSites)/sum(TEinfo$length)

EVEplot=ggplot(data = TEinfo )+theme_bw()+geom_histogram(aes(TEinfo$piReads.EVEProp),binwidth=.025)+xlab("TE Content")+ylab("Number of Contigs")

sum(TEinfo$piReads)
sum(TEinfo$piEVEReads)/sum(TEinfo$piReads)

interval<-strsplit2(EVEpiInfo$ID,"F_")[,2]
EVEpiInfo$start<-strsplit2(interval,"_")[,1]
EVEpiInfo$end<-strsplit2(interval,"_")[,2]
EVEpiInfo$EVEdir<-ifelse(EVEpiInfo$end>EVEpiInfo$start,"+","-")
EVEpiInfo$sense<-ifelse(EVEpiInfo$EVEdir=="-",EVEpiInfo$negreads,EVEpiInfo$posreads)
EVEpiInfo$antisense<-ifelse(EVEpiInfo$EVEdir=="-",EVEpiInfo$posreads,EVEpiInfo$negreads)


#compare EVEs by family
Reads<-dcast(EVEpiInfo,formula = family~.,value.var = "piEVEreads",fun.aggregate = sum)
senseReads<-dcast(EVEpiInfo,formula = family~.,value.var = "sense",fun.aggregate = sum)
antisenseReads<-dcast(EVEpiInfo,formula = family~.,value.var = "antisense",fun.aggregate = sum)
PropGenome<-dcast(EVEpiInfo,formula = family~.,value.var = "propGenome",fun.aggregate = sum)
Count<-dcast(EVEpiInfo,formula = family~.,value.var = "propGenome")

familyTable <-data.frame(Count,propGenome=PropGenome[,2],reads=Reads[,2],sense=senseReads[,2],antisense=antisenseReads[,2])

ggplot(familyTable,aes(propGenome,reads,fill=family))+
  geom_point(aes(cex = .),pch=21)+
  ylab("Reads")+
  xlab("Proportion of Genome")+
  scale_fill_brewer(palette = "Spectral")+
  scale_x_log10()+
  scale_y_log10()+
  scale_size(range = c(2,9))+
  geom_text(hjust=0,nudge_x = .05,aes(label=family))

ggplot(familyTable,aes(propGenome,sense,fill=family))+
  geom_point(aes(cex = .),pch=21)+
  ylab("Reads")+
  xlab("Proportion of Genome")+
  scale_fill_brewer(palette = "Spectral")+
  scale_x_log10()+
  scale_y_log10()+
  scale_size(range = c(2,9))+
  geom_text(hjust=0,nudge_x = .05,aes(label=family))

ggplot(familyTable,aes(propGenome,antisense,fill=family))+
  geom_point(aes(cex = .),pch=21)+
  ylab("Reads")+
  xlab("Proportion of Genome")+
  scale_fill_brewer(palette = "Spectral")+
  scale_x_log10()+
  scale_y_log10()+
  scale_size(range = c(2,9))+
  geom_text(hjust=0,nudge_x = .05,aes(label=family))

thing<-strsplit2(mergeCluster$`Binding sites`,")")
DFlist<-apply(thing,1,function(X){data.frame(strsplit2(X,split = " \\("))})
n=0
all<-data.frame()
for(i in DFlist){
  df<-i
  n=n+1
  sel<-df[df$X1!="",]
  if(nrow(sel)>0){
    sel$ClusterID<-mergeCluster$ClusterID[n]
    sel$TF<- stringr::str_trim(sel$X1)
    sel$seq<-strsplit2(as.character(sel$X2),split = "[0-9]")[,1]
    sel$position<- as.integer(stringr::str_replace_all(sel$X2,"[ATCG]",""))
    all<-rbind(all,sel)
  }
}
head(all)

#ggplot(all)+geom_point(aes(position,ClusterID))

superMerge<-merge(all,mergeCluster)
superMerge<-data.table(superMerge,stringsAsFactors = T)
superMerge[,adjStart:=(position + Start + allContigs$startpositions[as.character(allContigs$CONTIG)==as.character(Contig)[1]]),ClusterID]


ggplot(superMerge)+geom_point(aes(bias,adjStart/Length))+facet_wrap(~TF)
ggplot(superMerge)+geom_boxplot(aes(TF,count))+scale_y_log10()+facet_grid(~TF)

ggplot(superMerge)+geom_point()

superMerge$seq<-factor(superMerge$seq)
seqCount<-dcast(superMerge,formula=ClusterID+TF+EVEhits+piClust_reads+`Hits (normalized) per kb`+`Size [bp]`+`Predicted directionality`+seq~.)
seqCount$count<-seqCount$./8
seqCount$length<-seqCount$`Size [bp]`
seqCount$dir<-limma::strsplit2(seqCount$`Predicted directionality`," \\(")[,1]
ggsave(ggplot(seqCount)+geom_jitter(aes(cex=`Size [bp]`,piClust_reads/`Size [bp]`,count,color=dir),alpha=0.6)+facet_wrap(~TF+seq,drop = T)+scale_y_log10()+scale_x_log10(),filename = "plot.pdf",width = 20,height = 20)

ggplot(seqCount)+geom_jitter(aes(x=dir,y=seq ,cex=length,color=dir))

dcast(seqCount,TF+seq~dir,value.var = "count",sum)




SeqBypiDir<-dcast(seqCount,seq+TF~dir,value.var = "count",fun.aggregate = sum)

SeqBypiDir<-dcast(seqCount,seq+TF~dir,value.var = "count",fun.aggregate = sum)
SeqBypiDirL<-dcast(seqCount,seq+TF~dir,value.var = "length",fun.aggregate = sum)

SeqBypiDirEVE<-dcast(seqCount,seq+TF~dir+EVEhits>0,value.var = "count",fun.aggregate = sum)
SeqBypiDirEVEL<-dcast(seqCount,seq+TF~dir+EVEhits>0,value.var = "length",fun.aggregate = sum)

TFBypiDir<-dcast(seqCount,TF~dir,value.var = "count",fun.aggregate = sum)
TFBypiDirL<-dcast(seqCount,TF~dir,value.var = "length",fun.aggregate = sum)

TFByDirByEVE<-dcast(seqCount,TF~dir+EVEhits>0,value.var = "count",fun.aggregate = sum)
TFByDirByEVEL<-dcast(seqCount,TF~dir+EVEhits>0,value.var = "length",fun.aggregate = sum)

TFByEVE<-dcast(seqCount,TF~EVEhits>0,value.var = "count",fun.aggregate = sum)
TFByEVEL<-dcast(seqCount,TF~EVEhits>0,value.var = "length",fun.aggregate = sum)

SeqByEVEs<-dcast(seqCount,seq+TF~EVEhits>0,value.var = "count",fun.aggregate = length)
SeqByEVEsL<-dcast(seqCount,seq+TF~EVEhits>0,value.var = "length",fun.aggregate = sum)
#dirandeves
#

Chis<-function(Y,Z){
  pZ<-apply(Z,1,function(X){X/sum(X)})
  Y<-cbind(Y,t(pZ))
  output<-apply(Y,1,CHI)
}
chisq.test()

CHI<-function(Y){
  br<-length(Y)/2
  E<-sum(Y[1:br])*Y[(br+1):length(Y)]
  print(paste("E:",E))
  print(paste("O:",Y[1:br]))
  ChiStat=0
  for (i in 1:br){
    ChiStat<-sum((((Y[i]-E[i])^2)/E[i]),ChiStat,na.rm = T)
    print(ChiStat)
  }
  p<-pchisq(log.p = T,q = ChiStat,df = (br-1) ,lower.tail=FALSE)
  return(p)}


SeqBypiDir$bias<-ifelse(SeqBypiDir$`mono:minus`>SeqBypiDir$`mono:plus`,yes = 1,no=-1)
SeqBypiDir$pvals<-Chis(SeqBypiDir[,4:5],SeqBypiDirL[,4:5])


SeqByEVEs$bias<-ifelse(SeqByEVEs$`TRUE`>SeqByEVEs$`FALSE`,yes = 1,no=-1)
SeqByEVEs$pvals<-Chis(SeqByEVEs[,3:4],SeqByEVEsL[,3:4])

TFByDirByEVE$pvals<-Chis(TFByDirByEVE[,2:6],TFByDirByEVEL[,2:6],df = 5)

TFByEVE$bias<-ifelse(TFByEVE$`TRUE`>TFByEVE$`FALSE`,yes = 1,no=-1)
TFByEVE$pvals<-Chis(TFByEVE[,2:3],TFByEVEL[,2:3])

SeqByEVEs$bias<-ifelse(SeqByEVEs$`TRUE`>SeqByEVEs$`FALSE`,yes = 1,no=-1)
SeqByEVEs$pvals<-Chis(SeqByEVEs[,3:4],SeqByEVEsL[,3:4])

SeqBypiDirEVE$pvals<-Chis(SeqBypiDirEVE[,3:7],SeqBypiDirEVEL[,3:7])

levels(SeqBypiDir)

ggplot(SeqBypiDirEVE)+geom_bar(color="black",aes(x = reorder(seq,TF,c), y=pvals,fill=TF),stat = 'identity')+theme(axis.text.x = element_text(angle=90))
ggplot(SeqBypiDir)+geom_bar(color="black",aes(x = reorder(seq,bias*pvals,mean), y=bias*pvals,fill=TF),stat = 'identity')+theme(axis.text.x = element_text(angle=90))
SeqBypiDir$padj<-log10(p.adjust(10^(SeqBypiDir$pvals),method = 'hommel'))
ggplot(SeqBypiDir)+geom_hline(yintercept=c(log10(0.05),-log10(0.05)))+geom_bar(color="black",aes(x = reorder(seq,bias*padj,mean), y=bias*padj,fill=TF,alpha=padj<log10(0.05)),stat = 'identity')+theme(axis.text.x = element_text(angle=90))

ggplot(SeqBypiDir)+geom_abline(slope=1)+geom_point(aes(`mono:minus`,`mono:plus`,color=bias, cex=(-padj),fill=TF,alpha=padj<log10(0.05)),stat = 'identity')+theme(axis.text.x = element_text(angle=90))+scale_x_log10()+scale_y_log10()

ggplot(SeqBypiDir)+geom_bar(aes(x = reorder(seq,pvals,mean), y=pvals,fill=TF),stat = 'identity')+theme(axis.text.x = element_text(angle=90))
ggplot(SeqByEVEs)+geom_bar(aes(x = reorder(seq,TF,c), y=bias*pvals,fill=TF),stat = 'identity')+theme(axis.text.x = element_text(angle=90))
ggplot(TFByDirByEVE)+geom_bar(aes(x = TF, y=pvals,fill=TF),stat = 'identity')+theme(axis.text.x = element_text(angle=90))
ggplot(TFByEVE)+geom_bar(aes(x = TF, y=bias*pvals,fill=TF),stat = 'identity')+theme(axis.text.x = element_text(angle=90))
