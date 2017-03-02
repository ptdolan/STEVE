#TE EVE
library(limma)
library(ggplot2)
#byContig
TEinfo<-read.csv("~/GitHub/stEVE/EVE_TEpi_Analysis_byContig.csv")
#byEVE
EVEpiInfo<-read.csv("~/GitHub/stEVE/allStatsbyEVE.csv")
#bypiCluster
piClusterInfo<-read.csv("~/github/steve/allStats_byCluster.csv")

ggplot(piClusterInfo)+
  geom_bar(aes(with(piClusterInfo,reorder(X,piClust_reads,mean)),piClust_reads),stat="identity")+
  geom_bar(aes(with(piClusterInfo,reorder(X,piClust_reads,mean)),ClusterEVEpiReads),stat="identity",color='red')+scale_y_log10()

#View(EVEpiInfo)

TEcontent=ggplot(data = TEinfo)+geom_density(aes(teProp),fill="grey")+xlab("TE Content")+ylab("Density of Contigs")+geom_vline(aes(xintercept=sum(TEinfo$teSites)/sum(TEinfo$length)),color="darkred")
TEcontent=ggplot(data = TEinfo)+theme_bw()+geom_density(aes(teProp),fill="grey")+xlab("TE Content")+ylab("Number of Contigs")+geom_vline(aes(xintercept=sum(TEinfo$teSites)/sum(TEinfo$length)),color="darkred")
TEcontent=ggplot(data = TEinfo)+theme_bw()+geom_histogram(aes(teSites),binwidth=.025)+xlab("TE Content")+ylab("Number of Contigs")+geom_vline(aes(xintercept=sum(TEinfo$teSites)/sum(TEinfo$length)),color="darkred")

virus=as.factor(strsplit2(strsplit2(EVEpiInfo$ID,"\\[")[,2],"\\]")[,1])
piPropperEVE<-ggplot(data=EVEpiInfo)+geom_point(aes(piEVEProp,EVEpiCover,color=family))+scale_x_log10()+scale_y_log10()
piPosReads<-ggplot(data=EVEpiInfo)+geom_point(aes(piReadsPerPos,piEVEProp,size=piEVEProp,color=family))+scale_x_log10()+scale_y_log10()
ggplot(EVEpiInfo)+coord_polar()+geom_area(aes(x=ID, y=piEVEProp,fill=family),position='dodge',color="black",stat = 'identity')+theme(legend.position = "None")
#ggsave(TEcontent,filename = "TEcontent_byBase.pdf",width=5,height=3)

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

dev.off()

B+scale_fill_brewer(palette = "Spectral")+geom_point(pch=21,aes(EVEPos,piEVEProp,fill=family,cex=enrichGenome),color="black",alpha=0.7)+scale_x_log10()+ylab("Proportion of piRNA Reads")+xlab("EVE Length")+scale_size(trans='sqrt',breaks = c(0,10,100,1000))

B+geom_point(aes(EVEPos,piEVEProp,cex=piReadsPerPos),color='grey')+scale_x_log10()+
  geom_point(data=EVEpiInfo[EVEpiInfo$piEVEProp>.0005,],aes(EVEPos,piEVEProp,cex=piReadsPerPos,color=contig),fill="black")
  
B+geom_point(aes(EVEPos,piEVEProp,color=piReadsPerPos,size=piReadsPerPos))+scale_x_log10()

chart_meanLength<-
       ggplot(data=EVEpiInfo,aes(x=reorder(family,family,length)))+
       geom_bar(aes(y=),width=0.6,stat='identity',position='dodge')+scale_fill_brewer(palette = "Spectral")+xlab("EVE Virus Family")
dev.off()


piPosTotal<-sum(EVEpiInfo$piPos)
piReadsTotal<-sum(EVEpiInfo$piReads)
piEVEReadsTotal<-sum(EVEpiInfo$piEVEreads)
piEVEReadsTotal<-sum(EVEpiInfo$piEVEPos)
eveReadsTotal<-sum(EVEpiInfo$EVEPos)

sum(TEinfo$teSites)/sum(TEinfo$length)

EVEplot=ggplot(data = TEinfo )+theme_bw()+geom_histogram(aes(TEinfo$piReads.EVEProp),binwidth=.025)+xlab("TE Content")+ylab("Number of Contigs")

sum(TEinfo$piReads)
sum(TEinfo$piEVEReads)/sum(TEinfo$piReads)
