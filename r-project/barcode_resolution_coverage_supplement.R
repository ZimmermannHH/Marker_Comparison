############################################
#12.06.2023
#barcode resolution and taxonomic coverage
#author: Heike Zimmermann
############################################
#R v 4.1.3
library(stringr) #v. 1.5.0
library(tibble)  #v. 3.2.1
library(tidyr)   #v. 1.3.0
library(dplyr)   #v. 1.1.2
library(ggplot2) #v. 3.4.2 
library(phylotools) #to read-in fasta files v. 0.2.2
library(rlang)   # v. 1.1.0
library(patchwork) #v. 1.1.2

#read-in exotaxspecificity output
resolution.SAR<-read.table("ecotaxspec_SAR.txt", sep="\t", header=TRUE,quote="", comment.char="#")
resolution.fg<-read.table("ecotaxspec_focusgroups.txt", sep="\t", header=TRUE,quote="", comment.char="#")
str(resolution.SAR)
#subset for specific ranks
resolution.SAR.sub<-subset(resolution.SAR, rank=="family" |rank=="genus" |rank=="species")
resolution.fg.sub<-subset(resolution.fg, rank=="family" |rank=="genus" |rank=="species")

#from character to factor
resolution.SAR.sub$rank=as.factor(resolution.SAR.sub$rank)
resolution.SAR.sub$marker=as.factor(resolution.SAR.sub$marker)
resolution.SAR.sub$group=as.factor(resolution.SAR.sub$group)

resolution.fg.sub$rank=as.factor(resolution.fg.sub$rank)
resolution.fg.sub$marker=as.factor(resolution.fg.sub$marker)
resolution.fg.sub$group=as.factor(resolution.fg.sub$group)

#Taxonomic resolution for SAR group
pdf("resolution_SAR.pdf", height = 5, width = 6)
ggplot(resolution.SAR.sub, aes(x=marker, y=percent, fill=factor(group))) + 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(limits = rev) +
  theme(legend.position = "none",axis.text=element_text(size=8),
        axis.title=element_text(size=8))+
  scale_fill_manual(values=c( "#264653","#8d96a3","#2A9D8F"))+
  geom_bar(stat="identity", position=position_dodge(width=0.9),)+
  labs(y = "Percent (%)", x ="")+
  geom_text(aes(y=percent, label=round(percent, 0)), position=position_dodge(width=0.9), color="white", hjust=1.15, size=3.25) + 
  theme(legend.title=element_blank(), legend.position="bottom", 
        axis.title.x=element_text(face="bold"), axis.title.y=element_blank(), 
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
        panel.background=element_blank(), 
        panel.grid.minor=element_blank(), 
        #panel.grid.major=element_blank(),   
        strip.text.y = element_text(angle=0)) + 
  facet_wrap(~ rank,scales = "fixed", nrow = 1)
dev.off()

#taxonomic resolution for focus groups
pdf("resolution_focus_groups.pdf", height = 5, width = 6)
ggplot(resolution.fg.sub, aes(x=marker, y=percent, fill=factor(group))) + 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(limits = rev) +
  theme(legend.position = "none",axis.text=element_text(size=8),
        axis.title=element_text(size=8))+
  scale_fill_manual(values=c( "#E9C46A" ,"#264653","#8d96a3","#2A9D8F"))+
  geom_bar(stat="identity", position=position_dodge(width=0.9),)+
  labs(y = "Percent (%)", x ="")+
  geom_text(aes(y=percent, label=round(percent, 0)), position=position_dodge(width=0.9), color="white", hjust=1.15, size=3.25) + 
  theme(legend.title=element_blank(), legend.position="bottom", 
        axis.title.x=element_text(face="bold"), axis.title.y=element_blank(), 
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
        panel.background=element_blank(), 
        panel.grid.minor=element_blank(), 
        #panel.grid.major=element_blank(),   
        strip.text.y = element_text(angle=0)) + 
  facet_wrap(~ rank,scales = "fixed", nrow = 1)
dev.off()