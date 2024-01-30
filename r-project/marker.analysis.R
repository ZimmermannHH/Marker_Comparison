###################################
#12.06.2023
#marker analysis
#author: Heike Zimmermann
###################################
#R v 4.1.3
library(stringr) #v. 1.5.0
library(tibble)  #v. 3.2.1
library(tidyr)   #v. 1.3.0
library(dplyr)   #v. 1.1.2
library(ggplot2) #v. 3.4.2 
library(phylotools) #to read-in fasta files v. 0.2.2
library(rlang)   # v. 1.1.0
library(patchwork) #v. 1.1.2
#devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram") #v. 1.2.2


setwd("lineages")
#
#read-in lineage information
V9_Amaral1<-read.table("lineages.Amaral1_V9_e3.txt", sep="\t", header=FALSE,quote="", comment.char="#")
V9_1389F_1510R<-read.table("lineages.V9_1389F_1510R_e3.txt", sep="\t", header=FALSE,quote="", comment.char="#")
V9_Piredda<-read.table("lineages.Piredda_V9_e3.txt", sep="\t", header=FALSE,quote="", comment.char="#")
V9_Stoeck<-read.table("lineages.Stoeck_V9_e3.txt", sep="\t", header=FALSE,quote="", comment.char="#")
V7_1183mod_1443mod<-read.table("lineages.V7_1183mod_1443mod_e3.txt", sep="\t", header=FALSE,quote="", comment.char="#")
V7_960F_NSR1438<-read.table("lineages.V7_960F_NSR1438_e3.txt", sep="\t", header=FALSE,quote="", comment.char="#")
V7_18S_allshorts<-read.table("lineages.V7_18S_allshorts_correctprimer_e3.txt", sep="\t", header=FALSE,quote="", comment.char="#")
V7_18S_allshorts_modF<-read.table("lineages.V7_18S_allshorts_modF_e3.txt", sep="\t", header=FALSE,quote="", comment.char="#")
V3_SAR_V3<-read.table("lineages.SAR_V3_e3.txt", sep="\t", header=FALSE,quote="", comment.char="#")
V4<-read.table("lineages.V4_E572F_E1009R_e3.txt", sep="\t", header=FALSE,quote="", comment.char="#")

#read-in fasta files
V9_Amaral1.fa<-read.fasta("Amaral1_V9_e3_embl143_final_db.fasta")
V9_1389F_1510R.fa<-read.fasta("V9_1389F_1510R_e3_embl143_final_db.fasta")
V9_Piredda.fa<-read.fasta("Piredda_V9_e3_embl143_final_db.fasta")
V9_Stoeck.fa<-read.fasta("Stoeck_V9_e3_embl143_final_db.fasta")
V7_1183mod_1443mod.fa<-read.fasta("V7_1183mod_1443mod_e3_embl143_final_db.fasta")
V7_960F_NSR1438.fa<-read.fasta("V7_960F_NSR1438_e3_embl143_final_db.fasta")
V7_18S_allshorts.fa<-read.fasta("V7_18S_allshorts_correctprimer_e3_embl143_final_db.fasta")
V7_18S_allshorts_modF.fa<-read.fasta("V7_18S_allshorts_modF_e3_embl143_final_db.fasta")
V3_SAR_V3.fa<-read.fasta("SAR_V3_e3_embl143_final_db.fasta")
V4.fa<-read.fasta("V4_E572F_E1009R_e3_embl143_final_db.fasta")

setwd("..")

#for each marker split lineage info into separate columns, add marker name 
# remove unnecessary characters and combine datasets

list.markers<-c("V9_Amaral1", "V9_1389F_1510R", "V9_Piredda", "V9_Stoeck", 
                "V7_1183mod_1443mod", "V7_960F_NSR1438", "V7_18S_allshorts", "V7_18S_allshorts_modF",
                "V3_SAR_V3", "V4")
updated_lst_markers <- list()

for (i in list.markers) {
  tmp <- get(i)
  tmp.sep<-separate(data=tmp,col = V2, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|")
  tmp.sep$marker=paste(i)
  tmp.fasta=paste0(i, ".fa") #add the sequence from the corresponding fasta (name of i with .fa extension)
  tmp.fasta2=get(tmp.fasta)
  tmp.sep$sequence=tmp.fasta2$seq.text
  tmp.sep$Kingdom=gsub("k__", "",tmp.sep$Kingdom)
  tmp.sep$Phylum=gsub("p__", "",tmp.sep$Phylum)
  tmp.sep$Class=gsub("c__", "",tmp.sep$Class)
  tmp.sep$Order=gsub("o__", "",tmp.sep$Order)
  tmp.sep$Family=gsub("f__", "",tmp.sep$Family)
  tmp.sep$Genus=gsub("g__", "",tmp.sep$Genus)
  tmp.sep$Species=gsub("s__", "",tmp.sep$Species)
  updated_lst_markers[[i]]<- tmp.sep
}


marker.lineages.df<-data.table::rbindlist(updated_lst_markers)
head(marker.lineages.df)
names(marker.lineages.df)[1]="NCBI_taxid"
write.table(marker.lineages.df, "marker.lineages.df.txt", sep="\t")

#########################################
marker.lineages.df<-read.table("marker.lineages.df.txt", sep="\t", header=TRUE)
marker.lineages.df$NCBI_taxid=as.factor(marker.lineages.df$NCBI_taxid)
marker.lineages.df=subset(marker.lineages.df, marker!="V7_18S_allshorts_modF")
marker.lineages.df$marker=as.factor(marker.lineages.df$marker)
#replace primer names
marker.lineages.df$marker=gsub("V9_Amaral1", "V9_1380F-1510R", marker.lineages.df$marker)
marker.lineages.df$marker=gsub("V9_Piredda", "V9_1388F-1510R", marker.lineages.df$marker)
marker.lineages.df$marker=gsub("V9_Stoeck", "V9_1391F- EukB", marker.lineages.df$marker)
marker.lineages.df$marker=gsub("V4", "V4_E572F-E1009R", marker.lineages.df$marker)


#unique entries per rank
length(unique(marker.lineages.df$Phylum))
#61
length(unique(marker.lineages.df$Class))
#271
length(unique(marker.lineages.df$Family))
#5832
length(unique(marker.lineages.df$Genus))
#23464
length(unique(marker.lineages.df$Species))
#49639


#1) Phyla etc. per marker

no.phyla<-marker.lineages.df %>%
  group_by(marker) %>%
  summarise(n.phyla = length(unique(Phylum)))

no.class<-marker.lineages.df %>%
  group_by(marker) %>%
  summarise(n.class = length(unique(Class)))

no.family<-marker.lineages.df %>%
  group_by(marker) %>%
  summarise(n.family = length(unique(Family)))

no.genus<-marker.lineages.df %>%
  group_by(marker) %>%
  summarise(n.genus = length(unique(Genus)))

no.species<-marker.lineages.df %>%
  group_by(marker) %>%
  summarise(n.species = length(unique(Species)))

taxlevels2<-cbind(no.phyla, no.class$n.class, no.family$n.family, no.genus$n.genus, no.species$n.species)
names(taxlevels2)[2:6]=c("01.Phylum","02.Class", "03.Family", "04.Genus", "05.Species")
taxlevels.long2 <- gather(taxlevels2, taxlevel,frequency,  "01.Phylum":"05.Species", factor_key=TRUE)

taxlevels.long2$marker=as.factor(taxlevels.long2$marker)

ggplot(data=taxlevels.long2, aes(x=marker, y=frequency, fill=taxlevel)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+ 
  coord_flip()+ 
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#E76F51","#E9C46A", "#264653", "#8d96a3","#66a182"))+
  facet_wrap(~ taxlevel, scales = "free_x", nrow = 1)


#2) Amplicon length distribution per marker

marker.lineages.df$sequence.length=nchar(marker.lineages.df$sequence)

summary.size<-marker.lineages.df %>%
  group_by(marker) %>%
  summarise(length.mean = mean(nchar(sequence)), n = n(),
            length.median = median(nchar(sequence)),
            length.min = min(nchar(sequence)),
            length.max = max(nchar(sequence)))

seqlenmarker <-marker.lineages.df %>% 
  group_by(marker) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)
names(seqlenmarker)

pdf("seqlength.marker_violin.pdf", width = 6, height = 6, bg="transparent")
ggplot(seqlenmarker, aes(x=marker, y=sequence.length, fill=marker)) + 
  theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_violin(trim=TRUE)
dev.off()


#2) Amplicon length distribution per phyla per marker
#maybe grouping protists according to Burki et al. 2019: https://www.cell.com/trends/ecology-evolution/fulltext/S0169-5347%2819%2930257-5

no.phyla<-marker.lineages.df %>%
  group_by(marker) %>%
  summarise(n.phyla = length(unique(Phylum)))

seqlenmarker <-marker.lineages.df %>% 
  group_by(marker, Phylum) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)
names(seqlenmarker)

#Violin plot
ggplot(seqlenmarker, aes(x=Phylum, y=sequence.length, fill=marker)) + 
  theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none",axis.text = element_text(size = 8))+
  scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_violin(trim=TRUE)+
  facet_wrap(~ marker,scales = "free_x", nrow = 2)

#scatter plot amplicon size rank: phyla
#pdf("seqlength_phylum_scatter.pdf", height = 40, width = 40)
ggplot(seqlenmarker, aes(x=Phylum, y=sequence.length, fill=marker)) + 
  #theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none", axis.text = element_text(size = 14))+
  geom_point(size=.25, alpha=.25)+
  facet_wrap(~ marker,scales = "fixed", nrow = 3)
#dev.off()

#violin plot amplicon size rank: Kingdom
#pdf("ampliconsize_kingdom.pdf", height = 5, width = 10)
ggplot(seqlenmarker, aes(x=Kingdom, y=sequence.length, fill=marker)) + 
  #theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none",strip.text = element_text(size=7))+
  scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_violin(trim=TRUE)+
  labs(y="Amplicon length (bp)")+
  facet_wrap(~ marker,scales = "fixed", nrow = 2)
#dev.off()


################################################

#co-amplification 

#Bacteria

bac.sub<-subset(marker.lineages.df, Kingdom=="Bacteria")
dim(bac.sub)

bacteria.seq.per.phylum <-bac.sub %>% 
  group_by(marker) %>%
  distinct(NCBI_taxid, .keep_all = TRUE) %>% 
  count(Phylum)

ggplot(bacteria.seq.per.phylum, aes(x=Phylum, y=n, fill=marker)) + 
  coord_flip()+  
  scale_y_log10()+
  scale_x_discrete(limits = rev)+
  theme(legend.position = "none")+
  geom_point()+
  ylab("number of taxa") +
  facet_wrap(~ marker,scales = "fixed", nrow = 2)


#Archaea

ar.sub<-subset(marker.lineages.df, Kingdom=="Archaea")
archaea.seq.per.phylum <-ar.sub %>% 
  group_by(marker) %>%
  distinct(NCBI_taxid, .keep_all = TRUE) %>% 
  count(Phylum)

ggplot(archaea.seq.per.phylum, aes(x=Phylum, y=n, fill=marker)) + 
  coord_flip()+  
  scale_y_log10()+
  theme(legend.position = "none")+
  scale_x_discrete(limits = rev)+
  geom_point()+
  ylab("number of taxa") +
  facet_wrap(~ marker,scales = "fixed", nrow = 1)


#Viruses

vir.sub<-subset(marker.lineages.df, Kingdom=="Viruses")
dim(vir.sub)
viruses.seq.per.phylum <-vir.sub %>% 
  group_by(marker) %>%
  distinct(NCBI_taxid, .keep_all = TRUE) %>% 
  count(Phylum)

ggplot(viruses.seq.per.phylum, aes(x=Phylum, y=n, fill=marker)) + 
  coord_flip()+  
  scale_y_log10()+
  scale_x_discrete(limits = rev)+
  theme(legend.position = "none")+
  geom_point()+
  ylab("number of taxa") +
  facet_wrap(~ marker,scales = "fixed", nrow = 2)


#Fungi
fungi<-c("Ascomycota", "Basidiomycota","Mucoromycota","Zoopagomycota","Blastocladiomycota", "Olpidiomycota", "Chytridiomycota")
fu.sub<-subset(marker.lineages.df, Phylum %in% fungi)
fungi.seq.per.phylum <-fu.sub %>% 
  group_by(marker) %>%
  distinct(NCBI_taxid, .keep_all = TRUE) %>% 
  count(Phylum)

ggplot(fungi.seq.per.phylum, aes(x=Phylum, y=n, fill=marker)) + 
  coord_flip()+  
  scale_y_log10()+
  scale_x_discrete(limits = rev)+
  theme(legend.position = "none",strip.text = element_text(size=7))+
  geom_point()+
  ylab("number of taxa") +
  facet_wrap(~ marker,scales = "fixed", nrow = 2)


#Metazoa
metazoa<-c("Arthropoda", "Platyhelminthes",
           "Mollusca", "Nematoda", "Annelida",
           "Porifera", "Cnidaria",                   
           "Tardigrada", "Acanthocephala","Ctenophora",                 
           "Chordata", "Nematomorpha", 
           "Bryozoa",                    
           "Echinodermata", "Nemertea","Rotifera",                   
           "Brachiopoda" , "Picozoa" ,"Onychophora" ,               
           "Priapulida","Chaetognatha",
           "Hemichordata")

me.sub<-subset(marker.lineages.df, Phylum %in% metazoa)

metazoa.seq.per.phylum <-me.sub %>% 
  group_by(marker) %>%
  distinct(NCBI_taxid, .keep_all = TRUE) %>% 
  count(Phylum)

ggplot(metazoa.seq.per.phylum, aes(x=Phylum, y=n, fill=marker)) + 
  coord_flip()+  
  scale_y_log10()+
  theme(legend.position = "none")+
  scale_x_discrete(limits = rev)+
  geom_point()+
  ylab("number of taxa") +
  facet_wrap(~ marker,scales = "fixed", nrow = 2)


##############################################################
#taxonomic resolution
res<-read.table("ecotaxspecificity_markers.txt", sep="\t", header=TRUE)


#res.sub<-subset(res, rank=="phylum" | rank=="class" |rank=="order" |rank=="family" |rank=="genus" |rank=="species")
res.sub<-subset(res, rank=="family" |rank=="genus" |rank=="species")
str(res.sub)
res.sub$rank=as.factor(res.sub$rank)
res.sub$marker=as.factor(res.sub$marker)
#replace primer names
res.sub$marker=gsub("V9_Amaral1", "V9_1380F-1510R", res.sub$marker)
res.sub$marker=gsub("V9_Piredda", "V9_1388F-1510R", res.sub$marker)
res.sub$marker=gsub("V9_Stoeck", "V9_1391F- EukB", res.sub$marker)
res.sub$marker=gsub("V4", "V4_E572F-E1009R", res.sub$marker)


res1<-ggplot(res.sub, aes(x=marker, y=percent)) + 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(limits = rev) +
  theme(legend.position = "none",axis.text=element_text(size=8),
        axis.title=element_text(size=8))+
  geom_bar(stat="identity")+
  labs(y = "Percent (%)", x ="")+
  facet_wrap(~ rank,scales = "fixed", nrow = 1)

res2<-ggplot(res.sub, aes(x=marker, y=taxon_ok)) + 
  theme_minimal()+
  coord_flip()+
  scale_x_discrete(limits = rev) +
  theme(legend.position = "none",axis.text=element_text(size=8),
        axis.title=element_text(size=8))+
  #scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_bar(stat="identity")+
  labs(y = "Taxon resolved", x ="")+
  facet_wrap(~ rank,scales = "fixed", nrow = 1)

library(cowplot)
#png("ecotaxspecificity.png", height = 10, width = 20, res = 300, units = "cm")
plot_grid(res2, res1, labels = c('A', 'B'), nrow = 2)
#dev.off()


########################################################

#taxonomic overlap via Venn Diagrams
########################################################
set.seed(20190708)

#subset Protists see further down --> 

#use full lineage
V9=marker.lineages.df[grep("V9_", marker.lineages.df$marker), ]
V7=marker.lineages.df[grep("V7_", marker.lineages.df$marker), ]
V4=marker.lineages.df[grep("V4", marker.lineages.df$marker), ]
V3=marker.lineages.df[grep("V3_", marker.lineages.df$marker), ]
dim(V9)
#47507
dim(V7)
#113436     
dim(V4)
#57102    
dim(V3)
#4065   

V7$marker=as.character(V7$marker)
V9$marker=as.character(V9$marker)
V4$marker=as.character(V4$marker)
V3$marker=as.character(V3$marker)
names(V9)

V9 <- unique( V9[ , 1:9 ] )
V7 <- unique( V7[ , 1:9 ] )
V4 <- unique( V4[ , 1:9 ] )
V3 <- unique( V3[ , 1:9 ] )
dim(V9)
#32893
dim(V7)
#87137     
dim(V4)
#42073    
dim(V3)
#2643    

V4.V7=rbind(V4, V7)
V4.V9=rbind(V4, V9)
V4.V7.V3=rbind(V4, V7, V3)
V4.V9.V3=rbind(V4, V9, V3)
V4.V3=rbind(V4, V3)
V7.V3=rbind(V7, V3)
V9.V3=rbind(V9, V3)

#spec. names
# V7.list=split(V7$Species,V7$marker)
# V9.list=split(V9$Species,V9$marker)
# V4.V7.list=split(V4.V7$Species,V4.V7$marker)
# V4.V9.list=split(V4.V9$Species,V4.V9$marker)
# V4.V3.list=split(V4.V3$Species,V4.V3$marker)
# V7.V3.list=split(V7.V3$Species,V7.V3$marker)
# V9.V3.list=split(V9.V3$Species,V9.V3$marker)

#taxIDs
V7.list=split(as.character(V7$NCBI_taxid),V7$marker)
V9.list=split(as.character(V9$NCBI_taxid),V9$marker)
V4.V7.V3.list=split(as.character(V4.V7.V3$NCBI_taxid),V4.V7.V3$marker)
V4.V9.V3.list=split(as.character(V4.V9.V3$NCBI_taxid),V4.V9.V3$marker)
V4.V3.list=split(as.character(V4.V3$NCBI_taxid),V4.V3$marker)
V7.V3.list=split(as.character(V7.V3$NCBI_taxid),V7.V3$marker)
V9.V3.list=split(as.character(V9.V3$NCBI_taxid),V9.V3$marker)
V4.V7.list=split(as.character(V4.V7$NCBI_taxid),V4.V7$marker)
V4.V9.list=split(as.character(V4.V9$NCBI_taxid),V4.V9$marker)

a<-ggVennDiagram(V9.list, label_alpha = 0) +
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")

b<-ggVennDiagram(V7.list, label_alpha = 0)+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")

c<-ggVennDiagram(V4.V7.V3.list, label_alpha = 0,set_size = 5, label_size = 5, label = "count") +
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")

d<-ggVennDiagram(V4.V9.V3.list, label_alpha = 0,set_size = 5, label_size = 5, label = "count")+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")

e<-ggVennDiagram(V4.V3.list, label_alpha = 0)+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")

f<-ggVennDiagram(V7.V3.list, label_alpha = 0)+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")

g<-ggVennDiagram(V9.V3.list, label_alpha = 0)+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")

################
ggVennDiagram(V4.V7.list, label_alpha = 0) +
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")

ggVennDiagram(V4.V9.list, label_alpha = 0)+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")

################
#compare 18S_allshorts, TARA Oceans marker, SAR_V3 and V4
V7.18Sallshorts=marker.lineages.df[grep("V7_18S_allshorts", marker.lineages.df$marker), ]
V9_1389F_1510R=marker.lineages.df[grep("V9_1389F_1510R", marker.lineages.df$marker), ]
V3.2=marker.lineages.df[grep("V4", marker.lineages.df$marker), ]
V4.2=marker.lineages.df[grep("V3", marker.lineages.df$marker), ]

marker.comb<-rbind(V7.18Sallshorts, V9_1389F_1510R,V4.2, V3.2)
marker.comb$marker=as.character(marker.comb$marker)
marker.comb.list=split(marker.comb$Species,marker.comb$marker)

ggVennDiagram(marker.comb.list, label_alpha = 0,set_size = 5, label_size = 5, label = "count")+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")
ggVennDiagram(marker.comb.list, label_alpha = 0)+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")


c + d + aa +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = 'A')


############################################
#Filter protist groups by phyla

protists<-c("Apicomplexa", "Bacillariophyta", "Cercozoa", "Chlorophyta",
            "Ciliophora", "Discosea", "Endomyxa", "Euglenozoa",
            "Evosea", "Foraminifera", "Haptista", "Picozoa",
            "Prasinodermophyta", "Tubulinea")

euk="Eukaryota"
unwanted<-c("Arthropoda", "Platyhelminthes",
            "Mollusca", "Nematoda", "Annelida",
            "Porifera", "Cnidaria",                   
            "Tardigrada", "Acanthocephala","Ctenophora",                 
            "Chordata", "Nematomorpha", 
            "Bryozoa",                    
            "Echinodermata", "Nemertea","Rotifera",                   
            "Brachiopoda"  ,"Onychophora" ,               
            "Priapulida","Chaetognatha",
            "Hemichordata",
            "Ascomycota", "Basidiomycota", "Mucoromycota", "Zoopagomycota",
            "Olpidiomycota", "Chytridiomycota","Blastocladiomycota",
            "Streptophyta", "Rhodophyta")
unwanted.class=c("Entorrhizomycetes", "Phaeophyceae", "Ulvophyceae")
a.sub<-subset(marker.lineages.df, Kingdom %in% euk)
b.sub<-subset(a.sub, !Phylum %in% unwanted)
c.sub<-subset(b.sub, !Class %in% unwanted.class)

c.subclass<-unique(c.sub$Class)
pclass.sub<-subset(marker.lineages.df, Class %in% c.subclass)
dim(pclass.sub)
#53461    11

classes.clades<-read.table("protist-clades-classes.txt", sep="\t", header=TRUE)
Stramenopiles<-subset(classes.clades, Clade=="Stramenopiles")
Alveolata<-subset(classes.clades, Clade=="Alveolata")
Rhizaria<-subset(classes.clades, Clade=="Rhizaria")
Viridiplantae<-subset(classes.clades, Clade=="Viridiplantae")
others<-subset(classes.clades, Clade=="Discoba" |Clade=="Haptista" |
                 Clade=="incertae sedis" |Clade=="Eukaryota" |
                 Clade=="Opisthokonta"  )
Haptista<-subset(classes.clades, Clade=="Haptista")
Foraminifera<-subset(classes.clades, Clade2=="Foraminifera")

#sequence lengths per clade
marker.lineages.df$sequence.length=nchar(marker.lineages.df$sequence)

#create subsets of different clades
Stramenopiles.sub<-subset(marker.lineages.df, Class %in% Stramenopiles$Class)
Alveolata.sub<-subset(marker.lineages.df, Class %in% Alveolata$Class)
Rhizaria.sub<-subset(marker.lineages.df, Class %in% Rhizaria$Class)
Viridiplantae.sub<-subset(marker.lineages.df, Class %in% Viridiplantae$Class)
others.sub<-subset(marker.lineages.df, Class %in% others$Class)
Haptista.sub<-subset(marker.lineages.df, Class %in% Haptista$Class)
Chlorophyta.sub<-subset(marker.lineages.df, Class %in% others$Class)
Foraminifera.sub<-subset(marker.lineages.df, Class %in% Foraminifera$Class)

#######################################
#Stramenopiles
#######################################
seqlenmarker.s <-Stramenopiles.sub %>% 
  group_by(marker) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

ggplot(seqlenmarker.s, aes(x=marker, y=sequence.length, fill=marker)) + 
  theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_violin(trim=TRUE)

#classes within clade
seqlenmarker.sc <-Stramenopiles.sub %>% 
  group_by(marker, Class) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

seqlenmarker.sc$Class<- factor(seqlenmarker.sc$Class, levels=rev(seqlenmarker.sc(seqlenmarker.sc$Class)))

#scatter plot Stramenopiles per class
ggplot(seqlenmarker.sc, aes(x=Class, y=sequence.length, fill=marker)) + 
  coord_flip()+  
  theme(legend.position = "none", axis.text = element_text(size = 12))+
  geom_point(size=.75, alpha=.25)+
  scale_x_discrete(limits=rev)+
  facet_wrap(~ marker,scales = "fixed", nrow = 2)

#violin plot Stramenopiles
ggplot(seqlenmarker.sc, aes(x=marker, y=sequence.length, fill=marker)) + 
  theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_violin(trim=TRUE)

#######################################
#Alveolata
#######################################
seqlenmarker.a <-Alveolata.sub %>% 
  group_by(marker) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

#violin plot Alveolata
ggplot(seqlenmarker.a, aes(x=marker, y=sequence.length, fill=marker)) + 
  theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_violin(trim=TRUE)

#classes within clade
seqlenmarker.ac <-Alveolata.sub %>% 
  group_by(marker, Class) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

#scatter plot Alveolata
ggplot(seqlenmarker.ac, aes(x=Class, y=sequence.length, fill=marker)) + 
  theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none", axis.text = element_text(size = 10))+
  #scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_point(size=.75, alpha=.25)+
  facet_wrap(~ marker,scales = "fixed", nrow = 3)

#######################################
#Rhizaria
#######################################
seqlenmarker.r <-Rhizaria.sub %>% 
  group_by(marker) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

#violin plot
ggplot(seqlenmarker.r, aes(x=marker, y=sequence.length, fill=marker)) + 
  theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_violin(trim=TRUE)

#classes within clade
seqlenmarker.rc <-Rhizaria.sub %>% 
  group_by(marker, Class) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

ggplot(seqlenmarker.rc, aes(x=Class, y=sequence.length, fill=marker)) + 
  coord_flip()+  
  theme(legend.position = "none", axis.text = element_text(size = 10))+
  geom_point(size=.75, alpha=.25)+
  scale_x_discrete(limits=rev)+
  facet_wrap(~ marker,scales = "fixed", nrow = 3)

#######################################
#Viridiplantae
#######################################
seqlenmarker.v <-Viridiplantae.sub %>% 
  group_by(marker) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

ggplot(seqlenmarker.v, aes(x=marker, y=sequence.length, fill=marker)) + 
  theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_violin(trim=TRUE)


#classes within clade
seqlenmarker.vc <-Viridiplantae.sub %>% 
  group_by(marker, Class) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

ggplot(seqlenmarker.vc, aes(x=Class, y=sequence.length, fill=marker)) + 
  coord_flip()+  
  theme(legend.position = "none", axis.text = element_text(size = 10))+
  #scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_point(size=.75, alpha=.25)+
  scale_x_discrete(limits=rev)+
  facet_wrap(~ marker,scales = "fixed", nrow = 3)

#######################################
#Haptista
#######################################
seqlenmarker.h <-Haptista.sub %>% 
  group_by(marker) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

ggplot(seqlenmarker.h, aes(x=marker, y=sequence.length, fill=marker)) + 
  theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_violin(trim=TRUE)

#classes within clade
seqlenmarker.hc <-Haptista.sub %>% 
  group_by(marker, Class) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

ggplot(seqlenmarker.hc, aes(x=Class, y=sequence.length, fill=marker)) + 
  coord_flip()+  
  theme(legend.position = "none", axis.text = element_text(size = 10))+
  geom_point(size=.75, alpha=.25)+
  scale_x_discrete(limits=rev)+
  facet_wrap(~ marker,scales = "fixed", nrow = 3)

#######################################
#others
#######################################
seqlenmarker.o <-others.sub %>% 
  group_by(marker) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

ggplot(seqlenmarker.o, aes(x=marker, y=sequence.length, fill=marker)) + 
  theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_violin(trim=TRUE)

#classes within clade
seqlenmarker.oc <-others.sub %>% 
  group_by(marker, Class) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

ggplot(seqlenmarker.oc, aes(x=Class, y=sequence.length, fill=marker)) + 
  #theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none", axis.text = element_text(size = 12))+
  #scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_point(size=.75, alpha=.25)+
  scale_x_discrete(limits=rev)+
  facet_wrap(~ marker,scales = "fixed", nrow = 3)

#################################################
#Metazoa
#################################################
seqlenmarker.m <-me.sub %>% 
  group_by(marker) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

ggplot(seqlenmarker.m, aes(x=marker, y=sequence.length, fill=marker)) + 
  theme_minimal()+ 
  coord_flip()+  
  theme(legend.position = "none")+
  scale_fill_manual(values=c("#A16B56","#E76F51", "#F4A261", "#E9C46A" ,"#f4debe","#bfa082","#847b45", "#264653","#8d96a3","#2A9D8F"))+
  geom_violin(trim=TRUE)

#classes within clade
seqlenmarker.mp <-me.sub %>% 
  group_by(marker, Phylum) %>%
  distinct(NCBI_taxid, .keep_all = TRUE)

ggplot(seqlenmarker.mp, aes(x=Phylum, y=sequence.length, fill=marker)) + 
  coord_flip()+  
  theme(legend.position = "none", axis.text = element_text(size = 12))+
  geom_point(size=.75, alpha=.25)+
  scale_x_discrete(limits=rev)+
  facet_wrap(~ marker,scales = "fixed", nrow = 3)

#################################################
#number of genera per clade per marker
no.g.s<-Stramenopiles.sub %>%
  group_by(marker) %>%
  summarise(n.class = length(unique(Genus)))
no.g.a<-Alveolata.sub %>%
  group_by(marker) %>%
  summarise(n.genus = length(unique(Genus)))
no.g.r<-Rhizaria.sub %>%
  group_by(marker) %>%
  summarise(n.genus = length(unique(Genus)))
no.g.v<-Viridiplantae.sub %>%
  group_by(marker) %>%
  summarise(n.genus = length(unique(Genus)))
no.g.s<-others.sub %>%
  group_by(marker) %>%
  summarise(n.class = length(unique(Genus)))
no.g.h<-Haptista.sub %>%
  group_by(marker) %>%
  summarise(n.genus = length(unique(Genus)))


#################################################
#Venn diagram of classes within clades

#subsets of Protists
Stramenopiles.sub<-subset(marker.lineages.df, Class %in% Stramenopiles$Class)
Alveolata.sub<-subset(marker.lineages.df, Class %in% Alveolata$Class)
Rhizaria.sub<-subset(marker.lineages.df, Class %in% Rhizaria$Class)
Viridiplantae.sub<-subset(marker.lineages.df, Class %in% Viridiplantae$Class)
others.sub<-subset(marker.lineages.df, Class %in% others$Class)
Haptista.sub<-subset(marker.lineages.df, Class %in% Haptista$Class)

#choose subset and run from L775-L850
prot.sub=Rhizaria.sub

#use full lineage
V9=prot.sub[grep("V9_", prot.sub$marker), ]
V7=prot.sub[grep("V7_", prot.sub$marker), ]
V4=prot.sub[grep("V4", prot.sub$marker), ]
V3=prot.sub[grep("V3_", prot.sub$marker), ]

V7$marker=as.character(V7$marker)
V9$marker=as.character(V9$marker)
V4$marker=as.character(V4$marker)
V3$marker=as.character(V3$marker)

V9 <- unique( V9[ , 1:9 ] )
V7 <- unique( V7[ , 1:9 ] )
V4 <- unique( V4[ , 1:9 ] )
V3 <- unique( V3[ , 1:9 ] )

V4.V7=rbind(V4, V7)
V4.V9=rbind(V4, V9)
V4.V3=rbind(V4, V3)
V7.V3=rbind(V7, V3)
V9.V3=rbind(V9, V3)

V7.list=split(as.character(V7$NCBI_taxid),V7$marker)
V9.list=split(as.character(V9$NCBI_taxid),V9$marker)
V4.V7.list=split(as.character(V4.V7$NCBI_taxid),V4.V7$marker)
V4.V9.list=split(as.character(V4.V9$NCBI_taxid),V4.V9$marker)
V4.V3.list=split(as.character(V4.V3$NCBI_taxid),V4.V3$marker)
V7.V3.list=split(as.character(V7.V3$NCBI_taxid),V7.V3$marker)
V9.V3.list=split(as.character(V9.V3$NCBI_taxid),V9.V3$marker)

prot.sub$marker=as.character(prot.sub$marker)
V7.18Sallshorts=subset(prot.sub, marker=="V7_18S_allshorts")
V9_1389F_1510R=subset(prot.sub, marker=="V9_1389F_1510R")
V4=subset(prot.sub, marker=="V4")
SAR_V3=subset(prot.sub, marker=="V3_SAR_V3")

marker.comb<-rbind(V7.18Sallshorts, V9_1389F_1510R,V4, SAR_V3)
marker.comb$marker=as.character(marker.comb$marker)
marker.comb.list=split(marker.comb$Species,marker.comb$marker)

a<-ggVennDiagram(V9.list, label_alpha = 0, set_size = 3, label_size = 3, label = "count") +
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")+
  theme(text = element_text(size = 11))

b<-ggVennDiagram(V7.list, label_alpha = 0, set_size = 3, label_size = 3, label = "count")+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")+
  theme(text = element_text(size = 11))

c<-ggVennDiagram(V4.V7.list, label_alpha = 0, set_size = 3, label_size = 3, label = "count") +
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")+
  theme(text = element_text(size = 11))

d<-ggVennDiagram(V4.V9.list, label_alpha = 0, set_size = 3, label_size = 3, label = "count")+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")+
  theme(text = element_text(size = 11))

e<-ggVennDiagram(V4.V3.list, label_alpha = 0, set_size = 3, label_size = 3, label = "count")+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")+
  theme(text = element_text(size = 11))

f<-ggVennDiagram(V7.V3.list, label_alpha = 0, set_size = 3, label_size = 3, label = "count")+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")+
  theme(text = element_text(size = 11))

g<-ggVennDiagram(V9.V3.list, label_alpha = 0, set_size = 3, label_size = 3, label = "count")+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")+
  theme(text = element_text(size = 11))

h<-ggVennDiagram(marker.comb.list, label_alpha = 0,set_size = 3, label_size = 3,label = "count")+
  ggplot2::scale_fill_gradient(low="cornflowerblue",high = "#E9C46A")+
  theme(text = element_text(size = 11))

#plot
c + d + h+
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = 'A')
