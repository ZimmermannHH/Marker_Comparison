library(tidyr)   #v. 1.3.0
library(dplyr)   #v. 1.1.2
library(ggimage) #v. 
library(ggtree)  #v. 3.9.0.001
library(TDbook)  #v. 0.0.6
#remotes::install_github("arendsee/phylostratr")
library(phylostratr) #v. 0.2.1


#listing all Eukaryote phyla and classes in ncbi taxonomy using Taxonkit (on unix server)
#taxonkit list --ids 2759 | taxonkit reformat -I 1 -F -P -f "{k}|{p}|{c}|{o}|{f}|{g}|{s}|{t}" | csvtk pretty -H -t > ncbi_eukaryota_reformatted.txt

#import output from taxonkit into R
db<-read.table("ncbi_eukaryota_reformatted.txt", header=TRUE, sep="\t" ,quote = "", fill=TRUE,
               stringsAsFactors = FALSE, encoding="UTF-8")

db2=separate(data=db,col = lineage, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep = "\\|")
head(db2)
db2[2]=NULL
db2$Phylum=gsub("p__", "",db2$Phylum)
db2$Class=gsub("c__", "",db2$Class)
db2$Order=gsub("o__", "",db2$Order)
db2$Family=gsub("f__", "",db2$Family)
db2$Genus=gsub("g__", "",db2$Genus)
db2$Species=gsub("s__", "",db2$Species)
db2$Strain=gsub("t__", "",db2$Strain)

classes<-read.table("ncbi_eukaryota_classes.txt", header=FALSE, sep="\t" )
classes[3]=NULL
classes$V2=gsub("cellular organisms;Eukaryota;", "", classes$V2)


##############################################################################
#import subset of output which only contains protists
clades.classes<-read.table("ncbi_eukaryota_classes_cladeinfo.txt", header=TRUE, sep="\t" )
#Foraminifera --> only Globothalamea as class in NCBI taxonomy

taxids<-clades.classes$taxid

classes_label<-as.data.frame(cbind(clades.classes$taxid, clades.classes$Class))
names(classes_label)=c( "ncbi_taxid","Newick_label")

#build phylogenetic tree based on NCBI-Taxonomy
classes.tree<-ncbi_tree(taxids)

#replace taxids in tip.label by class name 
classes.tree$tip.label=clades.classes$Class

set.seed(123)


#prepare heatmap
marker.lineages.df<-read.table("marker.lineages.df.txt", sep="\t", header=TRUE)
#replace primer names
marker.lineages.df$marker=gsub("V9_Amaral1", "V9_1380F-1510R", marker.lineages.df$marker)
marker.lineages.df$marker=gsub("V9_Piredda", "V9_1388F-1510R", marker.lineages.df$marker)
marker.lineages.df$marker=gsub("V9_Stoeck", "V9_1391F- EukB", marker.lineages.df$marker)

protist.classes.sub<-subset(marker.lineages.df, Class %in% clades.classes$Class)
protist.classes.sub<-subset(marker.lineages.df, marker!="V7_18S_allshorts_modF")

heatmap.in.genus<-protist.classes.sub %>% 
  group_by(marker) %>% 
  distinct(Genus, .keep_all = TRUE) %>% 
  count(Class)

#prepare matrix as input for gheatmap
hm<-heatmap.in.genus %>%
  pivot_wider(names_from = marker, values_from = n, values_fill = 0)

hm=as.data.frame(hm)
row.names(hm)=hm$Class
hm[1]=NULL
hm.ma=as.matrix(hm)

#to highlight Classes for which no genera etc. are retrieved, use NA instead of 0 
#then the NAs will be excluded from the color gradient and better visible (as seen in
#heatmaps above)
hm.ma2 <- replace(hm.ma, hm.ma==0, NA) 

#prepare barplot: genera per class
db2.protists<-subset(db2, Class %in% clades.classes$Class)
dim(db2.protists)
#41563 entries

no.ncbi.genus<-db2.protists %>%
  group_by(Class) %>%
  summarise(n.genus = length(unique(Genus)))


#plot tree, heatmap and barplot

p<-ggtree(classes.tree) +  xlim(-.1, 8)+ 
  geom_tiplab(size=3, align=TRUE, linesize=.5) 

pdf("phylo-class.heatmap.gen_new.pdf", width = 10, height = 12) #offset 3.2 for pdf
#png("phylo-class.heatmap-genus2.png", width = 20, height = 17, res=300, units = "cm") #offset 4 for png
gheatmap(p, hm.ma2, offset=3.2, width=0.6,
         font.size=2, colnames_angle=-45, hjust=0,
         colnames=TRUE, legend_title="Genus") + 
  #geom_tiplab(as_ylab=TRUE)+
  scale_x_ggtree() + 
  scale_y_continuous(expand = c(0, .3))+
  scale_fill_gradient2(low="white", mid="#2A9D8F", high = "cornflowerblue", na.value = "grey50",
                       midpoint = 70, space = "Lab")+
  geom_facet(panel = "Trait", data = no.ncbi.genus, geom = geom_col, 
             aes(x = n.genus), orientation = 'y', width = .8)+
  scale_x_continuous(expand = c(0, .3))+
  theme_tree2(legend.position=c(.95, .85))
dev.off()

