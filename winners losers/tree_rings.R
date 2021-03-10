

#####
# Create Phylo-ring to show patterns on tips
#####

#load packages:
#for ggtree if you have r version 4 or above
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("ggtree")

library(ggtree)
library(ggplot2)
library(stringr)
library(plyr)
library(ape)
library(scico)

#for V.PhyloMaker
library(devtools)
#devtools::install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)

#set directory:
#my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/"
my.wd <- "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"
#my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"

#load data.
species.data<-read.table(paste(my.wd,"Species_DCiDiff_newtrts.csv",sep=""), header=T, sep=",")

#create table for the tree:
spp<-as.data.frame(unique(species.data$species_matched))
names(spp)[1]<-paste("species")
spp$genus<-word(spp$species, 1)

#load families for species and rearrange to create necessary fields for the phylogenies:
fam<-read.table(paste(my.wd, "species_families_2021.csv",sep=""), header=T, sep=",", fill = TRUE)
names(fam)[1]<-paste("species")
names(fam)[2]<-paste("family")
spp<-merge(spp, fam, by="species", all.x=T)
spp$species.relative <-  rep("",length(nrow(spp)))
spp$genus.relative <-  rep("",length(nrow(spp)))

#unify family names:
spp$family[spp$family=="Compositae"]<-"Asteraceae"
spp$family[spp$family=="Leguminosae"]<-"Fabaceae"

#get phylogenetic tree:
scorre.tree <- phylo.maker(sp.list = spp, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S3")

#remove from the original table non-vascular plants that couldn't be added to the tree:
species.data$species_matched<-gsub(" ", "_", species.data$species_matched) #unify nomenclature
in.data.not.tree <- setdiff(unique(species.data$species_matched), scorre.tree$scenario.3$tip.label)
species.data <- species.data[-which(species.data$species_matched %in% in.data.not.tree),] #only works if some species from the data are not in the tree

#rearrange original table:
trts<-levels(species.data$trt_type2)
final.trt<-as.data.frame(unique(species.data$species_matched))
names(final.trt)[1]<-paste("species_matched")

for (i in 1:length(trts))
{
  sub.trt<-subset(species.data, trt_type2==trts[i])[,c(1,4)]
  final.trt<-merge(final.trt, sub.trt, by="species_matched", all.x=T)
  colnames(final.trt)[i+1]<-trts[i]
}
rownames(final.trt)<-final.trt$species_matched
final.trt$species_matched<-NULL

#Plot tree:
p <- ggtree(scorre.tree$scenario.3, size=1.5, branch.length="none")+
  #geom_tiplab() +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=17.5),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

#Add heatmap for treatment:
p <- gheatmap(p, final.trt, offset=0.03, width=.5, colnames = T,
              colnames_angle=90) +
  scale_fill_scico(palette = "vik", limits = c(-1, 1) * max(abs(final.trt), na.rm = T), na.value="white")

#plot:
png("phylo_ring.png",
    res=300,height=30,width=20,units="in"); 
p
dev.off()
