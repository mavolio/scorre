

#####
# Create Phylo-ring to show patterns on tips
#####

#load packages:
library(ggtree)
library(ggplot2)
library(V.PhyloMaker)
library(stringr)
library(plyr)
library(ape)
library(scico)

#load species data:
species.data<-read.table("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/Species_DCiDiff.csv", header=T, sep=",")
species.data$species_matched <- revalue(species.data$species_matched, c("Aronia x"="Aronia x prunifolia"))
species.data<-subset(species.data, trt_type2=="N") #subset one treatment to make it easier (here N addition).

#create table for the tree:
spp<-as.data.frame(unique(species.data$species_matched))
names(spp)[1]<-paste("species")
spp$genus<-word(spp$species, 1)

#load families for species and rearrange to create necessary fields for the phylogenies:
fam<-read.table("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/species_families_2021.csv", header=T, sep=",", fill = TRUE)
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
df<-unique(species.data[,c(1,5)])
rownames(df)<-df$species_matched
df$species_matched<-NULL

#Plot tree:
p <- ggtree(scorre.tree$scenario.3, layout="circular", size=1.5, branch.length="none")+
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=17.5),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

#Add heatmap for treatment:
p <- gheatmap(p, df, offset=0.03, width=.05,colnames = F,
              colnames_angle=90, colnames_offset_y = .25) +
  scale_fill_scico(palette = "vik", limits = c(-1, 1) * max(abs(df$ave_diff)))

#plot:
png("phylo_ring2.png",
    res=300,height=25,width=25,units="in"); 
p
dev.off()
