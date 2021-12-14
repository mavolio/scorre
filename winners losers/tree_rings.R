

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
#my.wd <- "E:/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/"
my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/"

#load data.
species.data<-read.table(paste(my.wd,"Species_DCiDiff_Dec2021.csv",sep=""), header=T, sep=",")

#create table for the tree:
spp<-as.data.frame(unique(species.data$species_matched))
names(spp)[1]<-paste("species")
spp$genus<-word(spp$species, 1)

non.vascular <-  c("Andreaea obovata",            "Anthelia juratzkana" ,       "Aulacomnium turgidum",       
                   "Barbilophozia hatcheri",      "Barbilophozia kunzeana" ,     "Blepharostoma trichophyllum",
                   "Brachythecium albicans",      "Bryum arcticum"   ,           "Bryum pseudotriquetrum",     
                   "Campylium stellatum",         "Cyrtomnium hymenophyllum" ,   "Dicranoweisia crispula",     
                   "Dicranum brevifolium",        "Dicranum elongatum"  ,        "Dicranum fuscescens",        
                   "Dicranum groenlandicum",      "Dicranum scoparium" ,         "Distichium capillaceum",     
                   "Ditrichum flexicaule",        "Gymnomitrion concinnatum" ,   "Hamatocaulis vernicosus",    
                   "Homalothecium pinnatifidum",  "Hylocomium splendens",        "Hypnum cupressiforme",       
                   "Hypnum hamulosum",            "Isopterygiopsis pulchella",   "Kiaeria starkei",            
                   "Leiocolea heterocolpos",      "Marchantia polymorpha",       "Marsupella brevissima",      
                   "Meesia uliginosa",            "Myurella tenerrima",          "Oncophorus virens",         
                   "Oncophorus wahlenbergii",     "Pleurozium schreberi",        "Pogonatum urnigerum" ,       
                   "Pohlia cruda" ,               "Pohlia nutans",               "Polytrichastrum alpinum",    
                   "Polytrichum juniperinum",     "Polytrichum piliferum",       "Polytrichum strictum",       
                   "Preissia quadrata",           "Ptilidium ciliare",           "Racomitrium lanuginosum",    
                   "Rhytidium rugosum",           "Saelania glaucescens",        "Sanionia uncinata",          
                   "Schistidium apocarpum",       "Syntrichia ruralis",          "Tomentypnum nitens",         
                   "Tortella tortuosa",           "Tritomaria quinquedentata",   "Nephroma arcticum" , "Unknown NA",
                   "Campylopus flexuosus",        "Hypnum jutlandicum",          "Plagiothecium undulatum",    
                   "Polytrichum commune",         "Pseudoscleropodium purum",    "Rhytidiadelphus loreus",
                   "Rhytidiadelphus triquetrus",  "Thuidium tamariscinum")
spp <- spp[!spp$species %in% non.vascular, ] #remove

#load families for species and rearrange to create necessary fields for the phylogenies:
fam<-read.table(paste(my.wd, "species_families_trees_2021.csv",sep=""), header=T, sep=",", fill = TRUE)
names(fam)[1]<-paste("species")
names(fam)[2]<-paste("family")
spp<-merge(spp, fam, by="species", all.x=T)
spp$species.relative <-  rep("",length(nrow(spp)))
spp$genus.relative <-  rep("",length(nrow(spp)))
spp$tree.non.tree<-NULL

#unify family names:
spp$family[spp$family=="Compositae"]<-"Asteraceae"
spp$family[spp$family=="Leguminosae"]<-"Fabaceae"
spp$family[spp$family=="Viburnaceae"]<-"Adoxaceae"
spp$family[spp$family=="Polypodiaceae"]<-"Dryopteridaceae"
spp$family[spp$genus=="Blechnum"]<-"Blechnaceae"
spp$family[spp$genus=="Onoclea"]<-"Onocleaceae"
spp$family[spp$genus=="Thelypteris"]<-"Thelypteridaceae"

#get phylogenetic tree:
scorre.tree <- phylo.maker(sp.list = spp, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S3")

#remove from the original table non-vascular plants that couldn't be added to the tree:
species.data$species_matched<-gsub(" ", "_", species.data$species_matched) #unify nomenclature
in.data.not.tree <- setdiff(unique(species.data$species_matched), scorre.tree$scenario.3$tip.label)
species.data <- species.data[-which(species.data$species_matched %in% in.data.not.tree),] #only works if some species from the data are not in the tree

#save tree and subset table:
write.tree(scorre.tree$scenario.3, paste(my.wd, "scorre.tree.win.los.tre.dec2021", sep=""))
write.table(species.data, paste(my.wd, "Species_DCiDiff_filtered_Dec2021.csv", sep=""))







#####
# Create big phylogenetic ring (all treatments)
#####

#subset table with all treatments:
trt<-subset(species.data, trt_type2=="all mult")

#prune tree to get only
tree<-read.tree(paste(my.wd, "scorre.tree.win.los.tre", sep="")) #load tree
tree<-keep.tip(tree, unique(trt$species_matched)) #all species are in the tree

#get data frame with DCi:
final.trt<-trt[,c(1,4)]
rownames(final.trt)<-final.trt$species_matched
final.trt$species_matched<-NULL

#Plot tree:
p <- ggtree(tree, layout="circular", size=0.5, branch.length="none")+
  #geom_tiplab() +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=17.5),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

#Add heatmap for treatment:
p <- gheatmap(p, final.trt, offset=0.03, width=.05, colnames = F,
              colnames_angle=90) +
  scale_fill_scico(palette = "vik", limits = c(-1, 1) * max(abs(final.trt), na.rm = T), na.value="white")

#plot:
png("phylo_ring_all_mult.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()







#####
# Create big phylogenetic ring (CO2 + CO2 & others)
#####

#subset table with all treatments:
trt<-subset(species.data, trt_type2=="co2" | trt_type2=="co2_other")

#prune tree to get only:
tree<-read.tree(paste(my.wd, "scorre.tree.win.los.tre", sep="")) #load tree
tree<-keep.tip(tree, unique(trt$species_matched)) #all species are in the tree

#rearrange original table:
trts<-levels(as.factor(as.character(trt$trt_type2)))
final.trt<-as.data.frame(unique(trt$species_matched))
names(final.trt)[1]<-paste("species_matched")

for (i in 1:length(trts))
{
  sub.trt<-subset(trt, trt_type2==trts[i])[,c(1,4)]
  final.trt<-merge(final.trt, sub.trt, by="species_matched", all.x=T)
  colnames(final.trt)[i+1]<-trts[i]
}
rownames(final.trt)<-final.trt$species_matched
final.trt$species_matched<-NULL

#Plot tree:
p <- ggtree(tree, layout="circular", size=0.5, branch.length="none")+
  #geom_tiplab() +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=17.5),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

#Add heatmap for treatment:
p <- gheatmap(p, final.trt, offset=0.03, width=.05, colnames = F,
              colnames_angle=90) +
  scale_fill_scico(palette = "vik", limits = c(-1, 1) * max(abs(final.trt), na.rm = T), na.value="white")

#plot:
png("phylo_ring_co2.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()






#####
# Create big phylogenetic ring (drought + irrigation)
#####

#subset table with all treatments:
trt<-subset(species.data, trt_type2=="drought" | trt_type2=="drt_other" | trt_type2=="irrigation" | trt_type2=="irg_other")

#prune tree to get only:
tree<-read.tree(paste(my.wd, "scorre.tree.win.los.tre", sep="")) #load tree
tree<-keep.tip(tree, unique(trt$species_matched)) #all species are in the tree

#rearrange original table:
trts<-levels(as.factor(as.character(trt$trt_type2)))
final.trt<-as.data.frame(unique(trt$species_matched))
names(final.trt)[1]<-paste("species_matched")

for (i in 1:length(trts))
{
  sub.trt<-subset(trt, trt_type2==trts[i])[,c(1,4)]
  final.trt<-merge(final.trt, sub.trt, by="species_matched", all.x=T)
  colnames(final.trt)[i+1]<-trts[i]
}
rownames(final.trt)<-final.trt$species_matched
final.trt$species_matched<-NULL

#Plot tree:
p <- ggtree(tree, layout="circular", size=0.5, branch.length="none")+
  #geom_tiplab() +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=17.5),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

#Add heatmap for treatment:
p <- gheatmap(p, final.trt, offset=0.03, width=.1, colnames = F,
              colnames_angle=90) +
  scale_fill_scico(palette = "vik", limits = c(-1, 1) * max(abs(final.trt), na.rm = T), na.value="white")

#plot:
png("phylo_ring_drought_irrigation.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()






#####
# Create big phylogenetic ring (disturbance + herb removal)
#####

#subset table with all treatments:
trt<-subset(species.data, trt_type2=="disturbance" | trt_type2=="dist_other" | trt_type2=="herb_removal" | trt_type2=="herb_rem_other")

#prune tree to get only:
tree<-read.tree(paste(my.wd, "scorre.tree.win.los.tre", sep="")) #load tree
tree<-keep.tip(tree, unique(trt$species_matched)) #all species are in the tree

#rearrange original table:
trts<-levels(as.factor(as.character(trt$trt_type2)))
final.trt<-as.data.frame(unique(trt$species_matched))
names(final.trt)[1]<-paste("species_matched")

for (i in 1:length(trts))
{
  sub.trt<-subset(trt, trt_type2==trts[i])[,c(1,4)]
  final.trt<-merge(final.trt, sub.trt, by="species_matched", all.x=T)
  colnames(final.trt)[i+1]<-trts[i]
}
rownames(final.trt)<-final.trt$species_matched
final.trt$species_matched<-NULL

#Plot tree:
p <- ggtree(tree, layout="circular", size=0.5, branch.length="none")+
  #geom_tiplab() +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=17.5),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

#Add heatmap for treatment:
p <- gheatmap(p, final.trt, offset=0.03, width=.1, colnames = F,
              colnames_angle=90) +
  scale_fill_scico(palette = "vik", limits = c(-1, 1) * max(abs(final.trt), na.rm = T), na.value="white")

#plot:
png("phylo_ring_disturbance.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()








#####
# Create big phylogenetic ring (temperature)
#####

#subset table with all treatments:
trt<-subset(species.data, trt_type2=="temp" | trt_type2=="temp_other")

#prune tree to get only:
tree<-read.tree(paste(my.wd, "scorre.tree.win.los.tre", sep="")) #load tree
tree<-keep.tip(tree, unique(trt$species_matched)) #all species are in the tree

#rearrange original table:
trts<-levels(as.factor(as.character(trt$trt_type2)))
final.trt<-as.data.frame(unique(trt$species_matched))
names(final.trt)[1]<-paste("species_matched")

for (i in 1:length(trts))
{
  sub.trt<-subset(trt, trt_type2==trts[i])[,c(1,4)]
  final.trt<-merge(final.trt, sub.trt, by="species_matched", all.x=T)
  colnames(final.trt)[i+1]<-trts[i]
}
rownames(final.trt)<-final.trt$species_matched
final.trt$species_matched<-NULL

#Plot tree:
p <- ggtree(tree, layout="circular", size=0.5, branch.length="none")+
  #geom_tiplab() +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=17.5),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

#Add heatmap for treatment:
p <- gheatmap(p, final.trt, offset=0.03, width=.05, colnames = F,
              colnames_angle=90) +
  scale_fill_scico(palette = "vik", limits = c(-1, 1) * max(abs(final.trt), na.rm = T), na.value="white")

#plot:
png("phylo_ring_temp.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()





#####
# Create big phylogenetic ring (nutrients)
#####

#subset table with all treatments:
trt<-subset(species.data, trt_type2=="nutrients" | trt_type2=="nuts_other")

#prune tree to get only:
tree<-read.tree(paste(my.wd, "scorre.tree.win.los.tre", sep="")) #load tree
tree<-keep.tip(tree, unique(trt$species_matched)) #all species are in the tree

#rearrange original table:
trts<-levels(as.factor(as.character(trt$trt_type2)))
final.trt<-as.data.frame(unique(trt$species_matched))
names(final.trt)[1]<-paste("species_matched")

for (i in 1:length(trts))
{
  sub.trt<-subset(trt, trt_type2==trts[i])[,c(1,4)]
  final.trt<-merge(final.trt, sub.trt, by="species_matched", all.x=T)
  colnames(final.trt)[i+1]<-trts[i]
}
rownames(final.trt)<-final.trt$species_matched
final.trt$species_matched<-NULL

#Plot tree:
p <- ggtree(tree, layout="circular", size=0.5, branch.length="none")+
  #geom_tiplab() +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=17.5),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

#Add heatmap for treatment:
p <- gheatmap(p, final.trt, offset=0.03, width=.05, colnames = F,
              colnames_angle=90) +
  scale_fill_scico(palette = "vik", limits = c(-1, 1) * max(abs(final.trt), na.rm = T), na.value="white")

#plot:
png("phylo_ring_nuts.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()