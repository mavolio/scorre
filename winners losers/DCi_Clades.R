
#load libraries:
library(devtools)
devtools::install_github("GuangchuangYu/ggtree")

library(ggtree)
library(ggplot2)
library(stringr)
library(plyr)
library(ape)
library(scico)
library(phytools)


#set directory.
my.wd<-"/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"
my.wd<-"E:\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\"
my.wd<-"C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\"

###
# Load data
###

#load species data (filtered after removing mosses and species missing from the phylogeny):
species.data<-read.csv(paste(my.wd,"Species_DCiDiff_Nov2022.csv",sep=""))
species.data$species_matched<-gsub(" ", "_", species.data$species_matched) #unify nomenclature

#load phylogenetic tree:
tree<-read.tree(paste(my.wd, "scorre.tree.win.los.tre.Dec2021", sep=""))
tree<-read.tree(paste(my.wd, "scorre.phylo.tree.S3.tre", sep="")) #why not using the original tree? I detected that some species are missing from the other.

###
# Filter by treatment == all mult
###

dat<-subset(species.data, trt_type2=="all mult")[,c(1,2)] #select "all mult" treatment from original data
#dat<-aggregate(dat[, 2], list(dat$species_matched), mean, na.rm=T) #get mean DCi value per species - this is not necessary - is already the average.
dat<-dat[dat$species_matched %in% tree$tip.label, ] #make sure all species in the data are on the tree
rownames(dat)<-dat$species_matched #set species names as rownames
dat$species_matched<-NULL #and delete column with species names

#prune tree:
tree2<-keep.tip(tree, rownames(dat))

###
# Run function to calculate if each node has significantly higher or lower mean DCi than
# expected if phylogenetic relationships were at random (laod function in the "DCi_nodes_scorre.R" script)
###

res<-node.mean(tree2, dat, 999)
write.table(res, paste(my.wd, "res_phylo_all_mult_Nov22.csv", sep="")) #save the result
#res2<-subset(res, P_value<0.01) #this would tell you what nodes are significant with alpha < 0.01
#tips(tree2, 1543) #and this would tell you what species are found in that clade

###
# Clean-up result to highlight nodes in the tree
###

significant<-res #create a copy of the main result
significant$P_value[significant$P_value>0.05]<-NA #replace non-significant with NA
significant$P_value[significant$SR<3]<-NA #assign NA to nodes with 2 or less species

significant$P_value[1]<-NA #set the first node (the root node) to NA
significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
significant<-c(rep(NA, length(tree2$tip.label)), significant$P_value) #merge "tip nodes" with "inner" tree nodes
significant<-as.factor(significant) #convert values into factors

significant<-factor(significant, levels = c("neg.05", "pos.05" )) #change order of factors



###
# Assign families to nodes
###

#load family data and clean-up:
fam<-read.table(paste(my.wd, "species_families_trees_2021.csv",sep=""), header=T, sep=",", fill = TRUE)[,-c(3)] #load data
fam$species_matched<-gsub(" ", "_", fam$species_matched) #adapt species nomenclature
fam$family[fam$family=="Compositae"]<-"Asteraceae" #replace family name
fam$family[fam$family=="Leguminosae"]<-"Fabaceae" #replace family name
fam$family[fam$family=="Viburnaceae"]<-"Adoxaceae" #replace family name
fam<-fam[which(fam$species_matched %in% rownames(dat)),] #subset only species included in our treatment

#get table with ranked families based on their number of species:
famf<-as.data.frame(table(fam$family)) #create dataframe
famf<-famf[order(-famf$Freq),] #order it
row.names(famf) <- NULL #get rid of rownames

#create vector with unique names of species:
list.r <- unique(fam$family)

#get nodes for families:
list.nod<-NULL
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$species_matched) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- findMRCA(tree2, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)


###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=35)) #select the top 15 families with 5 for more species

#Plot tree 
#remember to change angle = "auto" everywheree to avoid overlap in names. Consider also unifying "barsize" (to 0.1, for example):
p <- 
  ggtree(tree2, layout="circular", size=0.5)+ # build circular tree
  geom_point(aes(colour=as.factor(significant)), size=2, alpha=1, show.legend = TRUE) + # highlight nodes
  scale_colour_manual(values=c("red", "deepskyblue"), labels=c("Lower DCi", "Higher DCi"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  geom_cladelabel(node=subset(famf, Var1=="Poaceae")$num, label="Poaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asteraceae")$num, label="Asteraceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Fabaceae")$num, label="Fabaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 15) +
  geom_cladelabel(node=subset(famf, Var1=="Cyperaceae")$num, label="Cyperaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Rosaceae")$num, label="Rosaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 354) +

  geom_cladelabel(node=subset(famf, Var1=="Apiaceae")$num, label="Apiaceae",  fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Brassicaceae")$num, label="Brassicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 320) +
  geom_cladelabel(node=subset(famf, Var1=="Plantaginaceae")$num, label="Plantaginaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 85) +
  geom_cladelabel(node=subset(famf, Var1=="Caryophyllaceae")$num, label="Caryophyllaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ranunculaceae")$num, label="Ranunculaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 294) +

  geom_cladelabel(node=subset(famf, Var1=="Lamiaceae")$num, label="Lamiaceae",  fontsize=2.5, barsize = 0.1, angle ="auto") +
  geom_cladelabel(node=subset(famf, Var1=="Onagraceae")$num, label="Onagraceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 320) +
  geom_cladelabel(node=subset(famf, Var1=="Polygonaceae")$num, label="Polygonaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 35) +
  geom_cladelabel(node=subset(famf, Var1=="Gentianaceae")$num, label="Gentianaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 70) +
  geom_cladelabel(node=subset(famf, Var1=="Orobanchaceae")$num, label="Orobanchaceae", fontsize=2.5, barsize = 0.1, angle ="auto") +

  geom_cladelabel(node=subset(famf, Var1=="Euphorbiaceae")$num, label="Euphorbiaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 330) +
  geom_cladelabel(node=subset(famf, Var1=="Boraginaceae")$num, label="Boraginaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 83) +
  geom_cladelabel(node=subset(famf, Var1=="Amaranthaceae")$num, label="Amaranthaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 42) +
  geom_cladelabel(node=subset(famf, Var1=="Ericaceae")$num, label="Ericaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Rubiaceae")$num, label="Rubiaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 65) +

  geom_cladelabel(node=subset(famf, Var1=="Juncaceae")$num, label="Juncaceae",  fontsize=2.5, barsize = 0.1,angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Solanaceae")$num, label="Solanaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 64) +
  geom_cladelabel(node=subset(famf, Var1=="Geraniaceae")$num, label="Geraniaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 320) +
  geom_cladelabel(node=subset(famf, Var1=="Amaryllidaceae")$num, label="Amaryllidaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Saxifragaceae")$num, label="Saxifragaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +

  geom_cladelabel(node=subset(famf, Var1=="Polemoniaceae")$num, label="Polemoniaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 51) +
  geom_cladelabel(node=subset(famf, Var1=="Caprifoliaceae")$num, label="Caprifoliaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Malvaceae")$num, label="Malvaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 315) +
  geom_cladelabel(node=subset(famf, Var1=="Primulaceae")$num, label="Primulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Salicaceae")$num, label="Salicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 325) +
  
  geom_cladelabel(node=subset(famf, Var1=="Iridaceae")$num, label="Iridaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Crassulaceae")$num, label="Crassulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Cactaceae")$num, label="Cactaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 40) +
  geom_cladelabel(node=subset(famf, Var1=="Equisetaceae")$num, label="Linaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="bottom")

#save output:
png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_all_multNov22.png",
    res=300,height=8,width=8,units="in"); 
#png("phylo_ring_all_mult_dec.png", res=300,height=8,width=8,units="in"); 
p
dev.off()

#clean-up:
#rm(list = ls())




###
# Filter by treatment = N
###

dat<-subset(species.data, trt_type2=="n")[,c(1,2)] #select "all mult" treatment from original data
dat<-aggregate(dat[, 2], list(dat$species_matched), mean, na.rm=T) #get mean DCi value per species
dat<-dat[dat$Group.1 %in% tree$tip.label, ] #make sure all species in the data are on the tree
rownames(dat)<-dat$Group.1 #set species names as rownames
dat$Group.1<-NULL #and delete column with species names

#prune tree:
tree2<-keep.tip(tree, rownames(dat))


###
# Run function to calculate if each node has significantly higher or lower mean DCi than
# expected if phylogenetic relationships were at random (laod function in the "DCi_nodes_scorre.R" script)
###

res<-node.mean(tree2, dat, 999)
write.table(res, paste(my.wd, "res_phylo_n.csv", sep="")) #save the result
res<-read.table(paste(my.wd, "res_phylo_n.csv", sep=""))
#res2<-subset(res, P_value<0.01) #this would tell you what nodes are significant with alpha < 0.01
#tips(tree2, 1543) #and this would tell you what species are found in that clade

###
# Clean-up result to highlight nodes in the tree
###

significant<-res #create a copy of the main result
significant$P_value[significant$P_value>0.05]<-NA #replace non-significant with NA
significant$P_value[significant$SR<3]<-NA #assign NA to nodes with 2 or less species
significant$P_value[1]<-NA #set the first node (the root node) to NA
significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
significant<-c(rep(NA, length(tree2$tip.label)), significant$P_value) #merge "tip nodes" with "inner" tree nodes
significant<-as.factor(significant) #convert values into factors
significant<-factor(significant, levels = c("neg.05", "pos.05" )) #change order of factors


###
# Assign families to nodes
###

#load family data and clean-up:
fam<-read.table(paste(my.wd, "species_families_trees_2021.csv",sep=""), header=T, sep=",", fill = TRUE) #load data
fam$species_matched<-gsub(" ", "_", fam$species_matched) #adapt species nomenclature
fam$family[fam$family=="Compositae"]<-"Asteraceae" #replace family name
fam$family[fam$family=="Leguminosae"]<-"Fabaceae" #replace family name
fam$family[fam$family=="Viburnaceae"]<-"Adoxaceae" #replace family name
fam<-fam[which(fam$species_matched %in% rownames(dat)),] #subset only species included in our treatment

#get table with ranked families based on their number of species:
famf<-as.data.frame(table(fam$family)) #create dataframe
famf<-famf[order(-famf$Freq),] #order it
row.names(famf) <- NULL #get rid of rownames

#create vector with unique names of species:
list.r <- unique(fam$family)

#get nodes for families:
list.nod<-NULL
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$species_matched) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- findMRCA(tree2, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)


###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=38)) #select the top 38 families with 5 or more species

#Plot tree for N 
#remember to change angle = "auto" everywheree to avoid overlap in names. Consider also unifying "barsize" (to 0.1, for example):
p <- 
  ggtree(tree2, layout="circular", size=0.5)+ # build circular tree
  geom_point(aes(colour=as.factor(significant)), size=2, alpha=1, show.legend = TRUE) + # highlight nodes
  scale_colour_manual(values=c("red", "deepskyblue"), labels=c("Lower DCi", "Higher DCi"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  geom_cladelabel(node=subset(famf, Var1=="Poaceae")$num, label="Poaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asteraceae")$num, label="Asteraceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Fabaceae")$num, label="Fabaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 5) +
  geom_cladelabel(node=subset(famf, Var1=="Cyperaceae")$num, label="Cyperaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Rosaceae")$num, label="Rosaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 354) +
  
  geom_cladelabel(node=subset(famf, Var1=="Apiaceae")$num, label="Apiaceae",  fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Brassicaceae")$num, label="Brassicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 320) +
  geom_cladelabel(node=subset(famf, Var1=="Plantaginaceae")$num, label="Plantaginaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 85) +
  geom_cladelabel(node=subset(famf, Var1=="Caryophyllaceae")$num, label="Caryophyllaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ranunculaceae")$num, label="Ranunculaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 294) +
  
  geom_cladelabel(node=subset(famf, Var1=="Lamiaceae")$num, label="Lamiaceae",  fontsize=2.5, barsize = 0.1, angle ="auto") +
  geom_cladelabel(node=subset(famf, Var1=="Onagraceae")$num, label="Onagraceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 307) +
  geom_cladelabel(node=subset(famf, Var1=="Polygonaceae")$num, label="Polygonaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 35) +
  geom_cladelabel(node=subset(famf, Var1=="Gentianaceae")$num, label="Gentianaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 80) +
  geom_cladelabel(node=subset(famf, Var1=="Orobanchaceae")$num, label="Orobanchaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 90) +
  
  geom_cladelabel(node=subset(famf, Var1=="Euphorbiaceae")$num, label="Euphorbiaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1=="Boraginaceae")$num, label="Boraginaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 83) +
  geom_cladelabel(node=subset(famf, Var1=="Amaranthaceae")$num, label="Amaranthaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 42) +
  geom_cladelabel(node=subset(famf, Var1=="Ericaceae")$num, label="Ericaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 50) +
  geom_cladelabel(node=subset(famf, Var1=="Rubiaceae")$num, label="Rubiaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 65) +
  
  geom_cladelabel(node=subset(famf, Var1=="Juncaceae")$num, label="Juncaceae",  fontsize=2.5, barsize = 0.1,angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Solanaceae")$num, label="Solanaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 64) +
  geom_cladelabel(node=subset(famf, Var1=="Geraniaceae")$num, label="Geraniaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 303) +
  geom_cladelabel(node=subset(famf, Var1=="Apocynaceae")$num, label="Apocynaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 70) +
  geom_cladelabel(node=subset(famf, Var1=="Violaceae")$num, label="Violaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 320) +
  
  geom_cladelabel(node=subset(famf, Var1=="Amaryllidaceae")$num, label="Amaryllidaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asparagaceae")$num, label="Asparagaceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Campanulaceae")$num, label="Campanulaceae", fontsize=2.5, barsize=0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Polemoniaceae")$num, label="Polemoniaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 51) +
  
  geom_cladelabel(node=subset(famf, Var1=="Caprifoliaceae")$num, label="Caprifoliaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Convolvulaceae")$num, label="Convolvulaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 60) +
  geom_cladelabel(node=subset(famf, Var1=="Malvaceae")$num, label="Malvaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 315) +
  geom_cladelabel(node=subset(famf, Var1=="Primulaceae")$num, label="Primulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Salicaceae")$num, label="Salicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 320) +
  
  geom_cladelabel(node=subset(famf, Var1=="Verbenaceae")$num, label="Verbenaceae",  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 90) +
  geom_cladelabel(node=subset(famf, Var1=="Oxalidaceae")$num, label="Oxalidaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 322) +
  geom_cladelabel(node=subset(famf, Var1=="Iridaceae")$num, label="Iridaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Anacardiaceae")$num, label="Anacardiaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 310) +
  geom_cladelabel(node=subset(famf, Var1=="Betulaceae")$num, label="Betulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 340) +
  
  geom_cladelabel(node=subset(famf, Var1=="Nyctaginaceae")$num, label="Nyctaginaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 35) +
  geom_cladelabel(node=subset(famf, Var1=="Orchidaceae")$num, label="Orchidaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Crassulaceae")$num, label="Crassulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Liliaceae")$num, label="Liliaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1=="Cactaceae")$num, label="Cactaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 40) +
  geom_cladelabel(node=subset(famf, Var1=="Cistaceae")$num, label="Cistaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 312) +
  geom_cladelabel(node=subset(famf, Var1=="Commelinaceae")$num, label="Commelinaceae", fontsize=2.5, barsize=0.5, angle = 'auto') +
  geom_cladelabel(node=subset(famf, Var1=="Equisetaceae")$num, label="Equisetaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Lycopodiaceae")$num, label="Lycopodiaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="bottom")

#save output:
png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_n.png",
    res=300,height=8,width=8,units="in"); 
#png("phylo_ring_n.png", res=300,height=8,width=8,units="in"); 
p
dev.off()

###
# Filter by treatment = N +OTHER
###

dat<-subset(species.data, trt_type2=="n_other")[,c(1,2)] #select "all mult" treatment from original data
dat<-aggregate(dat[, 2], list(dat$species_matched), mean, na.rm=T) #get mean DCi value per species
dat<-dat[dat$Group.1 %in% tree$tip.label, ] #make sure all species in the data are on the tree
rownames(dat)<-dat$Group.1 #set species names as rownames
dat$Group.1<-NULL #and delete column with species names

#prune tree:
tree2<-keep.tip(tree, rownames(dat))


###
# Run function to calculate if each node has significantly higher or lower mean DCi than
# expected if phylogenetic relationships were at random (laod function in the "DCi_nodes_scorre.R" script)
###

res<-node.mean(tree2, dat, 999)
#write.table(res, paste(my.wd, "res_phylo_n_other.csv", sep="")) #save the result
#res<-read.table(paste(my.wd, "res_phylo_n_other (3).csv", sep=""))
#res2<-subset(res, P_value<0.01) #this would tell you what nodes are significant with alpha < 0.01
#tips(tree2, 1543) #and this would tell you what species are found in that clade

###
# Clean-up result to highlight nodes in the tree
###

significant<-res #create a copy of the main result
significant$P_value[significant$P_value>0.05]<-NA #replace non-significant with NA
significant$P_value[significant$SR<3]<-NA #assign NA to nodes with 2 or less species
significant$P_value[1]<-NA #set the first node (the root node) to NA
significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
significant<-c(rep(NA, length(tree2$tip.label)), significant$P_value) #merge "tip nodes" with "inner" tree nodes
significant<-as.factor(significant) #convert values into factors
significant<-factor(significant, levels = c("neg.05", "pos.05" )) #change order of factors


###
# Assign families to nodes
###

#load family data and clean-up:
fam<-read.table(paste(my.wd, "species_families_2021.csv",sep=""), header=T, sep=",", fill = TRUE) #load data
fam$species_matched<-gsub(" ", "_", fam$species_matched) #adapt species nomenclature
fam$family[fam$family=="Compositae"]<-"Asteraceae" #replace family name
fam$family[fam$family=="Leguminosae"]<-"Fabaceae" #replace family name
fam<-fam[which(fam$species_matched %in% rownames(dat)),] #subset only species included in our treatment

#get table with ranked families based on their number of species:
famf<-as.data.frame(table(fam$family)) #create dataframe
famf<-famf[order(-famf$Freq),] #order it
row.names(famf) <- NULL #get rid of rownames

#create vector with unique names of species:
list.r <- unique(fam$family)

#get nodes for families:
list.nod<-NULL
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$species_matched) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- findMRCA(tree2, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)


###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=27)) #select the top 27 families with 5 or more species

#Plot tree for N +other
#remember to change angle = "auto" everywhere to avoid overlap in names. Consider also unifying "barsize" (to 0.1, for example):
p <- 
  ggtree(tree2, layout="circular", size=0.5)+ # build circular tree
  geom_point(aes(colour=as.factor(significant)), size=2, alpha=1, show.legend = TRUE) + # highlight nodes
  scale_colour_manual(values=c("orange", "purple"), labels=c("Worse with N + Other", "Better with N + Other"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  geom_cladelabel(node=subset(famf, Var1=="Poaceae")$num, label="Poaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asteraceae")$num, label="Asteraceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Fabaceae")$num, label="Fabaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 5) +
  geom_cladelabel(node=subset(famf, Var1=="Cyperaceae")$num, label="Cyperaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Rosaceae")$num, label="Rosaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 354) +
  
  geom_cladelabel(node=subset(famf, Var1=="Apiaceae")$num, label="Apiaceae",  fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Brassicaceae")$num, label="Brassicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 320) +
  geom_cladelabel(node=subset(famf, Var1=="Plantaginaceae")$num, label="Plantaginaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 85) +
  geom_cladelabel(node=subset(famf, Var1=="Caryophyllaceae")$num, label="Caryophyllaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ranunculaceae")$num, label="Ranunculaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 294) +
  
  geom_cladelabel(node=subset(famf, Var1=="Lamiaceae")$num, label="Lamiaceae",  fontsize=2.5, barsize = 0.1, angle ="auto") +
  geom_cladelabel(node=subset(famf, Var1=="Onagraceae")$num, label="Onagraceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 307) +
  geom_cladelabel(node=subset(famf, Var1=="Polygonaceae")$num, label="Polygonaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 35) +
  geom_cladelabel(node=subset(famf, Var1=="Gentianaceae")$num, label="Gentianaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 80) +
  geom_cladelabel(node=subset(famf, Var1=="Orobanchaceae")$num, label="Orobanchaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 90) +
  
  geom_cladelabel(node=subset(famf, Var1=="Euphorbiaceae")$num, label="Euphorbiaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1=="Boraginaceae")$num, label="Boraginaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 83) +
  geom_cladelabel(node=subset(famf, Var1=="Amaranthaceae")$num, label="Amaranthaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 42) +
  geom_cladelabel(node=subset(famf, Var1=="Ericaceae")$num, label="Ericaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 50) +
  geom_cladelabel(node=subset(famf, Var1=="Rubiaceae")$num, label="Rubiaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 65) +
  
  geom_cladelabel(node=subset(famf, Var1=="Juncaceae")$num, label="Juncaceae",  fontsize=2.5, barsize = 0.1,angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Solanaceae")$num, label="Solanaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 64) +
  geom_cladelabel(node=subset(famf, Var1=="Geraniaceae")$num, label="Geraniaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 303) +
  geom_cladelabel(node=subset(famf, Var1=="Apocynaceae")$num, label="Apocynaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 70) +
  geom_cladelabel(node=subset(famf, Var1=="Violaceae")$num, label="Violaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 320) +
  
  geom_cladelabel(node=subset(famf, Var1=="Amaryllidaceae")$num, label="Amaryllidaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asparagaceae")$num, label="Asparagaceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Campanulaceae")$num, label="Campanulaceae", fontsize=2.5, barsize=0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Polemoniaceae")$num, label="Polemoniaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 51) +
  
  geom_cladelabel(node=subset(famf, Var1=="Caprifoliaceae")$num, label="Caprifoliaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Convolvulaceae")$num, label="Convolvulaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 60) +
  geom_cladelabel(node=subset(famf, Var1=="Malvaceae")$num, label="Malvaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 315) +
  geom_cladelabel(node=subset(famf, Var1=="Primulaceae")$num, label="Primulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Salicaceae")$num, label="Salicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 320) +
  
  geom_cladelabel(node=subset(famf, Var1=="Verbenaceae")$num, label="Verbenaceae",  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 90) +
  geom_cladelabel(node=subset(famf, Var1=="Oxalidaceae")$num, label="Oxalidaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 322) +
  geom_cladelabel(node=subset(famf, Var1=="Iridaceae")$num, label="Iridaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Anacardiaceae")$num, label="Anacardiaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 310) +
  geom_cladelabel(node=subset(famf, Var1=="Betulaceae")$num, label="Betulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 340) +
  
  geom_cladelabel(node=subset(famf, Var1=="Nyctaginaceae")$num, label="Nyctaginaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 35) +
  geom_cladelabel(node=subset(famf, Var1=="Orchidaceae")$num, label="Orchidaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Crassulaceae")$num, label="Crassulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Liliaceae")$num, label="Liliaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1=="Cactaceae")$num, label="Cactaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 40) +
  geom_cladelabel(node=subset(famf, Var1=="Cistaceae")$num, label="Cistaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 312) +
  geom_cladelabel(node=subset(famf, Var1=="Commelinaceae")$num, label="Commelinaceae", fontsize=2.5, barsize=0.5, angle = 'auto') +
  geom_cladelabel(node=subset(famf, Var1=="Equisetaceae")$num, label="Equisetaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Lycopodiaceae")$num, label="Lycopodiaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="bottom")

#save output:
png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_n_other.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()




###
# Filter by treatment = P
###

dat<-subset(species.data, trt_type2=="p")[,c(1,2)] #select "all mult" treatment from original data
dat<-aggregate(dat[, 2], list(dat$species_matched), mean, na.rm=T) #get mean DCi value per species
dat<-dat[dat$Group.1 %in% tree$tip.label, ] #make sure all species in the data are on the tree
rownames(dat)<-dat$Group.1 #set species names as rownames
dat$Group.1<-NULL #and delete column with species names

#prune tree:
tree2<-keep.tip(tree, rownames(dat))


###
# Run function to calculate if each node has significantly higher or lower mean DCi than
# expected if phylogenetic relationships were at random (laod function in the "DCi_nodes_scorre.R" script)
###

res<-node.mean(tree2, dat, 999)
write.table(res, paste(my.wd, "res_phylo_p.csv", sep="")) #save the result
res<-read.table(paste(my.wd, "res_phylo_p.csv", sep=""))
#res2<-subset(res, P_value<0.01) #this would tell you what nodes are significant with alpha < 0.01
#tips(tree2, 1543) #and this would tell you what species are found in that clade

###
# Clean-up result to highlight nodes in the tree
###

significant<-res #create a copy of the main result
significant$P_value[significant$P_value>0.05]<-NA #replace non-significant with NA
significant$P_value[significant$SR<3]<-NA #assign NA to nodes with 2 or less species
significant$P_value[1]<-NA #set the first node (the root node) to NA
significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
significant<-c(rep(NA, length(tree2$tip.label)), significant$P_value) #merge "tip nodes" with "inner" tree nodes
significant<-as.factor(significant) #convert values into factors
significant<-factor(significant, levels = c("neg.05", "pos.05" )) #change order of factors


###
# Assign families to nodes
###

#load family data and clean-up:
fam<-read.table(paste(my.wd, "species_families_2021.csv",sep=""), header=T, sep=",", fill = TRUE) #load data
fam$species_matched<-gsub(" ", "_", fam$species_matched) #adapt species nomenclature
fam$family[fam$family=="Compositae"]<-"Asteraceae" #replace family name
fam$family[fam$family=="Leguminosae"]<-"Fabaceae" #replace family name
fam<-fam[which(fam$species_matched %in% rownames(dat)),] #subset only species included in our treatment

#get table with ranked families based on their number of species:
famf<-as.data.frame(table(fam$family)) #create dataframe
famf<-famf[order(-famf$Freq),] #order it
row.names(famf) <- NULL #get rid of rownames

#create vector with unique names of species:
list.r <- unique(fam$family)

#get nodes for families:
list.nod<-NULL
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$species_matched) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- findMRCA(tree2, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)


###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=31)) #select the top 31 families with 5 or more species

#Plot tree for P
#remember to change angle = "auto" everywheree to avoid overlap in names. Consider also unifying "barsize" (to 0.1, for example):
p <- 
  ggtree(tree2, layout="circular", size=0.5)+ # build circular tree
  geom_point(aes(colour=as.factor(significant)), size=2, alpha=1, show.legend = TRUE) + # highlight nodes
  scale_colour_manual(values=c("red", "deepskyblue"), labels=c("Lower DCi", "Higher DCi"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  geom_cladelabel(node=subset(famf, Var1=="Poaceae")$num, label="Poaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asteraceae")$num, label="Asteraceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Fabaceae")$num, label="Fabaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 5) +
  geom_cladelabel(node=subset(famf, Var1=="Cyperaceae")$num, label="Cyperaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Rosaceae")$num, label="Rosaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 354) +
  
  geom_cladelabel(node=subset(famf, Var1=="Apiaceae")$num, label="Apiaceae",  fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Brassicaceae")$num, label="Brassicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 320) +
  geom_cladelabel(node=subset(famf, Var1=="Plantaginaceae")$num, label="Plantaginaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 85) +
  geom_cladelabel(node=subset(famf, Var1=="Caryophyllaceae")$num, label="Caryophyllaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ranunculaceae")$num, label="Ranunculaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 294) +
  
  geom_cladelabel(node=subset(famf, Var1=="Lamiaceae")$num, label="Lamiaceae",  fontsize=2.5, barsize = 0.1, angle ="auto") +
  geom_cladelabel(node=subset(famf, Var1=="Onagraceae")$num, label="Onagraceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 307) +
  geom_cladelabel(node=subset(famf, Var1=="Polygonaceae")$num, label="Polygonaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 35) +
  geom_cladelabel(node=subset(famf, Var1=="Gentianaceae")$num, label="Gentianaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 80) +
  geom_cladelabel(node=subset(famf, Var1=="Orobanchaceae")$num, label="Orobanchaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 90) +
  
  geom_cladelabel(node=subset(famf, Var1=="Euphorbiaceae")$num, label="Euphorbiaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1=="Boraginaceae")$num, label="Boraginaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 83) +
  geom_cladelabel(node=subset(famf, Var1=="Amaranthaceae")$num, label="Amaranthaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 42) +
  geom_cladelabel(node=subset(famf, Var1=="Ericaceae")$num, label="Ericaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 50) +
  geom_cladelabel(node=subset(famf, Var1=="Rubiaceae")$num, label="Rubiaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 65) +
  
  geom_cladelabel(node=subset(famf, Var1=="Juncaceae")$num, label="Juncaceae",  fontsize=2.5, barsize = 0.1,angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Solanaceae")$num, label="Solanaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 64) +
  geom_cladelabel(node=subset(famf, Var1=="Geraniaceae")$num, label="Geraniaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 303) +
  geom_cladelabel(node=subset(famf, Var1=="Apocynaceae")$num, label="Apocynaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 70) +
  geom_cladelabel(node=subset(famf, Var1=="Violaceae")$num, label="Violaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 338) +
  
  geom_cladelabel(node=subset(famf, Var1=="Amaryllidaceae")$num, label="Amaryllidaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asparagaceae")$num, label="Asparagaceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Campanulaceae")$num, label="Campanulaceae", fontsize=2.5, barsize=0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Saxifragaceae")$num, label="Saxifragaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Polemoniaceae")$num, label="Polemoniaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 51) +
  
  geom_cladelabel(node=subset(famf, Var1=="Caprifoliaceae")$num, label="Caprifoliaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Convolvulaceae")$num, label="Convolvulaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 58) +
  geom_cladelabel(node=subset(famf, Var1=="Malvaceae")$num, label="Malvaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 315) +
  geom_cladelabel(node=subset(famf, Var1=="Primulaceae")$num, label="Primulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 50) +
  geom_cladelabel(node=subset(famf, Var1=="Salicaceae")$num, label="Salicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 335) +
  
  geom_cladelabel(node=subset(famf, Var1=="Verbenaceae")$num, label="Verbenaceae",  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 90) +
  geom_cladelabel(node=subset(famf, Var1=="Oxalidaceae")$num, label="Oxalidaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 322) +
  geom_cladelabel(node=subset(famf, Var1=="Iridaceae")$num, label="Iridaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Anacardiaceae")$num, label="Anacardiaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 310) +
  geom_cladelabel(node=subset(famf, Var1=="Betulaceae")$num, label="Betulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 347) +
  
  geom_cladelabel(node=subset(famf, Var1=="Fagaceae")$num, label="Fagaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 340) +
  geom_cladelabel(node=subset(famf, Var1=="Nyctaginaceae")$num, label="Nyctaginaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 39) +
  geom_cladelabel(node=subset(famf, Var1=="Orchidaceae")$num, label="Orchidaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Crassulaceae")$num, label="Crassulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Liliaceae")$num, label="Liliaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1=="Cactaceae")$num, label="Cactaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 40) +
  geom_cladelabel(node=subset(famf, Var1=="Cistaceae")$num, label="Cistaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 312) +
  geom_cladelabel(node=subset(famf, Var1=="Commelinaceae")$num, label="Commelinaceae", fontsize=2.5, barsize=0.5, angle = 'auto') +
  geom_cladelabel(node=subset(famf, Var1=="Equisetaceae")$num, label="Equisetaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Lycopodiaceae")$num, label="Lycopodiaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="bottom")
#save output:
png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_p.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()

###
# Filter by treatment == Co2
###

dat<-subset(species.data, trt_type2=="CO2")[,c(1,2)] #select "all mult" treatment from original data
dat<-aggregate(dat[, 2], list(dat$species_matched), mean, na.rm=T) #get mean DCi value per species - this is not necessary - is already the average.
dat<-dat[dat$Group.1 %in% tree$tip.label, ] #make sure all species in the data are on the tree
rownames(dat)<-dat$Group.1 #set species names as rownames
dat$Group.1<-NULL #and delete column with species names

#prune tree:
tree2<-keep.tip(tree, rownames(dat))


###
# Run function to calculate if each node has significantly higher or lower mean DCi than
# expected if phylogenetic relationships were at random (laod function in the "DCi_nodes_scorre.R" script)
###

res<-node.mean(tree2, dat, 999)
write.table(res, paste(my.wd, "res_phylo_co2.csv", sep="")) #save the result
res<-read.table(paste(my.wd, "res_phylo_co2.csv", sep=""))
#res2<-subset(res, P_value<0.01) #this would tell you what nodes are significant with alpha < 0.01
#tips(tree2, 1543) #and this would tell you what species are found in that clade

###
# Clean-up result to highlight nodes in the tree
###

significant<-res #create a copy of the main result
significant$P_value[significant$P_value>0.05]<-NA #replace non-significant with NA
significant$P_value[significant$SR<3]<-NA #assign NA to nodes with 2 or less species

significant$P_value[1]<-NA #set the first node (the root node) to NA
significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
significant<-c(rep(NA, length(tree2$tip.label)), significant$P_value) #merge "tip nodes" with "inner" tree nodes
significant<-as.factor(significant) #convert values into factors

significant<-factor(significant, levels = c("neg.05", "pos.05" )) #change order of factors



###
# Assign families to nodes
###

#load family data and clean-up:
fam<-read.table(paste(my.wd, "species_families_2021.csv",sep=""), header=T, sep=",", fill = TRUE) #load data
fam$species_matched<-gsub(" ", "_", fam$species_matched) #adapt species nomenclature
fam$family[fam$family=="Compositae"]<-"Asteraceae" #replace family name
fam$family[fam$family=="Leguminosae"]<-"Fabaceae" #replace family name
fam<-fam[which(fam$species_matched %in% rownames(dat)),] #subset only species included in our treatment

#get table with ranked families based on their number of species:
famf<-as.data.frame(table(fam$family)) #create dataframe
famf<-famf[order(-famf$Freq),] #order it
row.names(famf) <- NULL #get rid of rownames

#create vector with unique names of species:
list.r <- unique(fam$family)

#get nodes for families:
list.nod<-NULL
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$species_matched) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- findMRCA(tree2, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)


###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=9)) #select the top 9 families with 5 for more species

#Plot tree - for CO2
#remember to change angle = "auto" everywheree to avoid overlap in names. Consider also unifying "barsize" (to 0.1, for example):
p <- 
  ggtree(tree2, layout="circular", size=0.5)+ # build circular tree
  geom_point(aes(colour=as.factor(significant)), size=2, alpha=1, show.legend = TRUE) + # highlight nodes
  scale_colour_manual(values=c("red", "deepskyblue"), labels=c("Lower DCi", "Higher DCi"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  geom_cladelabel(node=subset(famf, Var1=="Poaceae")$num, label="Poaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asteraceae")$num, label="Asteraceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Fabaceae")$num, label="Fabaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 10) +
  geom_cladelabel(node=subset(famf, Var1=="Cyperaceae")$num, label="Cyperaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Rosaceae")$num, label="Rosaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 15) +
  
  geom_cladelabel(node=subset(famf, Var1=="Apiaceae")$num, label="Apiaceae",  fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Brassicaceae")$num, label="Brassicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1=="Plantaginaceae")$num, label="Plantaginaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 85) +
  geom_cladelabel(node=subset(famf, Var1=="Caryophyllaceae")$num, label="Caryophyllaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ranunculaceae")$num, label="Ranunculaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 294) +
  
  geom_cladelabel(node=subset(famf, Var1=="Lamiaceae")$num, label="Lamiaceae",  fontsize=2.5, barsize = 0.1, angle ="auto") +
  geom_cladelabel(node=subset(famf, Var1=="Onagraceae")$num, label="Onagraceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 307) +
  geom_cladelabel(node=subset(famf, Var1=="Polygonaceae")$num, label="Polygonaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 40) +
  geom_cladelabel(node=subset(famf, Var1=="Gentianaceae")$num, label="Gentianaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 80) +
  geom_cladelabel(node=subset(famf, Var1=="Orobanchaceae")$num, label="Orobanchaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 90) +
  
  geom_cladelabel(node=subset(famf, Var1=="Euphorbiaceae")$num, label="Euphorbiaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1=="Boraginaceae")$num, label="Boraginaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 83) +
  geom_cladelabel(node=subset(famf, Var1=="Amaranthaceae")$num, label="Amaranthaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ericaceae")$num, label="Ericaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Rubiaceae")$num, label="Rubiaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 65) +
  
  geom_cladelabel(node=subset(famf, Var1=="Juncaceae")$num, label="Juncaceae",  fontsize=2.5, barsize = 0.1,angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Solanaceae")$num, label="Solanaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 64) +
  geom_cladelabel(node=subset(famf, Var1=="Geraniaceae")$num, label="Geraniaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 303) +
  geom_cladelabel(node=subset(famf, Var1=="Apocynaceae")$num, label="Apocynaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 70) +
  geom_cladelabel(node=subset(famf, Var1=="Violaceae")$num, label="Violaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 340) +
  
  geom_cladelabel(node=subset(famf, Var1=="Amaryllidaceae")$num, label="Amaryllidaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asparagaceae")$num, label="Asparagaceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Campanulaceae")$num, label="Campanulaceae", fontsize=2.5, barsize=0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Saxifragaceae")$num, label="Saxifragaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Polemoniaceae")$num, label="Polemoniaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 51) +
  
  geom_cladelabel(node=subset(famf, Var1=="Caprifoliaceae")$num, label="Caprifoliaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Convolvulaceae")$num, label="Convolvulaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 70) +
  geom_cladelabel(node=subset(famf, Var1=="Malvaceae")$num, label="Malvaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 320) +
  geom_cladelabel(node=subset(famf, Var1=="Primulaceae")$num, label="Primulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Salicaceae")$num, label="Salicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 340) +
  
  geom_cladelabel(node=subset(famf, Var1=="Verbenaceae")$num, label="Verbenaceae",  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 90) +
  geom_cladelabel(node=subset(famf, Var1=="Oxalidaceae")$num, label="Oxalidaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 335) +
  geom_cladelabel(node=subset(famf, Var1=="Iridaceae")$num, label="Iridaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Anacardiaceae")$num, label="Anacardiaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 315) +
  geom_cladelabel(node=subset(famf, Var1=="Betulaceae")$num, label="Betulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 347) +
  
  geom_cladelabel(node=subset(famf, Var1=="Fagaceae")$num, label="Fagaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 340) +
  geom_cladelabel(node=subset(famf, Var1=="Nyctaginaceae")$num, label="Nyctaginaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 39) +
  geom_cladelabel(node=subset(famf, Var1=="Orchidaceae")$num, label="Orchidaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Crassulaceae")$num, label="Crassulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Liliaceae")$num, label="Liliaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1=="Cactaceae")$num, label="Cactaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 40) +
  geom_cladelabel(node=subset(famf, Var1=="Cistaceae")$num, label="Cistaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 312) +
  geom_cladelabel(node=subset(famf, Var1=="Commelinaceae")$num, label="Commelinaceae", fontsize=2.5, barsize=0.5, angle = 'auto') +
  geom_cladelabel(node=subset(famf, Var1=="Equisetaceae")$num, label="Equisetaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Lycopodiaceae")$num, label="Lycopodiaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="bottom")

#save output:
png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_co2.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()


###
# Filter by treatment == drought
###

dat<-subset(species.data, trt_type2=="drt")[,c(1,2)] #select "all mult" treatment from original data
dat<-aggregate(dat[, 2], list(dat$species_matched), mean, na.rm=T) #get mean DCi value per species - this is not necessary - is already the average.
dat<-dat[dat$Group.1 %in% tree$tip.label, ] #make sure all species in the data are on the tree
rownames(dat)<-dat$Group.1 #set species names as rownames
dat$Group.1<-NULL #and delete column with species names

#prune tree:
tree2<-keep.tip(tree, rownames(dat))


###
# Run function to calculate if each node has significantly higher or lower mean DCi than
# expected if phylogenetic relationships were at random (laod function in the "DCi_nodes_scorre.R" script)
###

res<-node.mean(tree2, dat, 999)
write.table(res, paste(my.wd, "res_phylo_drought.csv", sep="")) #save the result
res<-read.table(paste(my.wd, "res_phylo_drought.csv", sep=""))
#res2<-subset(res, P_value<0.01) #this would tell you what nodes are significant with alpha < 0.01
#tips(tree2, 1543) #and this would tell you what species are found in that clade

###
# Clean-up result to highlight nodes in the tree
###

significant<-res #create a copy of the main result
significant$P_value[significant$P_value>0.05]<-NA #replace non-significant with NA
significant$P_value[significant$SR<3]<-NA #assign NA to nodes with 2 or less species

significant$P_value[1]<-NA #set the first node (the root node) to NA
significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
significant<-c(rep(NA, length(tree2$tip.label)), significant$P_value) #merge "tip nodes" with "inner" tree nodes
significant<-as.factor(significant) #convert values into factors

significant<-factor(significant, levels = c("neg.05", "pos.05" )) #change order of factors



###
# Assign families to nodes
###

#load family data and clean-up:
fam<-read.table(paste(my.wd, "species_families_2021.csv",sep=""), header=T, sep=",", fill = TRUE) #load data
fam$species_matched<-gsub(" ", "_", fam$species_matched) #adapt species nomenclature
fam$family[fam$family=="Compositae"]<-"Asteraceae" #replace family name
fam$family[fam$family=="Leguminosae"]<-"Fabaceae" #replace family name
fam<-fam[which(fam$species_matched %in% rownames(dat)),] #subset only species included in our treatment

#get table with ranked families based on their number of species:
famf<-as.data.frame(table(fam$family)) #create dataframe
famf<-famf[order(-famf$Freq),] #order it
row.names(famf) <- NULL #get rid of rownames

#create vector with unique names of species:
list.r <- unique(fam$family)

#get nodes for families:
list.nod<-NULL
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$species_matched) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- findMRCA(tree2, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)


###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=23)) #select the top 23 families with 5 for more species

#Plot tree 
#remember to change angle = "auto" everywheree to avoid overlap in names. Consider also unifying "barsize" (to 0.1, for example):
p <- 
  ggtree(tree2, layout="circular", size=0.5)+ # build circular tree
  geom_point(aes(colour=as.factor(significant)), size=2, alpha=1, show.legend = TRUE) + # highlight nodes
  scale_colour_manual(values=c("red", "deepskyblue"), labels=c("Lower DCi", "Higher DCi"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  geom_cladelabel(node=subset(famf, Var1=="Poaceae")$num, label="Poaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asteraceae")$num, label="Asteraceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Fabaceae")$num, label="Fabaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 5) +
  geom_cladelabel(node=subset(famf, Var1=="Cyperaceae")$num, label="Cyperaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Rosaceae")$num, label="Rosaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 330) +
  
  geom_cladelabel(node=subset(famf, Var1=="Apiaceae")$num, label="Apiaceae",  fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Brassicaceae")$num, label="Brassicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Plantaginaceae")$num, label="Plantaginaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 85) +
  geom_cladelabel(node=subset(famf, Var1=="Caryophyllaceae")$num, label="Caryophyllaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ranunculaceae")$num, label="Ranunculaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1=="Lamiaceae")$num, label="Lamiaceae",  fontsize=2.5, barsize = 0.1, angle ="auto") +
  geom_cladelabel(node=subset(famf, Var1=="Onagraceae")$num, label="Onagraceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 275) +
  geom_cladelabel(node=subset(famf, Var1=="Polygonaceae")$num, label="Polygonaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 15) +
  geom_cladelabel(node=subset(famf, Var1=="Gentianaceae")$num, label="Gentianaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 75) +
  geom_cladelabel(node=subset(famf, Var1=="Orobanchaceae")$num, label="Orobanchaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 90) +
  
  geom_cladelabel(node=subset(famf, Var1=="Euphorbiaceae")$num, label="Euphorbiaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1=="Boraginaceae")$num, label="Boraginaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 83) +
  geom_cladelabel(node=subset(famf, Var1=="Amaranthaceae")$num, label="Amaranthaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ericaceae")$num, label="Ericaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Rubiaceae")$num, label="Rubiaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 65) +
  
  geom_cladelabel(node=subset(famf, Var1=="Juncaceae")$num, label="Juncaceae",  fontsize=2.5, barsize = 0.1,angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Solanaceae")$num, label="Solanaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 64) +
  geom_cladelabel(node=subset(famf, Var1=="Geraniaceae")$num, label="Geraniaceae", fontsize=2.5, barsize=0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Apocynaceae")$num, label="Apocynaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 70) +
  geom_cladelabel(node=subset(famf, Var1=="Violaceae")$num, label="Violaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 310) +
  
  geom_cladelabel(node=subset(famf, Var1=="Amaryllidaceae")$num, label="Amaryllidaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asparagaceae")$num, label="Asparagaceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Campanulaceae")$num, label="Campanulaceae", fontsize=2.5, barsize=0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Saxifragaceae")$num, label="Saxifragaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Polemoniaceae")$num, label="Polemoniaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 55) +
  
  geom_cladelabel(node=subset(famf, Var1=="Caprifoliaceae")$num, label="Caprifoliaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Convolvulaceae")$num, label="Convolvulaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 65) +
  geom_cladelabel(node=subset(famf, Var1=="Malvaceae")$num, label="Malvaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 290) +
  geom_cladelabel(node=subset(famf, Var1=="Primulaceae")$num, label="Primulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Salicaceae")$num, label="Salicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 340) +
  
  geom_cladelabel(node=subset(famf, Var1=="Verbenaceae")$num, label="Verbenaceae",  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 90) +
  geom_cladelabel(node=subset(famf, Var1=="Oxalidaceae")$num, label="Oxalidaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 305) +
  geom_cladelabel(node=subset(famf, Var1=="Iridaceae")$num, label="Iridaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Anacardiaceae")$num, label="Anacardiaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 315) +
  geom_cladelabel(node=subset(famf, Var1=="Betulaceae")$num, label="Betulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 347) +
  
  geom_cladelabel(node=subset(famf, Var1=="Fagaceae")$num, label="Fagaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 340) +
  geom_cladelabel(node=subset(famf, Var1=="Nyctaginaceae")$num, label="Nyctaginaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 20) +
  geom_cladelabel(node=subset(famf, Var1=="Orchidaceae")$num, label="Orchidaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Crassulaceae")$num, label="Crassulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Liliaceae")$num, label="Liliaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1=="Cactaceae")$num, label="Cactaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 40) +
  geom_cladelabel(node=subset(famf, Var1=="Cistaceae")$num, label="Cistaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 280) +
  geom_cladelabel(node=subset(famf, Var1=="Commelinaceae")$num, label="Commelinaceae", fontsize=2.5, barsize=0.5, angle = 'auto') +
  geom_cladelabel(node=subset(famf, Var1=="Equisetaceae")$num, label="Equisetaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Lycopodiaceae")$num, label="Lycopodiaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="bottom")

#save output:
png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_drought.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()

###
# Filter by treatment == irrigation
###

dat<-subset(species.data, trt_type2=="irg")[,c(1,2)] #select "all mult" treatment from original data
dat<-aggregate(dat[, 2], list(dat$species_matched), mean, na.rm=T) #get mean DCi value per species - this is not necessary - is already the average.
dat<-dat[dat$Group.1 %in% tree$tip.label, ] #make sure all species in the data are on the tree
rownames(dat)<-dat$Group.1 #set species names as rownames
dat$Group.1<-NULL #and delete column with species names

#prune tree:
tree2<-keep.tip(tree, rownames(dat))


###
# Run function to calculate if each node has significantly higher or lower mean DCi than
# expected if phylogenetic relationships were at random (laod function in the "DCi_nodes_scorre.R" script)
###

res<-node.mean(tree2, dat, 999)
write.table(res, paste(my.wd, "res_phylo_irrigation.csv", sep=""))
res<-read.table(paste(my.wd, "res_phylo_irrigation.csv", sep=""))#save the result
#res2<-subset(res, P_value<0.01) #this would tell you what nodes are significant with alpha < 0.01
#tips(tree2, 1543) #and this would tell you what species are found in that clade

###
# Clean-up result to highlight nodes in the tree
###

significant<-res #create a copy of the main result
significant$P_value[significant$P_value>0.05]<-NA #replace non-significant with NA
significant$P_value[significant$SR<3]<-NA #assign NA to nodes with 2 or less species

significant$P_value[1]<-NA #set the first node (the root node) to NA
significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
significant<-c(rep(NA, length(tree2$tip.label)), significant$P_value) #merge "tip nodes" with "inner" tree nodes
significant<-as.factor(significant) #convert values into factors

significant<-factor(significant, levels = c("neg.05", "pos.05" )) #change order of factors



###
# Assign families to nodes
###

#load family data and clean-up:
fam<-read.table(paste(my.wd, "species_families_2021.csv",sep=""), header=T, sep=",", fill = TRUE) #load data
fam$species_matched<-gsub(" ", "_", fam$species_matched) #adapt species nomenclature
fam$family[fam$family=="Compositae"]<-"Asteraceae" #replace family name
fam$family[fam$family=="Leguminosae"]<-"Fabaceae" #replace family name
fam<-fam[which(fam$species_matched %in% rownames(dat)),] #subset only species included in our treatment

#get table with ranked families based on their number of species:
famf<-as.data.frame(table(fam$family)) #create dataframe
famf<-famf[order(-famf$Freq),] #order it
row.names(famf) <- NULL #get rid of rownames

#create vector with unique names of species:
list.r <- unique(fam$family)

#get nodes for families:
list.nod<-NULL
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$species_matched) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- findMRCA(tree2, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)


###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=24)) #select the top 24 families with 5 for more species

#Plot tree 
#remember to change angle = "auto" everywheree to avoid overlap in names. Consider also unifying "barsize" (to 0.1, for example):
p <- 
  ggtree(tree2, layout="circular", size=0.5)+ # build circular tree
  geom_point(aes(colour=as.factor(significant)), size=2, alpha=1, show.legend = TRUE) + # highlight nodes
  scale_colour_manual(values=c("red", "deepskyblue"), labels=c("Lower DCi", "Higher DCi"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  geom_cladelabel(node=subset(famf, Var1=="Poaceae")$num, label="Poaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asteraceae")$num, label="Asteraceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Fabaceae")$num, label="Fabaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 5) +
  geom_cladelabel(node=subset(famf, Var1=="Cyperaceae")$num, label="Cyperaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Rosaceae")$num, label="Rosaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 330) +
  
  geom_cladelabel(node=subset(famf, Var1=="Apiaceae")$num, label="Apiaceae",  fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Brassicaceae")$num, label="Brassicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Plantaginaceae")$num, label="Plantaginaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 85) +
  geom_cladelabel(node=subset(famf, Var1=="Caryophyllaceae")$num, label="Caryophyllaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ranunculaceae")$num, label="Ranunculaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 280) +
  
  geom_cladelabel(node=subset(famf, Var1=="Lamiaceae")$num, label="Lamiaceae",  fontsize=2.5, barsize = 0.1, angle ="auto") +
  geom_cladelabel(node=subset(famf, Var1=="Onagraceae")$num, label="Onagraceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Polygonaceae")$num, label="Polygonaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 15) +
  geom_cladelabel(node=subset(famf, Var1=="Gentianaceae")$num, label="Gentianaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 75) +
  geom_cladelabel(node=subset(famf, Var1=="Orobanchaceae")$num, label="Orobanchaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 90) +
  
  geom_cladelabel(node=subset(famf, Var1=="Euphorbiaceae")$num, label="Euphorbiaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1=="Boraginaceae")$num, label="Boraginaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 83) +
  geom_cladelabel(node=subset(famf, Var1=="Amaranthaceae")$num, label="Amaranthaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ericaceae")$num, label="Ericaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Rubiaceae")$num, label="Rubiaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 65) +
  
  geom_cladelabel(node=subset(famf, Var1=="Juncaceae")$num, label="Juncaceae",  fontsize=2.5, barsize = 0.1,angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Solanaceae")$num, label="Solanaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 64) +
  geom_cladelabel(node=subset(famf, Var1=="Geraniaceae")$num, label="Geraniaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 290) +
  geom_cladelabel(node=subset(famf, Var1=="Apocynaceae")$num, label="Apocynaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 70) +
  geom_cladelabel(node=subset(famf, Var1=="Violaceae")$num, label="Violaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 310) +
  
  geom_cladelabel(node=subset(famf, Var1=="Amaryllidaceae")$num, label="Amaryllidaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asparagaceae")$num, label="Asparagaceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Campanulaceae")$num, label="Campanulaceae", fontsize=2.5, barsize=0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Saxifragaceae")$num, label="Saxifragaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Polemoniaceae")$num, label="Polemoniaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 55) +
  
  geom_cladelabel(node=subset(famf, Var1=="Caprifoliaceae")$num, label="Caprifoliaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Convolvulaceae")$num, label="Convolvulaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 65) +
  geom_cladelabel(node=subset(famf, Var1=="Malvaceae")$num, label="Malvaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 290) +
  geom_cladelabel(node=subset(famf, Var1=="Primulaceae")$num, label="Primulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Salicaceae")$num, label="Salicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 340) +
  
  geom_cladelabel(node=subset(famf, Var1=="Verbenaceae")$num, label="Verbenaceae",  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 90) +
  geom_cladelabel(node=subset(famf, Var1=="Oxalidaceae")$num, label="Oxalidaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 305) +
  geom_cladelabel(node=subset(famf, Var1=="Iridaceae")$num, label="Iridaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Anacardiaceae")$num, label="Anacardiaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 315) +
  geom_cladelabel(node=subset(famf, Var1=="Betulaceae")$num, label="Betulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 347) +
  
  geom_cladelabel(node=subset(famf, Var1=="Fagaceae")$num, label="Fagaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 340) +
  geom_cladelabel(node=subset(famf, Var1=="Nyctaginaceae")$num, label="Nyctaginaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 20) +
  geom_cladelabel(node=subset(famf, Var1=="Orchidaceae")$num, label="Orchidaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Crassulaceae")$num, label="Crassulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Liliaceae")$num, label="Liliaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1=="Cactaceae")$num, label="Cactaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 40) +
  geom_cladelabel(node=subset(famf, Var1=="Cistaceae")$num, label="Cistaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 280) +
  geom_cladelabel(node=subset(famf, Var1=="Commelinaceae")$num, label="Commelinaceae", fontsize=2.5, barsize=0.5, angle = 'auto') +
  geom_cladelabel(node=subset(famf, Var1=="Equisetaceae")$num, label="Equisetaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Lycopodiaceae")$num, label="Lycopodiaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="bottom")

#save output:
png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_irrigation.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()

###
# Filter by treatment == temp
###

dat<-subset(species.data, trt_type2=="temp")[,c(1,2)] #select "all mult" treatment from original data
dat<-aggregate(dat[, 2], list(dat$species_matched), mean, na.rm=T) #get mean DCi value per species - this is not necessary - is already the average.
dat<-dat[dat$Group.1 %in% tree$tip.label, ] #make sure all species in the data are on the tree
rownames(dat)<-dat$Group.1 #set species names as rownames
dat$Group.1<-NULL #and delete column with species names

#prune tree:
tree2<-keep.tip(tree, rownames(dat))


###
# Run function to calculate if each node has significantly higher or lower mean DCi than
# expected if phylogenetic relationships were at random (laod function in the "DCi_nodes_scorre.R" script)
###

res<-node.mean(tree2, dat, 999)
write.table(res, paste(my.wd, "res_phylo_temp.csv", sep=""))
res<-read.table(paste(my.wd, "res_phylo_temp.csv", sep=""))#save the result
#res2<-subset(res, P_value<0.01) #this would tell you what nodes are significant with alpha < 0.01
#tips(tree2, 1543) #and this would tell you what species are found in that clade

###
# Clean-up result to highlight nodes in the tree
###

significant<-res #create a copy of the main result
significant$P_value[significant$P_value>0.05]<-NA #replace non-significant with NA
significant$P_value[significant$SR<3]<-NA #assign NA to nodes with 2 or less species

significant$P_value[1]<-NA #set the first node (the root node) to NA
significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
significant<-c(rep(NA, length(tree2$tip.label)), significant$P_value) #merge "tip nodes" with "inner" tree nodes
significant<-as.factor(significant) #convert values into factors

significant<-factor(significant, levels = c("neg.05", "pos.05" )) #change order of factors



###
# Assign families to nodes
###

#load family data and clean-up:
fam<-read.table(paste(my.wd, "species_families_2021.csv",sep=""), header=T, sep=",", fill = TRUE) #load data
fam$species_matched<-gsub(" ", "_", fam$species_matched) #adapt species nomenclature
fam$family[fam$family=="Compositae"]<-"Asteraceae" #replace family name
fam$family[fam$family=="Leguminosae"]<-"Fabaceae" #replace family name
fam<-fam[which(fam$species_matched %in% rownames(dat)),] #subset only species included in our treatment

#get table with ranked families based on their number of species:
famf<-as.data.frame(table(fam$family)) #create dataframe
famf<-famf[order(-famf$Freq),] #order it
row.names(famf) <- NULL #get rid of rownames

#create vector with unique names of species:
list.r <- unique(fam$family)

#get nodes for families:
list.nod<-NULL
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$species_matched) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- findMRCA(tree2, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)


###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=24)) #select the top 24 families with 5 for more species

#Plot tree 
#remember to change angle = "auto" everywheree to avoid overlap in names. Consider also unifying "barsize" (to 0.1, for example):
p <- 
  ggtree(tree2, layout="circular", size=0.5)+ # build circular tree
  geom_point(aes(colour=as.factor(significant)), size=2, alpha=1, show.legend = TRUE) + # highlight nodes
  scale_colour_manual(values=c("red", "deepskyblue"), labels=c("Lower DCi", "Higher DCi"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  geom_cladelabel(node=subset(famf, Var1=="Poaceae")$num, label="Poaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asteraceae")$num, label="Asteraceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Fabaceae")$num, label="Fabaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 10) +
  geom_cladelabel(node=subset(famf, Var1=="Cyperaceae")$num, label="Cyperaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Rosaceae")$num, label="Rosaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 355) +
  
  geom_cladelabel(node=subset(famf, Var1=="Apiaceae")$num, label="Apiaceae",  fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Brassicaceae")$num, label="Brassicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 330) +
  geom_cladelabel(node=subset(famf, Var1=="Plantaginaceae")$num, label="Plantaginaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 85) +
  geom_cladelabel(node=subset(famf, Var1=="Caryophyllaceae")$num, label="Caryophyllaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ranunculaceae")$num, label="Ranunculaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 290) +
  
  geom_cladelabel(node=subset(famf, Var1=="Lamiaceae")$num, label="Lamiaceae",  fontsize=2.5, barsize = 0.1, angle ="auto") +
  geom_cladelabel(node=subset(famf, Var1=="Onagraceae")$num, label="Onagraceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 310) +
  geom_cladelabel(node=subset(famf, Var1=="Polygonaceae")$num, label="Polygonaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 20) +
  geom_cladelabel(node=subset(famf, Var1=="Gentianaceae")$num, label="Gentianaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 75) +
  geom_cladelabel(node=subset(famf, Var1=="Orobanchaceae")$num, label="Orobanchaceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1=="Euphorbiaceae")$num, label="Euphorbiaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 350) +
  geom_cladelabel(node=subset(famf, Var1=="Boraginaceae")$num, label="Boraginaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 83) +
  geom_cladelabel(node=subset(famf, Var1=="Amaranthaceae")$num, label="Amaranthaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ericaceae")$num, label="Ericaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Rubiaceae")$num, label="Rubiaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 65) +
  
  geom_cladelabel(node=subset(famf, Var1=="Juncaceae")$num, label="Juncaceae",  fontsize=2.5, barsize = 0.1,angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Solanaceae")$num, label="Solanaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 64) +
  geom_cladelabel(node=subset(famf, Var1=="Geraniaceae")$num, label="Geraniaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Apocynaceae")$num, label="Apocynaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 70) +
  geom_cladelabel(node=subset(famf, Var1=="Violaceae")$num, label="Violaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 340) +
  
  geom_cladelabel(node=subset(famf, Var1=="Amaryllidaceae")$num, label="Amaryllidaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asparagaceae")$num, label="Asparagaceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Campanulaceae")$num, label="Campanulaceae", fontsize=2.5, barsize=0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Saxifragaceae")$num, label="Saxifragaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Polemoniaceae")$num, label="Polemoniaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 55) +
  
  geom_cladelabel(node=subset(famf, Var1=="Caprifoliaceae")$num, label="Caprifoliaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Convolvulaceae")$num, label="Convolvulaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 65) +
  geom_cladelabel(node=subset(famf, Var1=="Malvaceae")$num, label="Malvaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1=="Primulaceae")$num, label="Primulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Salicaceae")$num, label="Salicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 340) +
  
  geom_cladelabel(node=subset(famf, Var1=="Verbenaceae")$num, label="Verbenaceae",  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 90) +
  geom_cladelabel(node=subset(famf, Var1=="Oxalidaceae")$num, label="Oxalidaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 335) +
  geom_cladelabel(node=subset(famf, Var1=="Iridaceae")$num, label="Iridaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Anacardiaceae")$num, label="Anacardiaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 315) +
  geom_cladelabel(node=subset(famf, Var1=="Betulaceae")$num, label="Betulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 347) +
  
  geom_cladelabel(node=subset(famf, Var1=="Fagaceae")$num, label="Fagaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 340) +
  geom_cladelabel(node=subset(famf, Var1=="Nyctaginaceae")$num, label="Nyctaginaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 20) +
  geom_cladelabel(node=subset(famf, Var1=="Orchidaceae")$num, label="Orchidaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Crassulaceae")$num, label="Crassulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Liliaceae")$num, label="Liliaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1=="Cactaceae")$num, label="Cactaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 40) +
  geom_cladelabel(node=subset(famf, Var1=="Cistaceae")$num, label="Cistaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 320) +
  geom_cladelabel(node=subset(famf, Var1=="Commelinaceae")$num, label="Commelinaceae", fontsize=2.5, barsize=0.5, angle = 'auto') +
  geom_cladelabel(node=subset(famf, Var1=="Equisetaceae")$num, label="Equisetaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Lycopodiaceae")$num, label="Lycopodiaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="bottom")

#save output:
png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_temp.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()

###
# Filter by treatment == herbivore removal
###

dat<-subset(species.data, trt_type2=="herb_removal")[,c(1,2)] #select "all mult" treatment from original data
dat<-aggregate(dat[, 2], list(dat$species_matched), mean, na.rm=T) #get mean DCi value per species - this is not necessary - is already the average.
dat<-dat[dat$Group.1 %in% tree$tip.label, ] #make sure all species in the data are on the tree
rownames(dat)<-dat$Group.1 #set species names as rownames
dat$Group.1<-NULL #and delete column with species names

#prune tree:
tree2<-keep.tip(tree, rownames(dat))


###
# Run function to calculate if each node has significantly higher or lower mean DCi than
# expected if phylogenetic relationships were at random (laod function in the "DCi_nodes_scorre.R" script)
###

res<-node.mean(tree2, dat, 999)
write.table(res, paste(my.wd, "res_phylo_herb_removal.csv", sep=""))
res<-read.table(paste(my.wd, "res_phylo_herb_removal.csv", sep=""))#save the result
#res2<-subset(res, P_value<0.01) #this would tell you what nodes are significant with alpha < 0.01
#tips(tree2, 1543) #and this would tell you what species are found in that clade

###
# Clean-up result to highlight nodes in the tree
###

significant<-res #create a copy of the main result
significant$P_value[significant$P_value>0.05]<-NA #replace non-significant with NA
significant$P_value[significant$SR<3]<-NA #assign NA to nodes with 2 or less species

significant$P_value[1]<-NA #set the first node (the root node) to NA
significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
significant<-c(rep(NA, length(tree2$tip.label)), significant$P_value) #merge "tip nodes" with "inner" tree nodes
significant<-as.factor(significant) #convert values into factors

significant<-factor(significant, levels = c("neg.05", "pos.05" )) #change order of factors



###
# Assign families to nodes
###

#load family data and clean-up:
fam<-read.table(paste(my.wd, "species_families_2021.csv",sep=""), header=T, sep=",", fill = TRUE) #load data
fam$species_matched<-gsub(" ", "_", fam$species_matched) #adapt species nomenclature
fam$family[fam$family=="Compositae"]<-"Asteraceae" #replace family name
fam$family[fam$family=="Leguminosae"]<-"Fabaceae" #replace family name
fam<-fam[which(fam$species_matched %in% rownames(dat)),] #subset only species included in our treatment

#get table with ranked families based on their number of species:
famf<-as.data.frame(table(fam$family)) #create dataframe
famf<-famf[order(-famf$Freq),] #order it
row.names(famf) <- NULL #get rid of rownames

#create vector with unique names of species:
list.r <- unique(fam$family)

#get nodes for families:
list.nod<-NULL
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$species_matched) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- findMRCA(tree2, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)


###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=24)) #select the top 24 families with 5 for more species

#Plot tree 
#remember to change angle = "auto" everywheree to avoid overlap in names. Consider also unifying "barsize" (to 0.1, for example):
p <- 
  ggtree(tree2, layout="circular", size=0.5)+ # build circular tree
  geom_point(aes(colour=as.factor(significant)), size=2, alpha=1, show.legend = TRUE) + # highlight nodes
  scale_colour_manual(values=c("red", "deepskyblue"), labels=c("Lower DCi", "Higher DCi"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  geom_cladelabel(node=subset(famf, Var1=="Poaceae")$num, label="Poaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asteraceae")$num, label="Asteraceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Fabaceae")$num, label="Fabaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 5) +
  geom_cladelabel(node=subset(famf, Var1=="Cyperaceae")$num, label="Cyperaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Rosaceae")$num, label="Rosaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 340) +
  
  geom_cladelabel(node=subset(famf, Var1=="Apiaceae")$num, label="Apiaceae",  fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Brassicaceae")$num, label="Brassicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 320) +
  geom_cladelabel(node=subset(famf, Var1=="Plantaginaceae")$num, label="Plantaginaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 85) +
  geom_cladelabel(node=subset(famf, Var1=="Caryophyllaceae")$num, label="Caryophyllaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ranunculaceae")$num, label="Ranunculaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 305) +
  
  geom_cladelabel(node=subset(famf, Var1=="Lamiaceae")$num, label="Lamiaceae",  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 90) +
  geom_cladelabel(node=subset(famf, Var1=="Onagraceae")$num, label="Onagraceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 310) +
  geom_cladelabel(node=subset(famf, Var1=="Polygonaceae")$num, label="Polygonaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 15) +
  geom_cladelabel(node=subset(famf, Var1=="Gentianaceae")$num, label="Gentianaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 75) +
  geom_cladelabel(node=subset(famf, Var1=="Orobanchaceae")$num, label="Orobanchaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 88) +
  
  geom_cladelabel(node=subset(famf, Var1=="Euphorbiaceae")$num, label="Euphorbiaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1=="Boraginaceae")$num, label="Boraginaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 83) +
  geom_cladelabel(node=subset(famf, Var1=="Amaranthaceae")$num, label="Amaranthaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ericaceae")$num, label="Ericaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Rubiaceae")$num, label="Rubiaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 65) +
  
  geom_cladelabel(node=subset(famf, Var1=="Juncaceae")$num, label="Juncaceae",  fontsize=2.5, barsize = 0.1,angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Solanaceae")$num, label="Solanaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 64) +
  geom_cladelabel(node=subset(famf, Var1=="Geraniaceae")$num, label="Geraniaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 308) +
  geom_cladelabel(node=subset(famf, Var1=="Apocynaceae")$num, label="Apocynaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 70) +
  geom_cladelabel(node=subset(famf, Var1=="Violaceae")$num, label="Violaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 335) +
  
  geom_cladelabel(node=subset(famf, Var1=="Amaryllidaceae")$num, label="Amaryllidaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asparagaceae")$num, label="Asparagaceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Campanulaceae")$num, label="Campanulaceae", fontsize=2.5, barsize=0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Saxifragaceae")$num, label="Saxifragaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Polemoniaceae")$num, label="Polemoniaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 55) +
  
  geom_cladelabel(node=subset(famf, Var1=="Caprifoliaceae")$num, label="Caprifoliaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Convolvulaceae")$num, label="Convolvulaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 65) +
  geom_cladelabel(node=subset(famf, Var1=="Malvaceae")$num, label="Malvaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 315) +
  geom_cladelabel(node=subset(famf, Var1=="Primulaceae")$num, label="Primulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Salicaceae")$num, label="Salicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 340) +
  
  geom_cladelabel(node=subset(famf, Var1=="Verbenaceae")$num, label="Verbenaceae",  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 80) +
  geom_cladelabel(node=subset(famf, Var1=="Oxalidaceae")$num, label="Oxalidaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1=="Iridaceae")$num, label="Iridaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Anacardiaceae")$num, label="Anacardiaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 315) +
  geom_cladelabel(node=subset(famf, Var1=="Betulaceae")$num, label="Betulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 347) +
  
  geom_cladelabel(node=subset(famf, Var1=="Fagaceae")$num, label="Fagaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 340) +
  geom_cladelabel(node=subset(famf, Var1=="Nyctaginaceae")$num, label="Nyctaginaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 20) +
  geom_cladelabel(node=subset(famf, Var1=="Orchidaceae")$num, label="Orchidaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Crassulaceae")$num, label="Crassulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Liliaceae")$num, label="Liliaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1=="Cactaceae")$num, label="Cactaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 40) +
  geom_cladelabel(node=subset(famf, Var1=="Cistaceae")$num, label="Cistaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 280) +
  geom_cladelabel(node=subset(famf, Var1=="Commelinaceae")$num, label="Commelinaceae", fontsize=2.5, barsize=0.5, angle = 'auto') +
  geom_cladelabel(node=subset(famf, Var1=="Equisetaceae")$num, label="Equisetaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Lycopodiaceae")$num, label="Lycopodiaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="bottom")

#save output:
png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_herb_removal.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()


###
# Filter by treatment == disturbance
###

dat<-subset(species.data, trt_type2=="dist")[,c(1,2)] #select "all mult" treatment from original data
dat<-aggregate(dat[, 2], list(dat$species_matched), mean, na.rm=T) #get mean DCi value per species - this is not necessary - is already the average.
dat<-dat[dat$Group.1 %in% tree$tip.label, ] #make sure all species in the data are on the tree
rownames(dat)<-dat$Group.1 #set species names as rownames
dat$Group.1<-NULL #and delete column with species names

#prune tree:
tree2<-keep.tip(tree, rownames(dat))


###
# Run function to calculate if each node has significantly higher or lower mean DCi than
# expected if phylogenetic relationships were at random (laod function in the "DCi_nodes_scorre.R" script)
###

res<-node.mean(tree2, dat, 999)
write.table(res, paste(my.wd, "res_phylo_disturbance.csv", sep=""))
res<-read.table(paste(my.wd, "res_phylo_disturbance.csv", sep=""))#save the result
#res2<-subset(res, P_value<0.01) #this would tell you what nodes are significant with alpha < 0.01
#tips(tree2, 1543) #and this would tell you what species are found in that clade

###
# Clean-up result to highlight nodes in the tree
###

significant<-res #create a copy of the main result
significant$P_value[significant$P_value>0.05]<-NA #replace non-significant with NA
significant$P_value[significant$SR<3]<-NA #assign NA to nodes with 2 or less species

significant$P_value[1]<-NA #set the first node (the root node) to NA
significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
significant<-c(rep(NA, length(tree2$tip.label)), significant$P_value) #merge "tip nodes" with "inner" tree nodes
significant<-as.factor(significant) #convert values into factors

significant<-factor(significant, levels = c("neg.05", "pos.05" )) #change order of factors



###
# Assign families to nodes
###

#load family data and clean-up:
fam<-read.table(paste(my.wd, "species_families_2021.csv",sep=""), header=T, sep=",", fill = TRUE) #load data
fam$species_matched<-gsub(" ", "_", fam$species_matched) #adapt species nomenclature
fam$family[fam$family=="Compositae"]<-"Asteraceae" #replace family name
fam$family[fam$family=="Leguminosae"]<-"Fabaceae" #replace family name
fam<-fam[which(fam$species_matched %in% rownames(dat)),] #subset only species included in our treatment

#get table with ranked families based on their number of species:
famf<-as.data.frame(table(fam$family)) #create dataframe
famf<-famf[order(-famf$Freq),] #order it
row.names(famf) <- NULL #get rid of rownames

#create vector with unique names of species:
list.r <- unique(fam$family)

#get nodes for families:
list.nod<-NULL
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$species_matched) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- findMRCA(tree2, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)


###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=24)) #select the top 24 families with 5 for more species

#Plot tree 
#remember to change angle = "auto" everywheree to avoid overlap in names. Consider also unifying "barsize" (to 0.1, for example):
p <- 
  ggtree(tree2, layout="circular", size=0.5)+ # build circular tree
  geom_point(aes(colour=as.factor(significant)), size=2, alpha=1, show.legend = TRUE) + # highlight nodes
  scale_colour_manual(values=c("red", "deepskyblue"), labels=c("Lower DCi", "Higher DCi"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  geom_cladelabel(node=subset(famf, Var1=="Poaceae")$num, label="Poaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asteraceae")$num, label="Asteraceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Fabaceae")$num, label="Fabaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 5) +
  geom_cladelabel(node=subset(famf, Var1=="Cyperaceae")$num, label="Cyperaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Rosaceae")$num, label="Rosaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 340) +
  
  geom_cladelabel(node=subset(famf, Var1=="Apiaceae")$num, label="Apiaceae",  fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Brassicaceae")$num, label="Brassicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 320) +
  geom_cladelabel(node=subset(famf, Var1=="Plantaginaceae")$num, label="Plantaginaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 85) +
  geom_cladelabel(node=subset(famf, Var1=="Caryophyllaceae")$num, label="Caryophyllaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ranunculaceae")$num, label="Ranunculaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 305) +
  
  geom_cladelabel(node=subset(famf, Var1=="Lamiaceae")$num, label="Lamiaceae",  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 90) +
  geom_cladelabel(node=subset(famf, Var1=="Onagraceae")$num, label="Onagraceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 310) +
  geom_cladelabel(node=subset(famf, Var1=="Polygonaceae")$num, label="Polygonaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 20) +
  geom_cladelabel(node=subset(famf, Var1=="Gentianaceae")$num, label="Gentianaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 75) +
  geom_cladelabel(node=subset(famf, Var1=="Orobanchaceae")$num, label="Orobanchaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 88) +
  
  geom_cladelabel(node=subset(famf, Var1=="Euphorbiaceae")$num, label="Euphorbiaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1=="Boraginaceae")$num, label="Boraginaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 83) +
  geom_cladelabel(node=subset(famf, Var1=="Amaranthaceae")$num, label="Amaranthaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ericaceae")$num, label="Ericaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Rubiaceae")$num, label="Rubiaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 65) +
  
  geom_cladelabel(node=subset(famf, Var1=="Juncaceae")$num, label="Juncaceae",  fontsize=2.5, barsize = 0.1,angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Solanaceae")$num, label="Solanaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 64) +
  geom_cladelabel(node=subset(famf, Var1=="Geraniaceae")$num, label="Geraniaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 308) +
  geom_cladelabel(node=subset(famf, Var1=="Apocynaceae")$num, label="Apocynaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 70) +
  geom_cladelabel(node=subset(famf, Var1=="Violaceae")$num, label="Violaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 335) +
  
  geom_cladelabel(node=subset(famf, Var1=="Amaryllidaceae")$num, label="Amaryllidaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asparagaceae")$num, label="Asparagaceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Campanulaceae")$num, label="Campanulaceae", fontsize=2.5, barsize=0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Saxifragaceae")$num, label="Saxifragaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Polemoniaceae")$num, label="Polemoniaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 55) +
  
  geom_cladelabel(node=subset(famf, Var1=="Caprifoliaceae")$num, label="Caprifoliaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Convolvulaceae")$num, label="Convolvulaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 65) +
  geom_cladelabel(node=subset(famf, Var1=="Malvaceae")$num, label="Malvaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 315) +
  geom_cladelabel(node=subset(famf, Var1=="Primulaceae")$num, label="Primulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Salicaceae")$num, label="Salicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 340) +
  
  geom_cladelabel(node=subset(famf, Var1=="Verbenaceae")$num, label="Verbenaceae",  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 80) +
  geom_cladelabel(node=subset(famf, Var1=="Oxalidaceae")$num, label="Oxalidaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1=="Iridaceae")$num, label="Iridaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Anacardiaceae")$num, label="Anacardiaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 315) +
  geom_cladelabel(node=subset(famf, Var1=="Betulaceae")$num, label="Betulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 347) +
  
  geom_cladelabel(node=subset(famf, Var1=="Fagaceae")$num, label="Fagaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 340) +
  geom_cladelabel(node=subset(famf, Var1=="Nyctaginaceae")$num, label="Nyctaginaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 20) +
  geom_cladelabel(node=subset(famf, Var1=="Orchidaceae")$num, label="Orchidaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Crassulaceae")$num, label="Crassulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Liliaceae")$num, label="Liliaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1=="Cactaceae")$num, label="Cactaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 40) +
  geom_cladelabel(node=subset(famf, Var1=="Cistaceae")$num, label="Cistaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 280) +
  geom_cladelabel(node=subset(famf, Var1=="Commelinaceae")$num, label="Commelinaceae", fontsize=2.5, barsize=0.5, angle = 'auto') +
  geom_cladelabel(node=subset(famf, Var1=="Equisetaceae")$num, label="Equisetaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Lycopodiaceae")$num, label="Lycopodiaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="bottom")

#save output:
png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_disturbance.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()

###
# Filter by treatment = temp_other
###

dat<-subset(species.data, trt_type2=="temp_other")[,c(1,2)] #select "all mult" treatment from original data
dat<-aggregate(dat[, 2], list(dat$species_matched), mean, na.rm=T) #get mean DCi value per species
dat<-dat[dat$Group.1 %in% tree$tip.label, ] #make sure all species in the data are on the tree
rownames(dat)<-dat$Group.1 #set species names as rownames
dat$Group.1<-NULL #and delete column with species names

#prune tree:
tree2<-keep.tip(tree, rownames(dat))


###
# Run function to calculate if each node has significantly higher or lower mean DCi than
# expected if phylogenetic relationships were at random (laod function in the "DCi_nodes_scorre.R" script)
###

res<-node.mean(tree2, dat, 999)
#write.table(res, paste(my.wd, "res_phylo_p_other.csv", sep="")) #save the result
#res<-read.table(paste(my.wd, "res_phylo_p_other.csv", sep=""))
#res2<-subset(res, P_value<0.01) #this would tell you what nodes are significant with alpha < 0.01
#tips(tree2, 1543) #and this would tell you what species are found in that clade

###
# Clean-up result to highlight nodes in the tree
###

significant<-res #create a copy of the main result
significant$P_value[significant$P_value>0.05]<-NA #replace non-significant with NA
significant$P_value[significant$SR<3]<-NA #assign NA to nodes with 2 or less species
significant$P_value[1]<-NA #set the first node (the root node) to NA
significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
significant<-c(rep(NA, length(tree2$tip.label)), significant$P_value) #merge "tip nodes" with "inner" tree nodes
significant<-as.factor(significant) #convert values into factors
significant<-factor(significant, levels = c("neg.05", "pos.05" )) #change order of factors


###
# Assign families to nodes
###

#load family data and clean-up:
fam<-read.table(paste(my.wd, "species_families_2021.csv",sep=""), header=T, sep=",", fill = TRUE) #load data
fam$species_matched<-gsub(" ", "_", fam$species_matched) #adapt species nomenclature
fam$family[fam$family=="Compositae"]<-"Asteraceae" #replace family name
fam$family[fam$family=="Leguminosae"]<-"Fabaceae" #replace family name
fam<-fam[which(fam$species_matched %in% rownames(dat)),] #subset only species included in our treatment

#get table with ranked families based on their number of species:
famf<-as.data.frame(table(fam$family)) #create dataframe
famf<-famf[order(-famf$Freq),] #order it
row.names(famf) <- NULL #get rid of rownames

#create vector with unique names of species:
list.r <- unique(fam$family)

#get nodes for families:
list.nod<-NULL
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$species_matched) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- findMRCA(tree2, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)


###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=38)) #select the top 38 families with 5 or more species

#Plot tree for N 
#remember to change angle = "auto" everywheree to avoid overlap in names. Consider also unifying "barsize" (to 0.1, for example):
p <- 
  ggtree(tree2, layout="circular", size=0.5)+ # build circular tree
  geom_point(aes(colour=as.factor(significant)), size=2, alpha=1, show.legend = TRUE) + # highlight nodes
  scale_colour_manual(values=c("red", "deepskyblue"), labels=c("Lower DCi", "Higher DCi"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  geom_cladelabel(node=subset(famf, Var1=="Poaceae")$num, label="Poaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asteraceae")$num, label="Asteraceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Fabaceae")$num, label="Fabaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 5) +
  geom_cladelabel(node=subset(famf, Var1=="Cyperaceae")$num, label="Cyperaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Rosaceae")$num, label="Rosaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 354) +
  
  geom_cladelabel(node=subset(famf, Var1=="Apiaceae")$num, label="Apiaceae",  fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Brassicaceae")$num, label="Brassicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 320) +
  geom_cladelabel(node=subset(famf, Var1=="Plantaginaceae")$num, label="Plantaginaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 85) +
  geom_cladelabel(node=subset(famf, Var1=="Caryophyllaceae")$num, label="Caryophyllaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 45) +
  geom_cladelabel(node=subset(famf, Var1=="Ranunculaceae")$num, label="Ranunculaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 294) +
  
  geom_cladelabel(node=subset(famf, Var1=="Lamiaceae")$num, label="Lamiaceae",  fontsize=2.5, barsize = 0.1, angle ="auto") +
  geom_cladelabel(node=subset(famf, Var1=="Onagraceae")$num, label="Onagraceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 307) +
  geom_cladelabel(node=subset(famf, Var1=="Polygonaceae")$num, label="Polygonaceae", fontsize=2.5, barsize=0.1, hjust= 1, angle = 35) +
  geom_cladelabel(node=subset(famf, Var1=="Gentianaceae")$num, label="Gentianaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 80) +
  geom_cladelabel(node=subset(famf, Var1=="Orobanchaceae")$num, label="Orobanchaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 90) +
  
  geom_cladelabel(node=subset(famf, Var1=="Euphorbiaceae")$num, label="Euphorbiaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1=="Boraginaceae")$num, label="Boraginaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 83) +
  geom_cladelabel(node=subset(famf, Var1=="Amaranthaceae")$num, label="Amaranthaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 42) +
  geom_cladelabel(node=subset(famf, Var1=="Ericaceae")$num, label="Ericaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 50) +
  geom_cladelabel(node=subset(famf, Var1=="Rubiaceae")$num, label="Rubiaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 65) +
  
  geom_cladelabel(node=subset(famf, Var1=="Juncaceae")$num, label="Juncaceae",  fontsize=2.5, barsize = 0.1,angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Solanaceae")$num, label="Solanaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 64) +
  geom_cladelabel(node=subset(famf, Var1=="Geraniaceae")$num, label="Geraniaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 303) +
  geom_cladelabel(node=subset(famf, Var1=="Apocynaceae")$num, label="Apocynaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 70) +
  geom_cladelabel(node=subset(famf, Var1=="Violaceae")$num, label="Violaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 320) +
  
  geom_cladelabel(node=subset(famf, Var1=="Amaryllidaceae")$num, label="Amaryllidaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Asparagaceae")$num, label="Asparagaceae", fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Campanulaceae")$num, label="Campanulaceae", fontsize=2.5, barsize=0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Polemoniaceae")$num, label="Polemoniaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 51) +
  
  geom_cladelabel(node=subset(famf, Var1=="Caprifoliaceae")$num, label="Caprifoliaceae",  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Convolvulaceae")$num, label="Convolvulaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 60) +
  geom_cladelabel(node=subset(famf, Var1=="Malvaceae")$num, label="Malvaceae", fontsize=2.5, barsize=0.5, hjust= 1, angle = 315) +
  geom_cladelabel(node=subset(famf, Var1=="Primulaceae")$num, label="Primulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 55) +
  geom_cladelabel(node=subset(famf, Var1=="Salicaceae")$num, label="Salicaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 320) +
  
  geom_cladelabel(node=subset(famf, Var1=="Verbenaceae")$num, label="Verbenaceae",  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 90) +
  geom_cladelabel(node=subset(famf, Var1=="Oxalidaceae")$num, label="Oxalidaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 322) +
  geom_cladelabel(node=subset(famf, Var1=="Iridaceae")$num, label="Iridaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Anacardiaceae")$num, label="Anacardiaceae", fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 310) +
  geom_cladelabel(node=subset(famf, Var1=="Betulaceae")$num, label="Betulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 340) +
  
  geom_cladelabel(node=subset(famf, Var1=="Nyctaginaceae")$num, label="Nyctaginaceae", fontsize=2.5, barsize = 0.1, hjust= 1, angle = 35) +
  geom_cladelabel(node=subset(famf, Var1=="Orchidaceae")$num, label="Orchidaceae", fontsize=2.5, barsize=0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Crassulaceae")$num, label="Crassulaceae", fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1=="Liliaceae")$num, label="Liliaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1=="Cactaceae")$num, label="Cactaceae",  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 40) +
  geom_cladelabel(node=subset(famf, Var1=="Cistaceae")$num, label="Cistaceae", fontsize=2.5, barsize = 0.5, hjust= 1, angle = 312) +
  geom_cladelabel(node=subset(famf, Var1=="Commelinaceae")$num, label="Commelinaceae", fontsize=2.5, barsize=0.5, angle = 'auto') +
  geom_cladelabel(node=subset(famf, Var1=="Equisetaceae")$num, label="Equisetaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1=="Lycopodiaceae")$num, label="Lycopodiaceae", fontsize=2.5, barsize = 0.5,  angle = "auto") +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="bottom")

#save output:
png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_temp_other.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()

