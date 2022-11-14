
####
# Script to calculate Phylogenetic Signal
###

#Clear work space
rm(list=ls())

#set directory:
#my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/"
my.wd <- "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"
my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/"
#my.wd <- "E:/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/"

#Load library()
library(geiger)
library(phangorn)
library(ape)
library(phytools)
library(ggplot2)
library(stringr)
library(plyr)
library(ape)
library(dplyr)

#sCorre species response DCiDiff on the phylogeny and phylogenetic signal#

#Load data
tree<-read.tree(paste(my.wd, "scorre.tree.win.los.tre.dec2021", sep="")) #load tree
dd<-read.csv(paste(my.wd,"Species_DCiDiff_Nov2022.csv",sep=""), header=T)

dd$species_matched<-gsub(" ", "_", dd$species_matched) #unify nomenclature
str(dd)
str(tree)

#get levels of treatments:
trt<-levels(as.factor(dd$trt_type2))


#re-arrange table:
final.dd<-as.data.frame(unique(dd$species_matched))
names(final.dd)[1]<-paste("species_matched")

for (i in 1:length(trt))
{
  sub.dd<-subset(dd, trt_type2==trt[i])[,c(1,2)]
  final.dd<-merge(final.dd, sub.dd, by="species_matched", all.x=T)
  colnames(final.dd)[i+1]<-trt[i]
}


rownames(final.dd)<-final.dd$species_matched
final.dd$species_matched<-NULL


######
#calculate phylogenetic signal:

#prepare data:
treat<-tree$tip.label
treat<-final.dd[match(treat, rownames(final.dd)),]
rownames(treat)<-NULL

output<-data.frame(treatment=trt, blombergK =rep(NA, length(trt)), blombergK_pval=rep(NA, length(trt)),
                   lambda =rep(NA, length(trt)), lambda_pval=rep(NA, length(trt)))
for (i in 1:length(trt))
{
  blomb<-phylosig(tree, treat[,i], method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                  control=list())
  output[i,2]<-blomb$K
  output[i,3]<-blomb$P
  blomb<-phylosig(tree, treat[,i], method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
                  control=list())
  output[i,4]<-blomb$lambda
  output[i,5]<-blomb$P
}

#write.csv(output, paste(my.wd,"phylo_signal.csv",sep=""), row.names = F)







###
# Repeat for input traits
###

tree<-read.tree(paste(my.wd, "scorre.tree.win.los.tre.dec2021", sep="")) #load tree
dd<-read.table(paste(my.wd,"Final_Input_sCorr.csv",sep=""), header=T, sep=",")
dd$species_matched<-gsub(" ", "_", dd$species_matched)

#remove from the original table non-vascular plants that couldn't be added to the tree:
in.data.not.tree <- setdiff(unique(dd$species_matched), tree$tip.label)
dd <- dd[-which(dd$species_matched %in% in.data.not.tree),]
tree<-keep.tip(tree, unique(dd$species_matched))

#get mean trait values:
#sapply(dd[,c(6:19)], function(x) sum(is.na(x)))/nrow(dd)
dd<- dd %>% group_by(species_matched) %>% summarise_at(vars("seed_dry_mass", "leaf_N", "LDMC",
                                                            "leaf_dry_mass", "plant_height_vegetative",
                                                            "SLA"), mean, na.rm=T)

#re-arrange data:
dd<-as.data.frame(dd)
rownames(dd)<-dd$species_matched
dd$species_matched<-NULL

######
#calculate phylogenetic signal:

#prepare data:
treat<-tree$tip.label
treat<-dd[match(treat, rownames(dd)),]
rownames(treat)<-NULL

output<-data.frame(treatment=colnames(treat), blombergK =rep(NA, ncol(treat)), blombergK_pval=rep(NA, ncol(treat)),
                   lambda =rep(NA, ncol(treat)), lambda_pval=rep(NA, ncol(treat)))
for (i in 1:ncol(treat))
{
  blomb<-phylosig(tree, treat[,i], method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                  control=list())
  output[i,2]<-blomb$K
  output[i,3]<-blomb$P
  blomb<-phylosig(tree, treat[,i], method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
                  control=list())
  output[i,4]<-blomb$lambda
  output[i,5]<-blomb$P
}

write.table(output, paste(my.wd,"phylo_signal_traits.csv",sep=""))


#clean-up:
rm(list = ls())

