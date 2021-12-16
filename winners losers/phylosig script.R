#sCorre species response DCiDiff on the phylogeny and phylogenetic signal#

#Clear work space
rm(list=ls())

#set directory:
#my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/"
my.wd <- "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"
my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/"
my.wd <- "E:/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/"

#Load library()
library(geiger)
library(phangorn)
library(ape)
library(phytools)
library(ggplot2)
library(stringr)
library(plyr)
library(ape)

#Load data
tree<-read.tree(paste(my.wd, "scorre.tree.win.los.tre.dec2021", sep="")) #load tree
dd<-read.table(paste(my.wd,"Species_DCiDiff_newtrts_filtered_Dec2021.csv",sep=""), header=T)
str(dd)
str(tree)

#get levels of treatments:
trt<-levels(as.factor(dd$trt_type2))

#[1] "all mult"       "co2"            "co2_other"      "dist_other"     "disturbance"    "drought"        "drt_other"     
#[8] "herb_rem_other" "herb_removal"   "irg_other"      "irrigation"     "n"              "n_other"        "p"             
#[15] "p_other"        "temp"           "temp_other"

#re-arrange table:
final.dd<-as.data.frame(unique(dd$species_matched))
names(final.dd)[1]<-paste("species_matched")

for (i in 1:length(trt))
{
  sub.dd<-subset(dd, trt_type2==trt[i])[,c(1,4)]
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

getwd()
write.table(output, paste(my.wd,"phylo_signal.csv",sep=""))

#clean-up:
rm(list = ls())

