#sCorre species response DCiDiff on the phylogeny and phylogenetic signal#

#Clear work space
rm(list=ls())

#Set working directory
setwd("~/Dropbox (iDiv)/My Mac (idivmac32.local)/Documents/Project sCorre/winners_losers")
setwd("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/")

#Load library()
library(geiger)
library(ape)

#Load data
tree<-read.tree("scorre.tree.S3.tre") #you can find in the Dropbox (Phylogenetic data)
dd<-read.csv("Species_DCiDiff.csv")
str(dd)
str(tree)

#####
# Create Phylo-ring to show patterns on tips
#####

#load packages:

library(phytools)
library(ggplot2)
library(stringr)
library(plyr)
library(ape)
library(scico)



#data manupulation:
dd$species_matched <- revalue(dd$species_matched, c("Aronia x"="Aronia x prunifolia"))
species.data<-subset(dd, trt_type2=="drought") #subset treatment to make it easier.

#Treatments:
#"overall"       "N"             "mult_nutrient" "irr"           "drought"


#remove from the original table non-vascular plants that couldn't be added to the tree:
species.data$species_matched<-gsub(" ", "_", species.data$species_matched) #unify nomenclature
in.data.not.tree <- setdiff(unique(species.data$species_matched), tree$tip.label)

#unlock the following line if "in.data.not.tree" is not empty (only applies to "drought").
#species.data <- species.data[-which(species.data$species_matched %in% in.data.not.tree),] 

#rearrange original table:
df<-unique(species.data[,c(1,5)])
rownames(df)<-df$species_matched
df$species_matched<-NULL

x<-name.check(tree, df)
tree2<-drop.tip(tree,c(x$tree_not_data))
x<-name.check(tree2, df)
scorre.tree<-tree2

################################################

#Plot on the tree
tree2<-ladderize(tree2)
plotTree.barplot(tree2,df[,1,drop=FALSE],
                 args.plotTree=list(ftype="off"),
                 args.barplot=list(xlab="DCiDiff",space=0.5))


df2<-df[,1]
names(df2)<-row.names(df)

#Phylosignal with Blombergs K
phylosig(tree2, df2, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
         control=list())

#Phylosignal with lambda
phylosig(tree2, df2, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
         control=list())

#Treatment overall
Phylogenetic signal K : 0.00295202 
P-value (based on 1000 randomizations) : 0.445 
Phylogenetic signal lambda : 0.139708 
logL(lambda) : 1979.19 
LR(lambda=0) : 8.50465 
P-value (based on LR test) : 0.0035424 

#Treatment N
Phylogenetic signal K : 0.00214985 
P-value (based on 1000 randomizations) : 0.923 
Phylogenetic signal lambda : 0.399009 
logL(lambda) : 910.299 
LR(lambda=0) : 28.2253 
P-value (based on LR test) : 1.07986e-07

#Treatment mult_nutrient
Phylogenetic signal K : 0.00422832 
P-value (based on 1000 randomizations) : 0.161 
Phylogenetic signal lambda : 0.239786 
logL(lambda) : 901.218 
LR(lambda=0) : 14.529 
P-value (based on LR test) : 0.000138019 

#Treatment irr
Phylogenetic signal K : 0.0132178 
P-value (based on 1000 randomizations) : 0.675 
Phylogenetic signal lambda : 6.13209e-05 
logL(lambda) : 438.909 
LR(lambda=0) : -0.0230268 
P-value (based on LR test) : 1 

#Treatment drought
Phylogenetic signal K : 0.00690762 
P-value (based on 1000 randomizations) : 0.96 
Phylogenetic signal lambda : 0.0252875 
logL(lambda) : 437.882 
LR(lambda=0) : 1.19607 
P-value (based on LR test) : 0.274108