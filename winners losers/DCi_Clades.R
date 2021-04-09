
#load libraries:
library(ggtree)
library(ggplot2)
library(stringr)
library(plyr)
library(ape)
library(scico)
library(phytools)

#set directory.
my.wd<-"/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"

###
# Load data
###

#load species data (filtered after removing mosses and species missing from the phylogeny):
species.data<-read.table(paste(my.wd,"Species_DCiDiff_newtrts_filtered.csv",sep=""), header=T, sep=" ")

#load phylogenetic tree:
tree<-read.tree(paste(my.wd, "scorre.tree.win.los.tre", sep=""))


###
# Filter by treatment
###

dat<-subset(species.data, trt_type2=="all mult")[,c(1,4)] #select "all mult" treatment from original data
dat<-aggregate(dat[, 2], list(dat$species_matched), mean, na.rm=T) #get mean DCi value per species
dat<-dat[dat$Group.1 %in% tree$tip.label, ] #make sure all species in the data are on the tree
rownames(dat)<-dat$Group.1 #set species names as rownames
dat$Group.1<-NULL #and delete column with species names

#prune tree:
tree2<-keep.tip(tree, rownames(dat))


###
# Run function to calculate if each node has significantly higher or lower median DCi than
# expected if phylogenetic relationships were at random (laod function in the "DCi_nodes_scorre.R" script)
###

res<-node.mean(tree2, dat, 999)
write.table(res, paste(my.wd, "res_phylo_all_mult.csv", sep="")) #save the result
#res2<-subset(res, P_value<0.01) #this would tell you what nodes are significant with alpha < 0.01
#tips(tree2, 1543) #and this would tell you what species are found in that clade

###
# Clean-up result to highlight nodes in the tree
###

significant<-res #create a copy of the main result
significant$P_value[significant$P_value>0.05 | significant$SD_Exp<0.001]<-NA #replace non-significant with NA
significant$P_value[significant$SR<3]<-NA #assign NA to nodes with 2 or less species
significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
significant<-c(rep(NA, length(tree2$tip.label)), significant$P_value) #merge "tip nodes" with "inner" tree nodes
significant<-as.factor(significant) #convert values into factors
significant<-factor(significant, levels = c("neg.05", "pos.05")) #change order of factors


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
toplot<-as.character(head(famf$Var1, n=48)) #select the top 48 families

#Plot tree 
#remember to change angle = "auto" everywheree to avoid overlap in names. Consider also unifying "barsize" (to 0.1, for example):
p <- ggtree(tree2, layout="circular", size=0.5, branch.letngth="none")+ # build circular tree
  geom_point(aes(colour=as.factor(significant)), size=2, alpha=1, show.legend = TRUE) + # highlight nodes
  scale_colour_manual(values=c("red", "deepskyblue"), labels=c(expression(atop("Lower DCi", italic("(P<0.05)"))),
                                                                                        expression(atop("Higher DCi", italic("(P<0.05)")))),
                      na.translate=FALSE)+ # set aesthetics for highlighted nodes
  geom_cladelabel(node=subset(famf, Var1==toplot[1])$num, label=toplot[1], fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[2])$num, label=toplot[2], fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[3])$num, label=toplot[3], fontsize=2.5, barsize=0.5, hjust= 1, angle = 16) +
  geom_cladelabel(node=subset(famf, Var1==toplot[4])$num, label=toplot[4], fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[5])$num, label=toplot[5], fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 354) +

  geom_cladelabel(node=subset(famf, Var1==toplot[6])$num, label=toplot[6],  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 323) +
  geom_cladelabel(node=subset(famf, Var1==toplot[7])$num, label=toplot[7], fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[8])$num, label=toplot[8], fontsize=2.5, barsize=0.1, hjust= 1, angle = 89) +
  geom_cladelabel(node=subset(famf, Var1==toplot[9])$num, label=toplot[9], fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 309) +
  geom_cladelabel(node=subset(famf, Var1==toplot[10])$num, label=toplot[10], fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 294) +

  geom_cladelabel(node=subset(famf, Var1==toplot[11])$num, label=toplot[11],  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 48) +
  geom_cladelabel(node=subset(famf, Var1==toplot[12])$num, label=toplot[12], fontsize=2.5, barsize = 0.5, hjust= 1, angle = 339) +
  geom_cladelabel(node=subset(famf, Var1==toplot[13])$num, label=toplot[13], fontsize=2.5, barsize=0.1, hjust= 1, angle = 35) +
  geom_cladelabel(node=subset(famf, Var1==toplot[14])$num, label=toplot[14], fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 81) +
  geom_cladelabel(node=subset(famf, Var1==toplot[15])$num, label=toplot[15], fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 76) +

  geom_cladelabel(node=subset(famf, Var1==toplot[16])$num, label=toplot[16],  fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[17])$num, label=toplot[17], fontsize=2.5, barsize = 0.5, hjust= 1, angle = 42) +
  geom_cladelabel(node=subset(famf, Var1==toplot[18])$num, label=toplot[18], fontsize=2.5, barsize=0.5, hjust= 1, angle = 54) +
  geom_cladelabel(node=subset(famf, Var1==toplot[19])$num, label=toplot[19], fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 68) +
  geom_cladelabel(node=subset(famf, Var1==toplot[20])$num, label=toplot[20], fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 64) +

  geom_cladelabel(node=subset(famf, Var1==toplot[21])$num, label=toplot[21],  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 304) +
  geom_cladelabel(node=subset(famf, Var1==toplot[22])$num, label=toplot[22], fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[23])$num, label=toplot[23], fontsize=2.5, barsize=0.5, hjust= 1, angle = 72) +
  geom_cladelabel(node=subset(famf, Var1==toplot[24])$num, label=toplot[24], fontsize=2.5, barsize = 0.1,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[25])$num, label=toplot[25], fontsize=2.5, barsize = 0.1,  angle = "auto") +

  geom_cladelabel(node=subset(famf, Var1==toplot[26])$num, label=toplot[26],  fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[27])$num, label=toplot[27], fontsize=2.5, barsize = 0.1, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[28])$num, label=toplot[28], fontsize=2.5, barsize=0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[29])$num, label=toplot[29], fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[30])$num, label=toplot[30], fontsize=2.5, barsize = 0.1,  angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[31])$num, label=toplot[31],  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 318) +
  geom_cladelabel(node=subset(famf, Var1==toplot[32])$num, label=toplot[32], fontsize=2.5, barsize = 0.1, hjust= 1, angle = 331) +
  geom_cladelabel(node=subset(famf, Var1==toplot[33])$num, label=toplot[33], fontsize=2.5, barsize=0.5, hjust= 1, angle = 59) +
  geom_cladelabel(node=subset(famf, Var1==toplot[34])$num, label=toplot[34], fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 37) +
  geom_cladelabel(node=subset(famf, Var1==toplot[35])$num, label=toplot[35], fontsize=2.5, barsize = 0.1, hjust= 1, angle = 57) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[36])$num, label=toplot[36],  fontsize=2.5, barsize = 0.1, hjust= 1, angle = 314) +
  geom_cladelabel(node=subset(famf, Var1==toplot[37])$num, label=toplot[37], fontsize=2.5, barsize = 0.1, hjust= 1, angle = 62) +
  geom_cladelabel(node=subset(famf, Var1==toplot[38])$num, label=toplot[38], fontsize=2.5, barsize=0.1, hjust= 1, angle = 345) +
  geom_cladelabel(node=subset(famf, Var1==toplot[39])$num, label=toplot[39], fontsize=2.5, barsize = 0.1,  hjust= 1, angle = 336) +
  geom_cladelabel(node=subset(famf, Var1==toplot[40])$num, label=toplot[40], fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 300) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[41])$num, label=toplot[41],  fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[42])$num, label=toplot[42], fontsize=2.5, barsize = 0.1, hjust= 1, angle = 39) +
  geom_cladelabel(node=subset(famf, Var1==toplot[43])$num, label=toplot[43], fontsize=2.5, barsize=0.1, hjust= 1, angle = 298) +
  geom_cladelabel(node=subset(famf, Var1==toplot[44])$num, label=toplot[44], fontsize=2.5, barsize = 0.5,  hjust= 1, angle = 346) +
  geom_cladelabel(node=subset(famf, Var1==toplot[45])$num, label=toplot[45], fontsize=2.5, barsize = 0.5,  angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[46])$num, label=toplot[46],  fontsize=2.5, barsize = 0.5, hjust= 1, angle = 333) +
  geom_cladelabel(node=subset(famf, Var1==toplot[47])$num, label=toplot[47], fontsize=2.5, barsize = 0.5, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[48])$num, label=toplot[48], fontsize=2.5, barsize=0.5, hjust= 1, angle = 316) +
  geom_cladelabel(node=subset(famf, Var1==toplot[49])$num, label=toplot[49], fontsize=2.5, barsize = 0.5,  angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[50])$num, label=toplot[50], fontsize=2.5, barsize = 0.5,  angle = "auto") +
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="bottom")

#save output:
png("phylo_ring_all_mult.png",
    res=300,height=8,width=8,units="in"); 
p
dev.off()

#clean-up:
#rm(list = ls())