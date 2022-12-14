
#load libraries:
library(devtools)
#devtools::install_github("GuangchuangYu/ggtree")

library(ggtree)
library(ggplot2)
library(stringr)
library(plyr)
library(ape)
library(scico)
library(phytools)
library(grid)

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

#res<-node.mean(tree2, dat, 999)
#write.table(res, paste(my.wd, "res_phylo_all_mult_Nov22.csv", sep="")) #save the result
res<-read.table(paste(my.wd, "res_phylo_all_mult_Nov22.csv", sep=""))
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

#To colour branches based on mean values:
svl <- as.matrix(dat)[,1]
fit <- phytools::fastAnc(tree2, svl, vars=TRUE, CI=TRUE)
td <- data.frame(node = nodeid(tree2, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
treef <- dplyr::full_join(tree2, d, by = 'node')

###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=50)) #select the top 9 families with 5 for more species

# Plot tree:
p <- 
  ggtree(treef, layout="circular", size=0.5)+ # build circular tree
  geom_tree(aes(color=trait), continuous = 'colour', show.legend = F) +
  scale_color_scico(palette = "vik", direction=-1, na.value = "gray48", limits = c(-1, 1)) +
  
  ggnewscale::new_scale("color") +
  
  geom_point(aes(color=as.factor(significant)), size=2, alpha=1, show.legend = F) + # highlight nodes
  scale_colour_manual(values=c("#8A6000", "#006FA4"), labels=c("Loss", "Gain"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  
  geom_cladelabel(node=subset(famf, Var1==toplot[1])$num, label=toplot[1], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[2])$num, label=toplot[2], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[3])$num, label=toplot[3], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[4])$num, label=toplot[4], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[5])$num, label=toplot[5], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[6])$num, label=toplot[6], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[7])$num, label=toplot[7], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[8])$num, label=toplot[8], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[9])$num, label=toplot[9], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[10])$num, label=toplot[10], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[11])$num, label=toplot[11], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[12])$num, label=toplot[12], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[13])$num, label=toplot[13], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[14])$num, label=toplot[14], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[15])$num, label=toplot[15], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[16])$num, label=toplot[16], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[17])$num, label=toplot[17], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[18])$num, label=toplot[18], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[19])$num, label=toplot[19], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[20])$num, label=toplot[20], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[21])$num, label=toplot[21], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[22])$num, label=toplot[22], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[23])$num, label=toplot[23], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[24])$num, label=toplot[24], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[25])$num, label=toplot[25], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[26])$num, label=toplot[26], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[27])$num, label=toplot[27], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[28])$num, label=toplot[28], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[29])$num, label=toplot[29], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[30])$num, label=toplot[30], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[31])$num, label=toplot[31], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[32])$num, label=toplot[32], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[33])$num, label=toplot[33], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[34])$num, label=toplot[34], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[35])$num, label=toplot[35], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[36])$num, label=toplot[36], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[37])$num, label=toplot[37], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[38])$num, label=toplot[38], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[39])$num, label=toplot[39], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[40])$num, label=toplot[40], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[41])$num, label=toplot[41], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[42])$num, label=toplot[42], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[43])$num, label=toplot[43], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[44])$num, label=toplot[44], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  #geom_cladelabel(node=subset(famf, Var1==toplot[45])$num, label=toplot[45], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  #geom_cladelabel(node=subset(famf, Var1==toplot[46])$num, label=toplot[46], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[47])$num, label=toplot[47], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[48])$num, label=toplot[48], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[49])$num, label=toplot[49], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[50])$num, label=toplot[50], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

# Add heatmap with tendency values:
p1 <- gheatmap(p, dat, offset=0.3, width=.03, colnames = F, color = NULL) +
  scale_fill_scico(palette = "vik", direction=-1, na.value = "gray48",
                   limits = c(-0.6, 0.6)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=15),
        legend.title=element_text(colour="white"))


#save output:
#png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_all_mult_withheat.png",
png("phylo_ring_all_mult_withheat.png", res=300,height=8,width=8,units="in"); 
p1
grid.text("Winners", x = unit(0.725, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#006FA4"))
grid.text("Losers", x = unit(0.35, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#8A6000"))

grid.text("(High DCi)", x = unit(0.725, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
grid.text("(Low DCi)", x = unit(0.35, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))

dev.off()

#clean-up:
#rm(list = ls())


###
# Filter by treatment == all nuts
###

dat<-subset(species.data, trt_type2=="all nuts")[,c(1,2)] #select "all nuts" treatment from original data
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

#res<-node.mean(tree2, dat, 999)
#write.table(res, paste(my.wd, "res_phylo_all_nuts_Nov22.csv", sep="")) #save the result
res<-read.table(paste(my.wd, "res_phylo_all_nuts_Nov22.csv", sep=""))
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


#To colour branches based on mean values:
svl <- as.matrix(dat)[,1]
fit <- phytools::fastAnc(tree2, svl, vars=TRUE, CI=TRUE)
td <- data.frame(node = nodeid(tree2, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
treef <- dplyr::full_join(tree2, d, by = 'node')

###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=50)) #select the top 9 families with 5 for more species

# Plot tree:
p <- 
  ggtree(treef, layout="circular", size=0.5)+ # build circular tree
  geom_tree(aes(color=trait), continuous = 'colour', show.legend = F) +
  scale_color_scico(palette = "vik", direction=-1, na.value = "gray48", limits = c(-1, 1)) +
  
  ggnewscale::new_scale("color") +
  
  geom_point(aes(color=as.factor(significant)), size=2, alpha=1, show.legend = F) + # highlight nodes
  scale_colour_manual(values=c("#8A6000", "#006FA4"), labels=c("Loss", "Gain"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  
  geom_cladelabel(node=subset(famf, Var1==toplot[1])$num, label=toplot[1], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[2])$num, label=toplot[2], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[3])$num, label=toplot[3], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[4])$num, label=toplot[4], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[5])$num, label=toplot[5], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[6])$num, label=toplot[6], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[7])$num, label=toplot[7], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[8])$num, label=toplot[8], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[9])$num, label=toplot[9], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[10])$num, label=toplot[10], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[11])$num, label=toplot[11], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[12])$num, label=toplot[12], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[13])$num, label=toplot[13], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[14])$num, label=toplot[14], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[15])$num, label=toplot[15], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[16])$num, label=toplot[16], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[17])$num, label=toplot[17], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[18])$num, label=toplot[18], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[19])$num, label=toplot[19], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[20])$num, label=toplot[20], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[21])$num, label=toplot[21], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[22])$num, label=toplot[22], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[23])$num, label=toplot[23], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[24])$num, label=toplot[24], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[25])$num, label=toplot[25], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[26])$num, label=toplot[26], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[27])$num, label=toplot[27], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[28])$num, label=toplot[28], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[29])$num, label=toplot[29], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[30])$num, label=toplot[30], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[31])$num, label=toplot[31], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[32])$num, label=toplot[32], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[33])$num, label=toplot[33], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[34])$num, label=toplot[34], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[35])$num, label=toplot[35], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[36])$num, label=toplot[36], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[37])$num, label=toplot[37], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[38])$num, label=toplot[38], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[39])$num, label=toplot[39], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[40])$num, label=toplot[40], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[41])$num, label=toplot[41], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[42])$num, label=toplot[42], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[43])$num, label=toplot[43], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[44])$num, label=toplot[44], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[45])$num, label=toplot[45], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[46])$num, label=toplot[46], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[47])$num, label=toplot[47], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[48])$num, label=toplot[48], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[49])$num, label=toplot[49], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[50])$num, label=toplot[50], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

# Add heatmap with tendency values:
p1 <- gheatmap(p, dat, offset=0.3, width=.03, colnames = F, color = NULL) +
  scale_fill_scico(palette = "vik", direction=-1, na.value = "gray48",
                   limits = c(-0.6, 0.6)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=15),
        legend.title=element_text(colour="white"))

#save output:
#png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_all_nutsNov22.png", 
#    res=300,height=8,width=8,units="in"); 
png("phylo_ring_all_nutsNov22.png", res=300,height=8,width=8,units="in"); 
p1
grid.text("Winners", x = unit(0.725, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#006FA4"))
grid.text("Losers", x = unit(0.35, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#8A6000"))

grid.text("(High DCi)", x = unit(0.725, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
grid.text("(Low DCi)", x = unit(0.35, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
dev.off()


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

#res<-node.mean(tree2, dat, 999)
#write.table(res, paste(my.wd, "res_phylo_n.csv", sep="")) #save the result
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


#To colour branches based on mean values:
svl <- as.matrix(dat)[,1]
fit <- phytools::fastAnc(tree2, svl, vars=TRUE, CI=TRUE)
td <- data.frame(node = nodeid(tree2, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
treef <- dplyr::full_join(tree2, d, by = 'node')

###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=50)) #select the top 9 families with 5 for more species

# Plot tree:
p <- 
  ggtree(treef, layout="circular", size=0.5)+ # build circular tree
  geom_tree(aes(color=trait), continuous = 'colour', show.legend = F) +
  scale_color_scico(palette = "vik", direction=-1, na.value = "gray48", limits = c(-1, 1)) +
  
  ggnewscale::new_scale("color") +
  
  geom_point(aes(color=as.factor(significant)), size=2, alpha=1, show.legend = F) + # highlight nodes
  scale_colour_manual(values=c("#8A6000", "#006FA4"), labels=c("Loss", "Gain"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  
  geom_cladelabel(node=subset(famf, Var1==toplot[1])$num, label=toplot[1], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[2])$num, label=toplot[2], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[3])$num, label=toplot[3], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[4])$num, label=toplot[4], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[5])$num, label=toplot[5], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[6])$num, label=toplot[6], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[7])$num, label=toplot[7], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[8])$num, label=toplot[8], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[9])$num, label=toplot[9], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[10])$num, label=toplot[10], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[11])$num, label=toplot[11], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[12])$num, label=toplot[12], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[13])$num, label=toplot[13], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[14])$num, label=toplot[14], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[15])$num, label=toplot[15], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[16])$num, label=toplot[16], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[17])$num, label=toplot[17], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[18])$num, label=toplot[18], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[19])$num, label=toplot[19], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[20])$num, label=toplot[20], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[21])$num, label=toplot[21], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[22])$num, label=toplot[22], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[23])$num, label=toplot[23], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[24])$num, label=toplot[24], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[25])$num, label=toplot[25], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[26])$num, label=toplot[26], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[27])$num, label=toplot[27], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[28])$num, label=toplot[28], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[29])$num, label=toplot[29], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[30])$num, label=toplot[30], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[31])$num, label=toplot[31], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[32])$num, label=toplot[32], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[33])$num, label=toplot[33], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[34])$num, label=toplot[34], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[35])$num, label=toplot[35], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[36])$num, label=toplot[36], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[37])$num, label=toplot[37], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  #geom_cladelabel(node=subset(famf, Var1==toplot[38])$num, label=toplot[38], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[39])$num, label=toplot[39], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[40])$num, label=toplot[40], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[41])$num, label=toplot[41], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[42])$num, label=toplot[42], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[43])$num, label=toplot[43], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[44])$num, label=toplot[44], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[45])$num, label=toplot[45], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[46])$num, label=toplot[46], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[47])$num, label=toplot[47], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  #geom_cladelabel(node=subset(famf, Var1==toplot[48])$num, label=toplot[48], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[49])$num, label=toplot[49], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[50])$num, label=toplot[50], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

# Add heatmap with tendency values:
p1 <- gheatmap(p, dat, offset=0.3, width=.03, colnames = F, color = NULL) +
  scale_fill_scico(palette = "vik", direction=-1, na.value = "gray48",
                   limits = c(-0.6, 0.6)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=15),
        legend.title=element_text(colour="white"))


#save output:
#png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_n_withheat.png",
#    res=300,height=8,width=8,units="in"); 
png("phylo_ring_n_withheat.png", res=300,height=8,width=8,units="in"); 
p1
grid.text("Winners", x = unit(0.725, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#006FA4"))
grid.text("Losers", x = unit(0.35, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#8A6000"))

grid.text("(High DCi)", x = unit(0.725, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
grid.text("(Low DCi)", x = unit(0.35, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
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

#res<-node.mean(tree2, dat, 999)
#write.table(res, paste(my.wd, "res_phylo_p.csv", sep="")) #save the result
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


#To colour branches based on mean values:
svl <- as.matrix(dat)[,1]
fit <- phytools::fastAnc(tree2, svl, vars=TRUE, CI=TRUE)
td <- data.frame(node = nodeid(tree2, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
treef <- dplyr::full_join(tree2, d, by = 'node')

###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=50)) #select the top 9 families with 5 for more species

# Plot tree:
p <- 
  ggtree(treef, layout="circular", size=0.5)+ # build circular tree
  geom_tree(aes(color=trait), continuous = 'colour', show.legend = F) +
  scale_color_scico(palette = "vik", direction=-1, na.value = "gray48", limits = c(-1, 1)) +

  ggnewscale::new_scale("color") +
  
  geom_point(aes(color=as.factor(significant)), size=2, alpha=1, show.legend = F) + # highlight nodes
  scale_colour_manual(values=c("#8A6000", "#006FA4"), labels=c("Loss", "Gain"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  #geom_text(aes(label=node), size=1)+
  
  geom_cladelabel(node=subset(famf, Var1==toplot[1])$num, label=toplot[1], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[2])$num, label=toplot[2], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[3])$num, label=toplot[3], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[4])$num, label=toplot[4], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[5])$num, label=toplot[5], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[6])$num, label=toplot[6], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[7])$num, label=toplot[7], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[8])$num, label=toplot[8], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[9])$num, label=toplot[9], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[10])$num, label=toplot[10], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[11])$num, label=toplot[11], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[12])$num, label=toplot[12], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[13])$num, label=toplot[13], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[14])$num, label=toplot[14], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[15])$num, label=toplot[15], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[16])$num, label=toplot[16], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[17])$num, label=toplot[17], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[18])$num, label=toplot[18], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[19])$num, label=toplot[19], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[20])$num, label=toplot[20], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[21])$num, label=toplot[21], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[22])$num, label=toplot[22], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[23])$num, label=toplot[23], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[24])$num, label=toplot[24], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[25])$num, label=toplot[25], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[26])$num, label=toplot[26], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[27])$num, label=toplot[27], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[28])$num, label=toplot[28], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[29])$num, label=toplot[29], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[30])$num, label=toplot[30], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[31])$num, label=toplot[31], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[32])$num, label=toplot[32], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[33])$num, label=toplot[33], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[34])$num, label=toplot[34], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[35])$num, label=toplot[35], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[36])$num, label=toplot[36], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[37])$num, label=toplot[37], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[38])$num, label=toplot[38], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[39])$num, label=toplot[39], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[40])$num, label=toplot[40], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[41])$num, label=toplot[41], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[42])$num, label=toplot[42], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[43])$num, label=toplot[43], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[44])$num, label=toplot[44], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[45])$num, label=toplot[45], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[46])$num, label=toplot[46], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[47])$num, label=toplot[47], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  #geom_cladelabel(node=subset(famf, Var1==toplot[48])$num, label=toplot[48], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[49])$num, label=toplot[49], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[50])$num, label=toplot[50], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

# Add heatmap with tendency values:
p1 <- gheatmap(p, dat, offset=0.3, width=.03, colnames = F, color = NULL) +
  scale_fill_scico(palette = "vik", direction=-1, na.value = "gray48",
                   limits = c(-0.6, 0.6)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=15),
        legend.title=element_text(colour="white"))


#save output:
#png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_p.png",
#    res=300,height=8,width=8,units="in"); 
png("phylo_ring_p_withheat.png", res=300,height=8,width=8,units="in"); 
p1
grid.text("Winners", x = unit(0.725, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#006FA4"))
grid.text("Losers", x = unit(0.35, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#8A6000"))

grid.text("(High DCi)", x = unit(0.725, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
grid.text("(Low DCi)", x = unit(0.35, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
dev.off()

#Get list of species for groups within families:

#Section Poaceae (increase) <-- node 1272 (unlock line in ggtree to see the name of the node)
caper::clade.members(1272, tree2, tip.labels = T, include.nodes=FALSE)


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

#res<-node.mean(tree2, dat, 999)
#write.table(res, paste(my.wd, "res_phylo_co2.csv", sep="")) #save the result
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

#To colour branches based on mean values:
svl <- as.matrix(dat)[,1]
fit <- phytools::fastAnc(tree2, svl, vars=TRUE, CI=TRUE)
td <- data.frame(node = nodeid(tree2, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
treef <- dplyr::full_join(tree2, d, by = 'node')

###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=50)) #select the top 9 families with 5 for more species

# Plot tree:
p <- 
  ggtree(treef, layout="circular", size=0.5)+ # build circular tree
  geom_tree(aes(color=trait), continuous = 'colour', show.legend = F) +
  scale_color_scico(palette = "vik", direction=-1, na.value = "gray48", limits = c(-1, 1)) +
  
  ggnewscale::new_scale("color") +
  
  geom_point(aes(color=as.factor(significant)), size=2, alpha=1, show.legend = F) + # highlight nodes
  scale_colour_manual(values=c("#8A6000", "#006FA4"), labels=c("Loss", "Gain"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  
  geom_cladelabel(node=subset(famf, Var1==toplot[1])$num, label=toplot[1], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[2])$num, label=toplot[2], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[3])$num, label=toplot[3], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[4])$num, label=toplot[4], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[5])$num, label=toplot[5], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[6])$num, label=toplot[6], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[7])$num, label=toplot[7], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[8])$num, label=toplot[8], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[9])$num, label=toplot[9], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[10])$num, label=toplot[10], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[11])$num, label=toplot[11], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[12])$num, label=toplot[12], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[13])$num, label=toplot[13], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[14])$num, label=toplot[14], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[15])$num, label=toplot[15], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[16])$num, label=toplot[16], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[17])$num, label=toplot[17], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[18])$num, label=toplot[18], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[19])$num, label=toplot[19], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[20])$num, label=toplot[20], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[21])$num, label=toplot[21], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[22])$num, label=toplot[22], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[23])$num, label=toplot[23], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[24])$num, label=toplot[24], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[25])$num, label=toplot[25], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[26])$num, label=toplot[26], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[27])$num, label=toplot[27], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[28])$num, label=toplot[28], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[29])$num, label=toplot[29], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[30])$num, label=toplot[30], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[31])$num, label=toplot[31], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[32])$num, label=toplot[32], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[33])$num, label=toplot[33], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[34])$num, label=toplot[34], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[35])$num, label=toplot[35], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[36])$num, label=toplot[36], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[37])$num, label=toplot[37], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[38])$num, label=toplot[38], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[39])$num, label=toplot[39], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[40])$num, label=toplot[40], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[41])$num, label=toplot[41], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[42])$num, label=toplot[42], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[43])$num, label=toplot[43], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[44])$num, label=toplot[44], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[45])$num, label=toplot[45], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[46])$num, label=toplot[46], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[47])$num, label=toplot[47], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[48])$num, label=toplot[48], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[49])$num, label=toplot[49], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[50])$num, label=toplot[50], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

# Add heatmap with tendency values:
p1 <- gheatmap(p, dat, offset=0.3, width=.03, colnames = F, color = NULL) +
  scale_fill_scico(palette = "vik", direction=-1, na.value = "gray48",
                   limits = c(-0.6, 0.6)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=15),
        legend.title=element_text(colour="white"))

#save output:
#png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_co2withheat.png",
#    res=300,height=8,width=8,units="in"); 
png("phylo_ring_co2_withheat.png", res=300,height=8,width=8,units="in"); 
p1
grid.text("Winners", x = unit(0.725, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#006FA4"))
grid.text("Losers", x = unit(0.35, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#8A6000"))

grid.text("(High DCi)", x = unit(0.725, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
grid.text("(Low DCi)", x = unit(0.35, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
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

#res<-node.mean(tree2, dat, 999)
#write.table(res, paste(my.wd, "res_phylo_drought.csv", sep="")) #save the result
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


#To colour branches based on mean values:
svl <- as.matrix(dat)[,1]
fit <- phytools::fastAnc(tree2, svl, vars=TRUE, CI=TRUE)
td <- data.frame(node = nodeid(tree2, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
treef <- dplyr::full_join(tree2, d, by = 'node')

###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=50)) #select the top 9 families with 5 for more species

# Plot tree:
p <- 
  ggtree(treef, layout="circular", size=0.5)+ # build circular tree
  geom_tree(aes(color=trait), continuous = 'colour', show.legend = F) +
  scale_color_scico(palette = "vik", direction=-1, na.value = "gray48", limits = c(-1, 1)) +
  
  ggnewscale::new_scale("color") +
  
  geom_point(aes(color=as.factor(significant)), size=2, alpha=1, show.legend = F) + # highlight nodes
  scale_colour_manual(values=c("#8A6000", "#006FA4"), labels=c("Loss", "Gain"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  #geom_text(aes(label=node), size=1)+
  
  geom_cladelabel(node=subset(famf, Var1==toplot[1])$num, label=toplot[1], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[2])$num, label=toplot[2], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[3])$num, label=toplot[3], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[4])$num, label=toplot[4], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[5])$num, label=toplot[5], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[6])$num, label=toplot[6], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[7])$num, label=toplot[7], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[8])$num, label=toplot[8], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[9])$num, label=toplot[9], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[10])$num, label=toplot[10], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[11])$num, label=toplot[11], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[12])$num, label=toplot[12], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[13])$num, label=toplot[13], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[14])$num, label=toplot[14], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[15])$num, label=toplot[15], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[16])$num, label=toplot[16], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[17])$num, label=toplot[17], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[18])$num, label=toplot[18], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[19])$num, label=toplot[19], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[20])$num, label=toplot[20], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[21])$num, label=toplot[21], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[22])$num, label=toplot[22], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[23])$num, label=toplot[23], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[24])$num, label=toplot[24], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[25])$num, label=toplot[25], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[26])$num, label=toplot[26], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[27])$num, label=toplot[27], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[28])$num, label=toplot[28], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[29])$num, label=toplot[29], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[30])$num, label=toplot[30], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[31])$num, label=toplot[31], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[32])$num, label=toplot[32], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[33])$num, label=toplot[33], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[34])$num, label=toplot[34], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[35])$num, label=toplot[35], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[36])$num, label=toplot[36], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[37])$num, label=toplot[37], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[38])$num, label=toplot[38], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[39])$num, label=toplot[39], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[40])$num, label=toplot[40], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[41])$num, label=toplot[41], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[42])$num, label=toplot[42], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[43])$num, label=toplot[43], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[44])$num, label=toplot[44], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[45])$num, label=toplot[45], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[46])$num, label=toplot[46], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[47])$num, label=toplot[47], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[48])$num, label=toplot[48], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[49])$num, label=toplot[49], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[50])$num, label=toplot[50], offset=6, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

# Add heatmap with tendency values:
p1 <- gheatmap(p, dat, offset=0.3, width=.03, colnames = F, color = NULL) +
  scale_fill_scico(palette = "vik", direction=-1, na.value = "gray48",
                   limits = c(-0.6, 0.6)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=15),
        legend.title=element_text(colour="white"))

#save output:
#png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_drought.png",
#    res=300,height=8,width=8,units="in"); 
png("phylo_ring_drought.png",  res=300,height=8,width=8,units="in"); 
p1
grid.text("Winners", x = unit(0.725, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#006FA4"))
grid.text("Losers", x = unit(0.35, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#8A6000"))

grid.text("(High DCi)", x = unit(0.725, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
grid.text("(Low DCi)", x = unit(0.35, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
dev.off()

#Section Fabaceae (decrease) <-- node 775 (unlock line in ggtree to see the name of the node)
caper::clade.members(775, tree2, tip.labels = T, include.nodes=FALSE)

#Section Poaceae (decrease) <-- node 971 (unlock line in ggtree to see the name of the node)
caper::clade.members(971, tree2, tip.labels = T, include.nodes=FALSE)




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

#res<-node.mean(tree2, dat, 999)
#write.table(res, paste(my.wd, "res_phylo_irrigation.csv", sep=""))
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


#To colour branches based on mean values:
svl <- as.matrix(dat)[,1]
fit <- phytools::fastAnc(tree2, svl, vars=TRUE, CI=TRUE)
td <- data.frame(node = nodeid(tree2, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
treef <- dplyr::full_join(tree2, d, by = 'node')

###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=50)) #select the top 9 families with 5 for more species

# Plot tree:
p <- 
  ggtree(treef, layout="circular", size=0.5)+ # build circular tree
  geom_tree(aes(color=trait), continuous = 'colour', show.legend = F) +
  scale_color_scico(palette = "vik", direction=-1, na.value = "gray48", limits = c(-1, 1)) +
  
  ggnewscale::new_scale("color") +
  
  geom_point(aes(color=as.factor(significant)), size=2, alpha=1, show.legend = F) + # highlight nodes
  scale_colour_manual(values=c("#8A6000", "#006FA4"), labels=c("Loss", "Gain"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  #geom_text(aes(label=node), size=1)+
  
  geom_cladelabel(node=subset(famf, Var1==toplot[1])$num, label=toplot[1], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[2])$num, label=toplot[2], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[3])$num, label=toplot[3], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[4])$num, label=toplot[4], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[5])$num, label=toplot[5], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[6])$num, label=toplot[6], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[7])$num, label=toplot[7], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[8])$num, label=toplot[8], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[9])$num, label=toplot[9], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[10])$num, label=toplot[10], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[11])$num, label=toplot[11], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[12])$num, label=toplot[12], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[13])$num, label=toplot[13], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[14])$num, label=toplot[14], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[15])$num, label=toplot[15], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[16])$num, label=toplot[16], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[17])$num, label=toplot[17], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[18])$num, label=toplot[18], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[19])$num, label=toplot[19], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[20])$num, label=toplot[20], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[21])$num, label=toplot[21], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[22])$num, label=toplot[22], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[23])$num, label=toplot[23], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[24])$num, label=toplot[24], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[25])$num, label=toplot[25], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[26])$num, label=toplot[26], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[27])$num, label=toplot[27], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[28])$num, label=toplot[28], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[29])$num, label=toplot[29], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[30])$num, label=toplot[30], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[31])$num, label=toplot[31], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[32])$num, label=toplot[32], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[33])$num, label=toplot[33], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[34])$num, label=toplot[34], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[35])$num, label=toplot[35], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[36])$num, label=toplot[36], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[37])$num, label=toplot[37], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[38])$num, label=toplot[38], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[39])$num, label=toplot[39], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[40])$num, label=toplot[40], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[41])$num, label=toplot[41], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[42])$num, label=toplot[42], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[43])$num, label=toplot[43], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[44])$num, label=toplot[44], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[45])$num, label=toplot[45], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[46])$num, label=toplot[46], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[47])$num, label=toplot[47], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[48])$num, label=toplot[48], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[49])$num, label=toplot[49], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[50])$num, label=toplot[50], offset=12,  fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

# Add heatmap with tendency values:
p1 <- gheatmap(p, dat, offset=0.3, width=.03, colnames = F, color = NULL) +
  scale_fill_scico(palette = "vik", direction=-1, na.value = "gray48",
                   limits = c(-0.6, 0.6)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=15),
        legend.title=element_text(colour="white"))


#save output:
#png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_irrigation_withheat.png",
#    res=300,height=8,width=8,units="in"); 
png("phylo_ring_irrigation_withheat.png", res=300,height=8,width=8,units="in"); 
p1
grid.text("Winners", x = unit(0.725, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#006FA4"))
grid.text("Losers", x = unit(0.35, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#8A6000"))

grid.text("(High DCi)", x = unit(0.725, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
grid.text("(Low DCi)", x = unit(0.35, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
dev.off()

#Get list of species for groups within families:

#Section Fabaceae (increase) <-- node 978 (unlock line in ggtree to see the name of the node)
caper::clade.members(978, tree2, tip.labels = T, include.nodes=FALSE)

#Section Asteraceae (increase) <-- node 797 (unlock line in ggtree to see the name of the node)
caper::clade.members(797, tree2, tip.labels = T, include.nodes=FALSE)




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

#res<-node.mean(tree2, dat, 999)
#write.table(res, paste(my.wd, "res_phylo_temp.csv", sep=""))
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


#To colour branches based on mean values:
svl <- as.matrix(dat)[,1]
fit <- phytools::fastAnc(tree2, svl, vars=TRUE, CI=TRUE)
td <- data.frame(node = nodeid(tree2, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
treef <- dplyr::full_join(tree2, d, by = 'node')

###
# Plot phylogenetic tree hihglighting nodes and families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=50)) #select the top 9 families with 5 for more species

# Plot tree:
p <- 
  ggtree(treef, layout="circular", size=0.5)+ # build circular tree
  geom_tree(aes(color=trait), continuous = 'colour', show.legend = F) +
  scale_color_scico(palette = "vik", direction=-1, na.value = "gray48", limits = c(-1, 1)) +
  
  ggnewscale::new_scale("color") +
  
  geom_point(aes(color=as.factor(significant)), size=2, alpha=1, show.legend = F) + # highlight nodes
  scale_colour_manual(values=c("#8A6000", "#006FA4"), labels=c("Loss", "Gain"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  
  geom_cladelabel(node=subset(famf, Var1==toplot[1])$num, label=toplot[1], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[2])$num, label=toplot[2], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[3])$num, label=toplot[3], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[4])$num, label=toplot[4], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[5])$num, label=toplot[5], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[6])$num, label=toplot[6], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[7])$num, label=toplot[7], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[8])$num, label=toplot[8], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[9])$num, label=toplot[9], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[10])$num, label=toplot[10], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[11])$num, label=toplot[11], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[12])$num, label=toplot[12], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[13])$num, label=toplot[13], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[14])$num, label=toplot[14], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[15])$num, label=toplot[15], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[16])$num, label=toplot[16], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[17])$num, label=toplot[17], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[18])$num, label=toplot[18], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[19])$num, label=toplot[19], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[20])$num, label=toplot[20], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[21])$num, label=toplot[21], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[22])$num, label=toplot[22], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[23])$num, label=toplot[23], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[24])$num, label=toplot[24], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[25])$num, label=toplot[25], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[26])$num, label=toplot[26], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[27])$num, label=toplot[27], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[28])$num, label=toplot[28], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[29])$num, label=toplot[29], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[30])$num, label=toplot[30], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[31])$num, label=toplot[31], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[32])$num, label=toplot[32], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[33])$num, label=toplot[33], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[34])$num, label=toplot[34], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[35])$num, label=toplot[35], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[36])$num, label=toplot[36], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[37])$num, label=toplot[37], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[38])$num, label=toplot[38], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[39])$num, label=toplot[39], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[40])$num, label=toplot[40], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[41])$num, label=toplot[41], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[42])$num, label=toplot[42], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[43])$num, label=toplot[43], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[44])$num, label=toplot[44], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[45])$num, label=toplot[45], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[46])$num, label=toplot[46], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[47])$num, label=toplot[47], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[48])$num, label=toplot[48], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[49])$num, label=toplot[49], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[50])$num, label=toplot[50], offset=12, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

# Add heatmap with tendency values:
p1 <- gheatmap(p, dat, offset=0.3, width=.03, colnames = F, color = NULL) +
  scale_fill_scico(palette = "vik", direction=-1, na.value = "gray48",
                   limits = c(-0.6, 0.6)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=15),
        legend.title=element_text(colour="white"))

#save output:
#png("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\Figs Dec 2021\\phylo_ring_temp.png",
#    res=300,height=8,width=8,units="in"); 
png("phylo_ring_temp.png", res=300,height=8,width=8,units="in"); 
p1
grid.text("Winners", x = unit(0.725, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#006FA4"))
grid.text("Losers", x = unit(0.35, "npc"), y = unit(0.09, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#8A6000"))

grid.text("(High DCi)", x = unit(0.725, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
grid.text("(Low DCi)", x = unit(0.35, "npc"), y = unit(0.06, "npc"), gp=gpar(fontsize=10))
dev.off()

#clean-up:
rm(list = ls())
