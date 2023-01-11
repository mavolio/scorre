
####
# Script to build final tree and highlight the nodes
####

# Load libraries:
library(phytools)
library(plyr)
library(ggtree)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(ggplotify)
library(grid)

#set directory.
my.wd<-"/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"
#my.wd<-"E:\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\"
#my.wd<-"C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\"

# Load the tree:
tree<-read.tree(paste(my.wd, "scorre.phylo.tree.S3.tre", sep="")) #why not using the original tree? I detected that some species are missing from the other.

# Load list of species with their families:
fam<-read.table(paste(my.wd, "species_families_trees_2021.csv", sep=""), header=T, sep=",", fill = TRUE)
fam$species_matched<-gsub(" ", "_", fam$species_matched) #unify nomenclature

# Get list of taxa that are not in the tree but they are in the table (only two):
in_tree_not_fam<-setdiff(tree$tip.label, fam$species_matched)
tree<-drop.tip(tree, in_tree_not_fam)

###
#get number of nodes per family:
list.r<-c("Poaceae", "Brassicaceae", "Solanaceae", "Cyperaceae", "Polemoniaceae",
        "Gentianaceae", "Plantaginaceae", "Euphorbiaceae", "Amaranthaceae",
        "Orchidaceae", "Fabaceae", "Gentianaceae", "Orobanchaceae", "Lamiaceae")
  
#get nodes for families:
list.nod1<-NULL
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$species_matched) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- findMRCA(tree, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod1<-rbind(list.nod1, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod1)[1]<-paste("Var1")


###
# Assign all families to nodes
###

#load family data and clean-up:
fam<-read.table(paste(my.wd, "species_families_trees_2021.csv",sep=""), header=T, sep=",", fill = TRUE)[,-c(3)] #load data
fam$species_matched<-gsub(" ", "_", fam$species_matched) #adapt species nomenclature
fam$family[fam$family=="Compositae"]<-"Asteraceae" #replace family name
fam$family[fam$family=="Leguminosae"]<-"Fabaceae" #replace family name
fam$family[fam$family=="Viburnaceae"]<-"Adoxaceae" #replace family name
fam<-fam[which(fam$species_matched %in% tree$tip.label),] #subset only species included in our treatment

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
    num <- findMRCA(tree, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)

###
# Plot phylogenetic tree with families
###

#get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=70)) #select the top 9 families with 5 for more species

#get groups by nodes:
tree2 <- groupClade(tree, c(list.nod1$num, 2597, 2381, 2725)) #the two extra nodes are from Asteraceae 

# Plot tree:
p <- 
  ggtree(tree2, aes(color=group), layout="circular", size=0.5)+ # build circular tree
  scale_color_manual(values=c("grey80", rep("black", 16)))+
  #geom_text(aes(label=node), size=1, col="red")+
  
  geom_cladelabel(node=subset(famf, Var1==toplot[1])$num, label=toplot[1], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[2])$num, label=toplot[2], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[3])$num, label=toplot[3], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, hjust=1, angle = 7) +
  geom_cladelabel(node=subset(famf, Var1==toplot[4])$num, label=toplot[4], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[5])$num, label=toplot[5], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 349) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[6])$num, label=toplot[6], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, hjust=1, angle = 313) +
  geom_cladelabel(node=subset(famf, Var1==toplot[7])$num, label=toplot[7], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[8])$num, label=toplot[8], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 285) +
  geom_cladelabel(node=subset(famf, Var1==toplot[9])$num, label=toplot[9], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 42) +
  geom_cladelabel(node=subset(famf, Var1==toplot[10])$num, label=toplot[10], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[11])$num, label=toplot[11], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 24) +
  geom_cladelabel(node=subset(famf, Var1==toplot[12])$num, label=toplot[12], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, hjust=1, angle = 85) +
  geom_cladelabel(node=subset(famf, Var1==toplot[13])$num, label=toplot[13], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 300) +
  geom_cladelabel(node=subset(famf, Var1==toplot[14])$num, label=toplot[14], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, hjust=1, angle = 325) +
  geom_cladelabel(node=subset(famf, Var1==toplot[15])$num, label=toplot[15], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, hjust=1, angle = 35) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[16])$num, label=toplot[16], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 78) +
  geom_cladelabel(node=subset(famf, Var1==toplot[17])$num, label=toplot[17], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, hjust=1, angle = 73) +
  geom_cladelabel(node=subset(famf, Var1==toplot[18])$num, label=toplot[18], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[19])$num, label=toplot[19], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 50) +
  geom_cladelabel(node=subset(famf, Var1==toplot[20])$num, label=toplot[20], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, hjust=1, angle = 62) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[21])$num, label=toplot[21], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 66) +
  geom_cladelabel(node=subset(famf, Var1==toplot[22])$num, label=toplot[22], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 334) +
  geom_cladelabel(node=subset(famf, Var1==toplot[23])$num, label=toplot[23], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[24])$num, label=toplot[24], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, hjust=1, angle = 56) +
  geom_cladelabel(node=subset(famf, Var1==toplot[25])$num, label=toplot[25], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 295) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[26])$num, label=toplot[26], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 69) +
  geom_cladelabel(node=subset(famf, Var1==toplot[27])$num, label=toplot[27], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 308) +
  geom_cladelabel(node=subset(famf, Var1==toplot[28])$num, label=toplot[28], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[29])$num, label=toplot[29], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[30])$num, label=toplot[30], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[31])$num, label=toplot[31], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 59) +
  geom_cladelabel(node=subset(famf, Var1==toplot[32])$num, label=toplot[32], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 54) +
  geom_cladelabel(node=subset(famf, Var1==toplot[33])$num, label=toplot[33], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 337) +
  geom_cladelabel(node=subset(famf, Var1==toplot[34])$num, label=toplot[34], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 32) +
  geom_cladelabel(node=subset(famf, Var1==toplot[35])$num, label=toplot[35], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[36])$num, label=toplot[36], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 331) +
  geom_cladelabel(node=subset(famf, Var1==toplot[37])$num, label=toplot[37], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 291) +
  geom_cladelabel(node=subset(famf, Var1==toplot[38])$num, label=toplot[38], offset=1, fontsize=2.8, fontface="bold", barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[39])$num, label=toplot[39], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 334) +
  geom_cladelabel(node=subset(famf, Var1==toplot[40])$num, label=toplot[40], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[41])$num, label=toplot[41], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[42])$num, label=toplot[42], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 28) +
  geom_cladelabel(node=subset(famf, Var1==toplot[43])$num, label=toplot[43], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 321) +
  geom_cladelabel(node=subset(famf, Var1==toplot[44])$num, label=toplot[44], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[45])$num, label=toplot[45], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[46])$num, label=toplot[46], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[47])$num, label=toplot[47], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 46) +
  geom_cladelabel(node=subset(famf, Var1==toplot[48])$num, label=toplot[48], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[49])$num, label=toplot[49], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 304) +
  geom_cladelabel(node=subset(famf, Var1==toplot[50])$num, label=toplot[50], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 290) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[51])$num, label=toplot[51], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[52])$num, label=toplot[52], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  #geom_cladelabel(node=subset(famf, Var1==toplot[53])$num, label=toplot[53], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[54])$num, label=toplot[54], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[55])$num, label=toplot[55], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 329) +
  
  #geom_cladelabel(node=subset(famf, Var1==toplot[56])$num, label=toplot[56], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[57])$num, label=toplot[57], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 355) +
  geom_cladelabel(node=subset(famf, Var1==toplot[58])$num, label=toplot[58], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 342) +
  geom_cladelabel(node=subset(famf, Var1==toplot[59])$num, label=toplot[59], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  #geom_cladelabel(node=subset(famf, Var1==toplot[60])$num, label=toplot[60], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[61])$num, label=toplot[61], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 306) +
  geom_cladelabel(node=subset(famf, Var1==toplot[62])$num, label=toplot[62], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[63])$num, label=toplot[63], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 20) +
  geom_cladelabel(node=subset(famf, Var1==toplot[64])$num, label=toplot[64], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 339) +
  geom_cladelabel(node=subset(famf, Var1==toplot[65])$num, label=toplot[65], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 29) +
  
  #geom_cladelabel(node=subset(famf, Var1==toplot[66])$num, label=toplot[66], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[67])$num, label=toplot[67], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 293) +
  geom_cladelabel(node=subset(famf, Var1==toplot[68])$num, label=toplot[68], offset=1, fontsize=2.8, barsize = 0.2, hjust=1, angle = 319) +
  geom_cladelabel(node=subset(famf, Var1==toplot[69])$num, label=toplot[69], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  #geom_cladelabel(node=subset(famf, Var1==toplot[70])$num, label=toplot[70], offset=1, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="none")



####
#run pie2_function first

trt<-c("CO2", "Drt.","Irg.", "Temp.", "N", "P", " Mult.","Inter.")
trt2<-c(expression("CO"[2]), "Drt.","Irg.", "Temp.", "N", "P", " Mult.","Inter.")
value<-c(1,1,1,1,1,1,1,1)
pie<-data.frame(trt, value)

##doing this in baseR
collst = RColorBrewer::brewer.pal(n = 8, name = "Dark2")

#make Legend
lg<-as.grob(~pie2(pie$value, labels=trt2, col=collst, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for orchidaceae
labellistorc<-c("  ", "  ","  ","  ", "  ","  ","-","  ")
colactive<-collst
colactive[labellistorc==""]="darkgray"
colactive[labellistorc=="  "]=NA
orc<-as.grob(~pie2(pie$value, labels=labellistorc, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for cyperaceae
labellistcyp<-c("", "","","", "+","","+","")
colactive<-collst
colactive[labellistcyp==""]="darkgray"
colactive[labellistcyp=="  "]=NA
cyp<-as.grob(~pie2(pie$value, labels=labellistcyp, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for poaceae overall
labellistpoa<-c("", "","","", "","","+","+")
colactive<-collst
colactive[labellistpoa==""]="darkgray"
colactive[labellistpoa=="  "]=NA
pot<-as.grob(~pie2(pie$value, labels=labellistpoa, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for poaceae sub1
labellistpoa<-c("  ", "","","", "","+","","")
colactive<-collst
colactive[labellistpoa==""]="darkgray"
colactive[labellistpoa=="  "]=NA
po1<-as.grob(~pie2(pie$value, labels=labellistpoa, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for poaceae sub2
labellistpoa<-c("", "-","","", "","","","")
colactive<-collst
colactive[labellistpoa==""]="darkgray"
colactive[labellistpoa=="  "]=NA
po2<-as.grob(~pie2(pie$value, labels=labellistpoa, col=colactive, border="white",line_length = 1, text_center = 0.8, textcol = "white"))

#for Brassicaceae
labellistbra<-c("", "","","", "","","+","")
colactive<-collst
colactive[labellistbra==""]="darkgray"
colactive[labellistbra=="  "]=NA
bra<-as.grob(~pie2(pie$value, labels=labellistbra, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for Euphorbiaceae
labellisteup<-c("  ", "","","+", "","","","")
colactive<-collst
colactive[labellisteup==""]="darkgray"
colactive[labellisteup=="  "]=NA
eup<-as.grob(~pie2(pie$value, labels=labellisteup, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for Fabaceae
labellistfab<-c("", "-","+","", "-","","-","")
colactive<-collst
colactive[labellistfab==""]="darkgray"
colactive[labellistfab=="  "]=NA
fab<-as.grob(~pie2(pie$value, labels=labellistfab, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for section Fabaceae
labellistfab<-c("", "-","+", "", "", "", "","")
colactive<-collst
colactive[labellistfab==""]="darkgray"
colactive[labellistfab=="  "]=NA
fab1<-as.grob(~pie2(pie$value, labels=labellistfab, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for asteraceae sub1
labellistast<-c("  ", "  ","","", "-","","","")
colactive<-collst
colactive[labellistast==""]="darkgray"
colactive[labellistast=="  "]=NA
ast1<-as.grob(~pie2(pie$value, labels=labellistast, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for asteraceae sub2
labellistast<-c("  ", "","+","", "","","","")
colactive<-collst
colactive[labellistast==""]="darkgray"
colactive[labellistast=="  "]=NA
ast2<-as.grob(~pie2(pie$value, labels=labellistast, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for Lamiacea and Orobanaceae and Polimoniaceae
labellistLO<-c("  ", "","","", "-","","","")
colactive<-collst
colactive[labellistLO==""]="darkgray"
colactive[labellistLO=="  "]=NA
lop<-as.grob(~pie2(pie$value, labels=labellistLO, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for Plantaginaceae
labellistpla<-c("", "","","", "","","","-")
colactive<-collst
colactive[labellistpla==""]="darkgray"
colactive[labellistpla=="  "]=NA
plt<-as.grob(~pie2(pie$value, labels=labellistpla, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for Gentinaceae
labellistgen<-c("  ", "  ","","", "-","","-","-")
colactive<-collst
colactive[labellistgen==""]="darkgray"
colactive[labellistgen=="  "]=NA
gen<-as.grob(~pie2(pie$value, labels=labellistgen, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for Solaneaceae
labellistsol<-c("  ", "","","", "+","","+","")
colactive<-collst
colactive[labellistsol==""]="darkgray"
colactive[labellistsol=="  "]=NA
sol<-as.grob(~pie2(pie$value, labels=labellistsol, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))

#for Amaranth
labellistama<-c("", "+","-","", "+","-","","")
colactive<-collst
colactive[labellistama==""]="darkgray"
colactive[labellistama=="  "]=NA
ama<-as.grob(~pie2(pie$value, labels=labellistama, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white"))


# Create layout to print:
h <- ggdraw(p)

# Add pie charts to the tree:
p2<- h + draw_grob(lg, x=.38, y=.5, width=.35, height=.35, hjust=0.5, vjust=0.5)+ #central one
         draw_grob(orc, x=.76, y=.57, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(cyp, x=.756, y=.655, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(pot, x=.70, y=.71, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(po2, x=.59, y=.80, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(po1, x=.49, y=.835, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(bra, x=.31, y=.71, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(eup, x=.26, y=.65, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(fab, x=.24, y=.50, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(fab1, x=.20, y=.44, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(ama, x=.27, y=.315, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(lop, x=.31, y=.26, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(sol, x=.372, y=.235, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(gen, x=.428, y=.197, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(plt, x=.494, y=.204, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(lop, x=.554, y=.177, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(ast2, x=.72, y=.25, width=.12, height=.12, hjust=0.5, vjust=0.5)+
         draw_grob(ast1, x=.765, y=.30, width=.12, height=.12, hjust=0.5, vjust=0.5)

#Save output:
png("phylo_ring_main.png", res=600,height=8,width=8,units="in"); 
p2
grid.text("Tribes Trifolieae,\nFabeae & Galegeae", x = unit(0.12, "npc"), y = unit(0.4, "npc"), rot=13, hjust=1, gp=gpar(fontsize=6.5, col="gray23"))
grid.text("Tribe Cardueae", x = unit(0.787, "npc"), y = unit(0.125, "npc"), rot=305, hjust=1, gp=gpar(fontsize=6.5, col="gray23"))
grid.text("Tribe Gnaphalieae", x = unit(0.862, "npc"), y = unit(0.25, "npc"), rot=328,  gp=gpar(fontsize=6.5, col="gray23"))
grid.text("Tribe Poeae", x = unit(0.49, "npc"), y = unit(0.925, "npc"), rot=90, gp=gpar(fontsize=6.5, col="gray23"))
grid.text("Subfamily Pooidae", x = unit(0.62, "npc"), y = unit(0.925, "npc"), rot=75, gp=gpar(fontsize=6.5, col="gray23"))
dev.off()

#clean-up:
rm(list = ls())

