
####
# Script to build final tree and highlight the nodes
####

# Load libraries:
library(phytools)

#set directory.
my.wd<-"/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"
my.wd<-"E:\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\"
my.wd<-"C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\"

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
faml<-c("Poaceae", "Asteraceae", "Brassicaceae", "Solanaceae", "Cyperaceae",
        "Gentianaceae", "Plantaginaceae", "Euphorbiaceae", "Amaranthaceae",
        "Orchidaceae", "Fabaceae", "Gentianaceae", "Orobanchaceae", "Lamiaceae", 
        "Polimonaceae")
  
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

