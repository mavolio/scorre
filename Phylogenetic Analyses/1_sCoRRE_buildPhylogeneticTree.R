################################################################################
##  build_phylo_tree.R: Building a phylogenetic tree for the species in the CoRRE database.
##
##  Author: Josep Padulles Cubino
##  Date created: December 19, 2021
################################################################################

#load libraries:
library(V.PhyloMaker) #reference paper: https://doi.org/10.1111/ecog.04434
library(Taxonstand)
library(rlist)
library(plyr)

#set directory:
my.wd <- "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"

#read data:
spp<-read.table(paste(my.wd, "FullList_Nov2021.csv", sep=""), header=T, sep=",", fill = TRUE)
spp<-subset(spp, type != "moss/lichen") #filter out mosses and lichens
spp$species_matched <- trimws(spp$species_matched, which="right") #delete empty spaces to the right of plant names

#manually correct one mistake:
spp$species_matched <- revalue(spp$species_matched, c("Aronia x"="Aronia x prunifolia"))

#Some species that could not be matched to the phylogeny still need to be removed.
#They are all briophytes. I will send the list to Kim. Meanwhile, I remove them manually.
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
                   "Tortella tortuosa",           "Tritomaria quinquedentata", "Nephroma arcticum" ,"Unknown NA",
                   "Campylopus flexuosus",       "Hypnum jutlandicum",         "Plagiothecium undulatum",    
                   "Polytrichum commune",  "Pseudoscleropodium purum",   "Rhytidiadelphus loreus",
                   "Rhytidiadelphus triquetrus", "Thuidium tamariscinum")
spp<-subset(spp, !grepl("(sp.)$", species_matched)) #remove species identified at the genus level
non.vascular<-gsub("_", " ", non.vascular)
spp <- spp[!spp$species_matched %in% non.vascular, ] #remove

#get families for species and rearrange to create necessary fields for the phylogenies:
fam<-read.table(paste(my.wd, "species_families_trees_2021.csv", sep=""), header=T, sep=",", fill = TRUE)
fam$species_matched <- revalue(fam$species_matched, c("Aronia x"="Aronia x prunifolia"))
spp<-merge(spp, fam, by="species_matched", all.x=T)
spp$genus <- gsub("([A-Za-z]+).*", "\\1", spp$species) #column with genus
spp$species.relative <- rep(NA, nrow(spp))
spp$genus.relative <- rep(NA, nrow(spp))


#########
#We create different phylogenies for questions 1 and 2.

# Question 1) Using Scenario 3.
spp1<-spp
#spp1<-subset(spp, type == "identified species") #remove species recorded at the genus level.
names(spp1)[1]<-paste("species")
spp1<-spp1[c(1,8,6,9,10)] #rearrange
spp1<-unique(spp1) #remove duplicated rows

#get list of species not included in the tree:
spp.not.tree <- setdiff((gsub(" ", "_", spp1$species)), tips.info$species)

#unify family names:
spp1$family[spp1$family=="Compositae"]<-"Asteraceae"
spp1$family[spp1$family=="Leguminosae"]<-"Fabaceae"
spp1$family[spp1$family=="Viburnaceae"]<-"Adoxaceae"
spp1$family[spp1$family=="Polypodiaceae"]<-"Dryopteridaceae"
spp1$family[spp1$genus=="Athyrium"]<-"Athyriaceae"
spp1$family[spp1$genus=="Blechnum"]<-"Blechnaceae"
spp1$family[spp1$genus=="Gymnocarpium"]<-"Cystopteridaceae"
spp1$family[spp1$genus=="Matteuccia"]<-"Onocleaceae"
spp1$family[spp1$genus=="Onoclea"]<-"Onocleaceae"
spp1$family[spp1$genus=="Phegopteris"]<-"Thelypteridaceae"
spp1$family[spp1$genus=="Thelypteris"]<-"Thelypteridaceae"

#build phylo tree based on Scenario 3. Check pag 4 in https://doi.org/10.1111/ecog.04434 for further details.
#this approach produces the same phylogenetic tree every time. This might be useful for Question 1, where we
#need to plot "winners" and "losers" to a fix phylogeny.
scorre.tree <- phylo.maker(sp.list = spp1, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S3")

#save tree:
write.tree(scorre.tree$scenario.3, paste(my.wd, "scorre.phylo.tree.S3.tre", sep=""))
#rm(scorre.tree, spp, sppt, sppt2)

# Question 2) Produce a X number of trees randomly placing tips on the phylogeny. Scenario 2.

#This operation might be time consuming. In my laptop, creating 1 tree takes about 76 seconds.
#This means that 100 trees would take 2 hours. I'm setting it to 10 trees, but we'll have to change it accordingly.

#note that some taxa might not be included because their families are not included in the tree.
scorre.trees<-list()
for (i in 1:10) #replace 10 by 100 to produce 100 trees.
{
  print(i)
  scorre.trees[[i]] <- phylo.maker(sp.list = spp1, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S2")
}
list.save(scorre.trees, paste(my.wd, 'scorre.tree.S2.rdata', sep=""))

#clean-up:
#rm(list = ls())
