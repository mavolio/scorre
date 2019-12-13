
#####
# Script to build phylogenies for sCoRre
#####

#Note: The following script does not include species identified at the genus level.
#Including these species would be useful to calculate phylogenetic diversity.

#load libraries:
library(V.PhyloMaker) #reference paper: https://doi.org/10.1111/ecog.04434
library(Taxonstand)
library(rlist)

#read data:
spp<-read.table("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/CoRRE_TRY_species_list.csv", header=T, sep=",", fill = TRUE)
spp<-unique(spp$species_matched) #get unique names of species

#check species list against taxonstand to get families:
sppt<-TPL(spp) #from package Taxonstand
#write.table(sppt, "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/TPL_spp_scorre.csv") #save output
sppt <- read.table("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/TPL_spp_scorre.csv")

#rearrange dataset to create table in the right format for V.PhyloMaker:
sppt2<-sppt[c(1,13)]
names(sppt2)[1]<-paste("species")
names(sppt2)[2]<-paste("family")
sppt2$genus <- gsub("([A-Za-z]+).*", "\\1", sppt2$species)
sppt2$species.relative <- rep(NA, nrow(sppt2))
sppt2$genus.relative <- rep(NA, nrow(sppt2))
sppt2<-sppt2[c(1,3,2,4,5)]

#get list of species not included in the phylogeny: about 29%
spp.not.tree <- setdiff((gsub(" ", "_", sppt2$species)), tips.info$species)

#build phylo tree:
EVA.tree <- phylo.maker(sp.list = sppt2, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S3")

#ups! looks like some species have wrong family names. We fix it manually:
sppt2$family <- as.character(sppt2$family)
sppt2$family[sppt2$family == "Compositae"] <- "Asteraceae"
sppt2$family[sppt2$family == "Leguminosae"] <- "Fabaceae"
sppt2$family[sppt2$family == "Xanthorrhoeaceae"] <- "Asphodelaceae"
sppt2$family[sppt2$family == "Centrolepidaceae"] <- "Restionaceae"
sppt2[1888,3]<-"Montiaceae"
sppt2[1057,3]<-"Asphodelaceae"
sppt2[1414,3]<-"Mazaceae"

#create list of briophytes to exclude from original list of species:
non.vascular <-  c("Anthelia_juratzkana", "Cyrtomnium_hymenophyllum", "Distichium_capillaceum",
                   "Dicranum_brevifolium", "Dicranum_elongatum", "Dicranum_fuscescens", 
                   "Dicranum_groenlandicum", "Dicranum_scoparium", "Gymnomitrion_concinnatum",
                   "Hylocomium_splendens", "Kiaeria_starkei", "Leiocolea_heterocolpos",
                   "Marsupella_brevissima", "Aulacomnium_turgidum", "Pleurozium_schreberi",
                   "Pogonatum_urnigerum", "Pohlia_cruda", "Pohlia_nutans",
                   "Polytrichastrum_alpinum", "Polytrichum_juniperinum", "Polytrichum_piliferum",
                   "Polytrichum_strictum", "Ptilidium_ciliare", "Racomitrium_lanuginosum",
                   "Rhytidium_rugosum", "Tomentypnum_nitens", "Tritomaria_quinquedentata",
                   "Barbilophozia_hatcheri", "Barbilophozia_kunzeana", "Blepharostoma_trichophyllum",
                   "Bryum_pseudotriquetrum", "Hypnum_cupressiforme", "Andreaea_obovata",
                   "Hamatocaulis_vernicosus", "Hypnum_hamulosum", "Isopterygiopsis_pulchella",  
                   "Oncophorus_virens", "Schistidium_apocarpum", "Syntrichia_ruralis",         
                   "Bryum_arcticum", "Dicranoweisia_crispula", "Preissia_quadrata",         
                   "Marchantia_polymorpha", "Brachythecium_albicans", "Campylium_stellatum",        
                   "Ditrichum_flexicaule", "Meesia_uliginosa", "Myurella_tenerrima",         
                   "Oncophorus_wahlenbergii", "Saelania_glaucescens", "Sanionia_uncinata",          
                   "Tortella_tortuosa", "Homalothecium_pinnatifidum")
non.vascular<-gsub("_", " ", non.vascular)
sppt2 <- sppt2[!sppt2$species %in% non.vascular, ] #remove
#write.table(sppt2, "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/CoRRE_tax_rank.csv") #save output
#sppt2<-read.table("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/CoRRE_tax_rank.csv")

#get list of species not included in the phylogeny: about 28%
spp.not.tree <- setdiff((gsub(" ", "_", sppt2$species)), tips.info$species)

#build phylo tree based on Scenario 3. Check pag 4 in https://doi.org/10.1111/ecog.04434 for further details.
#this approach produces the same phylogenetic tree every time. This might be useful for Question 1, where we
#need to plot "winners" and "losers" to a fix phylogeny.
scorre.tree <- phylo.maker(sp.list = sppt2, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S3")

#save tree:
write.tree(scorre.tree$scenario.3, "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/scorre.tree.S3.tre")
#rm(scorre.tree, spp, sppt, sppt2)

#Alternatively, we can also produce a determined number of trees randomly placing tips on the phylogeny.
#This approach might be useful to calculate phylogenetic diversity. It corresponds to Scenario 2

#This operation might be time consuming. In my laptop, creating 1 tree takes about 76 seconds.
#This means that 1000 trees would take 21 hours. I'm setting it to 3 trees, but we'll have to change it accordingly.

scorre.trees<-list()
for (i in 1:3) #replace 3 by 1000 to produce 1000 trees.
{
  scorre.trees[[i]] <- phylo.maker(sp.list = sppt2, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S2")
}
list.save(scorre.trees, '/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/scorre.tree.S2.rdata')
scorre.trees<-list.load("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/scorre.tree.S2.rdata")
#rm(trees, sppt2, i)

#end of code