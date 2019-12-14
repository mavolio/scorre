
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
spp<-subset(spp, type != "moss/lichen") #filter out mosses and lichens
spp$species_matched <- trimws(spp$species_matched, which="right") #delete spaces on plant names

#get families for species:
fam<-read.table("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/species_families.csv", header=T, sep=";", fill = TRUE)
spp<-merge(spp, fam, by="species_matched", all.x=T)
spp$genus <- gsub("([A-Za-z]+).*", "\\1", spp$species) #column with genus
spp$species.relative <- rep(NA, nrow(spp))
spp$genus.relative <- rep(NA, nrow(spp))

#Some species that could not be linked to the phylogeny have to be removed.
#I quickly checked them some of them, and they were mosses. Will send the list
#to Kim to classify them accordingly in the original file. Meanwhile, just remove manually.
non.vascular <-  c("Andreaea_obovata",            "Anthelia_juratzkana" ,       "Aulacomnium_turgidum",       
                   "Barbilophozia_hatcheri",      "Barbilophozia_kunzeana" ,     "Blepharostoma_trichophyllum",
                   "Brachythecium_albicans",      "Bryum_arcticum"   ,           "Bryum_pseudotriquetrum",     
                   "Campylium_stellatum",         "Cyrtomnium_hymenophyllum" ,   "Dicranoweisia_crispula",     
                   "Dicranum_brevifolium",        "Dicranum_elongatum"  ,        "Dicranum_fuscescens",        
                   "Dicranum_groenlandicum",      "Dicranum_scoparium" ,         "Distichium_capillaceum",     
                   "Ditrichum_flexicaule",        "Gymnomitrion_concinnatum" ,   "Hamatocaulis_vernicosus",    
                   "Homalothecium_pinnatifidum",  "Hylocomium_splendens",        "Hypnum_cupressiforme",       
                   "Hypnum_hamulosum",            "Isopterygiopsis_pulchella",   "Kiaeria_starkei",            
                   "Leiocolea_heterocolpos",      "Marchantia_polymorpha",       "Marsupella_brevissima",      
                   "Meesia_uliginosa",            "Myurella_tenerrima",          "Oncophorus_virens",         
                   "Oncophorus_wahlenbergii",     "Pleurozium_schreberi",        "Pogonatum_urnigerum" ,       
                   "Pohlia_cruda" ,               "Pohlia_nutans",               "Polytrichastrum_alpinum",    
                   "Polytrichum_juniperinum",     "Polytrichum_piliferum",       "Polytrichum_strictum",       
                   "Preissia_quadrata",           "Ptilidium_ciliare",           "Racomitrium_lanuginosum",    
                   "Rhytidium_rugosum",           "Saelania_glaucescens",        "Sanionia_uncinata",          
                   "Schistidium_apocarpum",       "Syntrichia_ruralis",          "Tomentypnum_nitens",         
                   "Tortella_tortuosa",           "Tritomaria_quinquedentata", "Barbilophozia_sp.",
                   "Cynodontium_sp.",   "Dicranella_sp.",    "Encalypta_sp.",     "Hypnum_sp.",       
                   "Pohlia_sp.",        "Polytrichum_sp.",   "Racomitrium_sp.",   "Scapania_sp.",
                   "Syntrichia_sp.")
non.vascular<-gsub("_", " ", non.vascular)
spp <- spp[!spp$species_matched %in% non.vascular, ] #remove

#we create different phylogenies for questions 1 and 2.

# 1) Phylogeny for Q1 with Scenario 3
spp1<-subset(spp, type == "identified species") #remove species recorded at the genus level.
names(spp1)[1]<-paste("species")
spp1<-spp1[c(1,6,5,7,8)] #rearrange
spp1<-unique(spp1) #remove duplicated rows

#build phylo tree based on Scenario 3. Check pag 4 in https://doi.org/10.1111/ecog.04434 for further details.
#this approach produces the same phylogenetic tree every time. This might be useful for Question 1, where we
#need to plot "winners" and "losers" to a fix phylogeny.
EVA.tree <- phylo.maker(sp.list = spp1, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S3")

#save tree:
write.tree(scorre.tree$scenario.3, "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/scorre.tree.S3.tre")
#rm(scorre.tree, spp, sppt, sppt2)

# 2) Phylogeny for Q2.  Produce a X number of trees randomly placing tips on the phylogeny.

#This operation might be time consuming. In my laptop, creating 1 tree takes about 76 seconds.
#This means that 1000 trees would take 21 hours. I'm setting it to 3 trees, but we'll have to change it accordingly.
spp1<-spp
names(spp1)[1]<-paste("species")
spp1<-spp1[c(1,6,5,7,8)] #rearrange
spp1<-unique(spp1) #remove duplicated rows

#note that some taxa might not be included because their families are not included in the tree.
for (i in 1:3) #replace 3 by 1000 to produce 1000 trees.
{
  scorre.trees[[i]] <- phylo.maker(sp.list = spp1, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S2")
}
list.save(scorre.trees, '/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/scorre.tree.S2.rdata')
scorre.trees<-list.load("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/scorre.tree.S2.rdata")
#rm(trees, sppt2, i)

#end of code