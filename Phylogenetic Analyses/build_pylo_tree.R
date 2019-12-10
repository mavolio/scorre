
#load libraries:
library(V.PhyloMaker)
library(Taxonstand)

#read data:
spp<-read.table("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/CoRRE_TRY_species_list.csv", header=T, sep=",", fill = TRUE)
spp<-unique(spp$species_matched) #get unique names of species

#check species list against taxonstand to get families:
sppt<-TPL(spp)
#write.table(sppt, "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/TPL_spp_scorre.csv") #save output
sppt <- read.table("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/TPL_spp_scorre.csv")

#rearrange dataset:
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

#get list of species not included in the phylogeny: about 28%
spp.not.tree <- setdiff((gsub(" ", "_", sppt2$species)), tips.info$species)

#build phylo tree:
scorre.tree <- phylo.maker(sp.list = sppt2, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S3")

#save tree:
write.tree(scorre.tree$scenario.3,"/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/scorre.tree.tre")
#rm(scorre.tree, spp, sppt, sppt2)
