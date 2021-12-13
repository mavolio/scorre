
######
# Script to calculate phylogenetic diversity/structure
#####

#load packages:
library(tidyverse)
library(reshape2)
library(PhyloMeasures) #reference paper: https://onlinelibrary.wiley.com/doi/10.1111/ecog.01814
library(V.PhyloMaker) #reference paper: https://doi.org/10.1111/ecog.04434
library(rlist)
library(matrixStats)
library(plyr)
library(magrittr)

#set directory:
my.wd <- "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"

#read data:
comm<-read.table(paste(my.wd, "CoRRE_relative_abundance_Feb2021.csv", sep=""), header=T, sep=",", fill = TRUE)
spp<-read.table(paste(my.wd, "CoRRE_TRY_species_list.csv", sep=""), header=T, sep=",", fill = TRUE)

#reduce spp to original and new name:
spp<-subset(spp, type != "moss/lichen") #filter out mosses and lichens

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
                   "Tortella tortuosa",           "Tritomaria quinquedentata", "Barbilophozia sp.",
                   "Cynodontium sp.",   "Dicranella sp.",    "Encalypta sp.",     "Hypnum sp.",       
                   "Pohlia sp.",        "Polytrichum sp.",   "Racomitrium sp.",   "Scapania sp.",
                   "Syntrichia sp.", "Nephroma arcticum", "Unknown NA")
spp <- spp[!spp$species_matched %in% non.vascular, ] #remove
spp<-unique(spp[c(2,3)]) #get unique list of species

#do some preliminary cleaning to remove empty spaces on names:
comm$genus_species <- trimws(comm$genus_species, which="right")
spp$species_matched <- trimws(spp$species_matched, which="right")

#merge with original dataset:
comm<-merge(comm, spp, by="genus_species", all.x=T)
comm<-comm[complete.cases(comm), ] #by now, we remove cases with NA to get rid of false species (i.e, unknowns,
#taxa recorded at family, level, etc.). However, we'll have to check this is 
#OK and we're not deleting accepted taxa.

#create new column with unique plot identifier:
comm <- comm %>% mutate(plot_id2 = paste(site_code, project_name, community_type,
                                         treatment_year, plot_id, sep = "::"))

#remove sites with only 1 species:
comm <- comm[comm$plot_id2 %in% comm$plot_id2[duplicated(comm$plot_id2) | duplicated(comm$plot_id2, fromLast=TRUE)], ]

#create list of sites:
sites <- unique(comm$site_code)




####
#Example to calculate phylogenetic diversity metrics using one single tree (scenario 3) (non-weighted by abundances):

#load trees:
scorre.tree<-read.tree(paste(my.wd, "scorre.phylo.tree.S3.tre", sep=""))

pd.all<-NULL
for (i in 1:length(sites)) #loop to calculate metrics for each site independently
{
  print(i*100/length(sites))
  comm2<-subset(comm, site_code == sites[i]) #subset plot within each site
  comm2<-comm2[c(16,15,14)] #reorder table
  comm2<-comm2[!duplicated(comm2[c(1,2)]),] #remove duplicated rows
  comm2<-subset(comm2, relcov>0) #remove species with zero cover
  comm2 <- dcast(comm2, plot_id2 ~ species_matched, value.var="relcov") # apply dcast to turn into community matrix
  rownames(comm2)<-comm2$plot_id2 #assign rownames
  comm2$plot_id2 <-NULL #and detele original columns
  comm2[is.na(comm2)]<-0 #turn NAs into 0s
  colnames(comm2)<-gsub(" ", "_", colnames(comm2))
  comm2[comm2>0]<-1 #replace all positive values with 1
  
  #Prune tree with only species in our site:
  tree<-keep.tip(scorre.tree, colnames(comm2))
  
  #calculate Faith's diversity (PD):
  pd.raw<-as.data.frame(pd.query(tree, comm2,  standardize = F))
  pd.ses<-as.data.frame(pd.query(tree, comm2,  null.model="uniform", reps=1000, standardize = T))
  pdp<-as.data.frame(1-(pd.pvalues(tree, comm2,  null.model="uniform", reps=1000)))
  pd<-cbind(pd.raw, pd.ses, pdp)
  
  #calculate Mean Pairwise Distances (MPD):
  mpd.raw<-as.data.frame(mpd.query(tree, comm2,  standardize = F))
  mpd.ses<-as.data.frame(mpd.query(tree, comm2,  null.model="uniform", reps=1000, standardize = T))
  mpdp<-as.data.frame(1-(mpd.pvalues(tree, comm2,  null.model="uniform", reps=1000)))
  mpd<-cbind(mpd.raw, mpd.ses, mpdp)
  
  #calculate Mean Nearest Taxon Distances (MNTD):
  mntd.raw<-as.data.frame(mntd.query(tree, comm2,  standardize = F))
  mntd.ses<-as.data.frame(mntd.query(tree, comm2,  null.model="uniform", reps=1000, standardize = T))
  mntdp<-as.data.frame(1-(mntd.pvalues(tree, comm2,  null.model="uniform", reps=1000)))
  mntd<-cbind(mntd.raw, mntd.ses, mntdp)
  
  #merge:
  pd<-cbind(pd, mpd, mntd)
  rownames(pd)<-rownames(comm2)
  colnames(pd) <- c("pd.raw", "pd.ses", "pd.pval", "mpd.raw", "mpd.ses", "mpd.pval", "mntd.raw", "mntd.ses", "mntd.pval")
  pd.all<-rbind(pd.all, pd)
}

#save output:
write.table(pd.all, paste(my.wd, "CoRRE_pd_metrics_non_weighted.csv", sep=""))





####
#Example to calculate phylogenetic diversity metrics using one single tree (scenario 3) (weighted by abundances):

#load trees:
scorre.tree<-read.tree(paste(my.wd, "scorre.phylo.tree.S3.tre", sep=""))

pd.all<-NULL
for (i in 1:length(sites)) #loop to calculate metrics for each site independently
{
  print(i*100/length(sites))
  comm2<-subset(comm, site_code == sites[i]) #subset plot within each site
  comm2<-comm2[c(16,15,14)] #reorder table
  comm2<-comm2[!duplicated(comm2[c(1,2)]),] #remove duplicated rows
  comm2<-subset(comm2, relcov>0) #remove species with zero cover
  comm2 <- dcast(comm2, plot_id2 ~ species_matched, value.var="relcov") # apply dcast to turn into community matrix
  rownames(comm2)<-comm2$plot_id2 #assign rownames
  comm2$plot_id2 <-NULL #and detele original columns
  comm2[is.na(comm2)]<-0 #turn NAs into 0s
  colnames(comm2)<-gsub(" ", "_", colnames(comm2))
  
  #create vector for abundances:
  weights <- colSums(comm2)
  comm2[comm2>0]<-1 #replace all positive values with 1

  #Prune tree with only species in our site:
  tree<-keep.tip(scorre.tree, colnames(comm2))
  
  #calculate Faith's diversity (PD):
  pd.raw<-as.data.frame(pd.query(tree, comm2,  abundance.weights=weights, standardize = F))
  pd.ses<-as.data.frame(pd.query(tree, comm2,  null.model="sequential", abundance.weights=weights, reps=1000, standardize = T))
  pdp<-as.data.frame(1-(pd.pvalues(tree, comm2,  null.model="sequential", abundance.weights=weights, reps=1000)))
  pd<-cbind(pd.raw, pd.ses, pdp)

  #calculate Mean Pairwise Distances (MPD):
  mpd.raw<-as.data.frame(mpd.query(tree, comm2,  abundance.weights=weights, standardize = F))
  mpd.ses<-as.data.frame(mpd.query(tree, comm2,  null.model="sequential", abundance.weights=weights, reps=1000, standardize = T))
  mpdp<-as.data.frame(1-(mpd.pvalues(tree, comm2,  null.model="sequential", abundance.weights=weights, reps=1000)))
  mpd<-cbind(mpd.raw, mpd.ses, mpdp)
  
  #calculate Mean Nearest Taxon Distances (MNTD):
  mntd.raw<-as.data.frame(mntd.query(tree, comm2,  abundance.weights=weights, standardize = F))
  mntd.ses<-as.data.frame(mntd.query(tree, comm2,  null.model="sequential", abundance.weights=weights, reps=1000, standardize = T))
  mntdp<-as.data.frame(1-(mntd.pvalues(tree, comm2,  null.model="sequential", abundance.weights=weights, reps=1000)))
  mntd<-cbind(mntd.raw, mntd.ses, mntdp)
  
  #merge:
  pd<-cbind(pd, mpd, mntd)
  rownames(pd)<-rownames(comm2)
  colnames(pd) <- c("pd.raw", "pd.ses", "pd.pval", "mpd.raw", "mpd.ses", "mpd.pval", "mntd.raw", "mntd.ses", "mntd.pval")
  pd.all<-rbind(pd.all, pd)
}

#save output:
write.table(pd.all, paste(my.wd, "CoRRE_pd_metrics_weighted.csv", sep=","))

#clean-up:
#rm(list = ls())
