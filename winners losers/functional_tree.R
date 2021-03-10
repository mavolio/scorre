
######
# Script to produce trait dendrogram
######

#load packages:
library(gawdis)
library(data.table)
library(magrittr)

my.wd<-"/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"

#load continuous traits & clean-up:
cont<-fread(paste(my.wd, "TRY_sCorre.csv", sep=""), dec=",")
cont %<>% dplyr::select( -X, -ObservationID, -Genus, -Family) %>% unique
cont <- aggregate(cont[, 3:34], list(cont$Species), mean, na.rm=T) #aggregate by species
names(cont)[1]<-paste("species_matched")

#load categorical traits & clean-up:
cate<-fread(paste(my.wd, "sCoRRE categorical trait data - traits_complete_pre spot check_03082021.csv", sep=""), sep=",")
cate <- cate[-c(2386:2389), ] #remove empty rows
cate %<>% dplyr::select( -spotcheck_assignment, -done, -assigned, -family, -leaf_type_source, -leaf_compoundness_source,
                         -growth_form_source, -photosynthetic_pathway_source, -lifespan_source, -stem_support_source,
                         -clonal_source, -pollination_source, -dispersal_mode_source, -mycorrhizal_source, 
                         -mycorrhizal_type_source, -n_fixation_source, -rhizobial_source, -actinorhizal_source,
                         -NOTES, -pollination, -dispersal_mode, -mycorrhizal, -mycorrhizal_type, -n_fixation,
                         -rhizobial, -actinorhizal)
cate$`Double-Checked By`<-NULL
cate[cate=="CHECK"]<-NA #In the categorical trait matrix there are a lot of CHECKS. I'll change them to NAs.

#merge both categorical and numerical:
trait<-as.data.frame(merge(cate, cont, by="species_matched", all.x=T))
rownames(trait)<-trait$species_matched
trait$species_matched<-NULL
trait <- trait[complete.cases(trait), ] # get only complete cases

#calculate functional dissimilarities between species using GAWDIS:
trait.dis<-gawdis(trait) #can take a while

#convert distance object into phylotree using the Ward's algorithm (https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust):
ftree<-as.phylo(hclust(trait.dis, method="ward.D2"))

#save output:
write.tree(ftree, paste(my.wd, "ftree.scorre.tre", sep=""))

#clean-up:
rm(list = ls())
