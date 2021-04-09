
######
# Script to produce trait dendrogram
######

#load packages:
library(gawdis)
library(data.table)
library(magrittr)
library(dplyr)
library(plyr)

#set directory.
my.wd<-"/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"

#load data.
species.data<-read.table(paste(my.wd,"Species_DCiDiff_newtrts_filtered.csv",sep=""), header=T, sep=" ")

#load continuous traits & clean-up:
cont<-fread(paste(my.wd, "TRY_new.csv", sep=""), dec=",")
cont %<>% dplyr::select( -X, -ObservationID, -Genus, -Family) %>% unique #remove unnecessary columns
cont[,c(3:35)] %<>% mutate_if(is.character,as.numeric) #convert columns from character to numeric
cont <- aggregate(cont[, 3:35], list(cont$Species), mean, na.rm=T) #get mean values by species
names(cont)[1]<-paste("species_matched") #change column name
cont$species_matched<-gsub(" ", "_", cont$species_matched) #unify nomenclature

#Subset only species in sCoRRe and the traits we are interested in:
cont %<>% dplyr::select(species_matched, X11, X47, X3112, X55, X26, X3120, X14, X163, X146, X15) %>% unique #subset traits
#X11=SLA; X47=LDMC; X3112= LEAF AREA; X55=LEAF DRY MASS; X3120=LEAF WATER CONTENT; X14=LEAF N DRY MASS; X163=LEAF FRESH MASS;
#X146=LEAF C/N; X15: LEAF P.
cont <- cont[cont$species_matched %in% unique(species.data$species_matched), ] #subset species from sCoRRe

#load categorical traits & clean-up:
cate<-fread(paste(my.wd, "sCoRRE categorical trait data - traits_complete_pre spot check_03082021.csv", sep=""), sep=",")
cate <- cate[-c(2386:2389), ] #remove empty rows
cate %<>% dplyr::select( -spotcheck_assignment, -done, -assigned, -family, -leaf_type_source, -leaf_compoundness_source,
                         -growth_form_source, -photosynthetic_pathway_source, -lifespan_source, -stem_support_source,
                         -clonal_source, -pollination_source, -dispersal_mode_source, -mycorrhizal_source, 
                         -mycorrhizal_type_source, -n_fixation_source, -rhizobial_source, -actinorhizal_source,
                         -NOTES, -pollination, -dispersal_mode, -mycorrhizal, -mycorrhizal_type, -n_fixation,
                         -rhizobial, -actinorhizal, -leaf_type, -leaf_compoundness, -growth_form, -stem_support)
cate$`Double-Checked By`<-NULL
cate[cate=="CHECK"]<-NA #In the categorical trait matrix there are a lot of CHECKS. I'll change them to NAs.
cate$species_matched<-gsub(" ", "_", cate$species_matched) #unify nomenclature

#subset sCoRRe species and redefine categories:
cate <- cate[cate$species_matched %in% unique(species.data$species_matched), ] #subset species from sCoRRe

#Photosynthetic pathway:
cate$photosynthetic_pathway<-as.factor(cate$photosynthetic_pathway)
cate$photosynthetic_pathway<-revalue(cate$photosynthetic_pathway, c("C3"="1", "C4"="0", "CAM"="0"))
names(cate)[2]<-paste("C3")

#Lifespan:
cate$lifespan<-as.factor(cate$lifespan)
cate$lifespan<-revalue(cate$lifespan, c("annual"="1", "biennial"="0", "perennial"="0"))
names(cate)[3]<-paste("annual")

#Clonality:
cate$clonal<-as.factor(cate$clonal)
cate$clonal<-revalue(cate$clonal, c("no"="0", "yes"="1"))

#change binary variable to numeric:
cate[,c(2:4)] %<>% mutate_if(is.factor, as.character)
cate[,c(2:4)] %<>% mutate_if(is.character, as.numeric)

#merge both categorical and numerical:
trait<-as.data.frame(merge(cont, cate, by="species_matched", all.x=T))
rownames(trait)<-trait$species_matched
trait$species_matched<-NULL
trait <- trait[complete.cases(trait), ] # get only complete cases

#log transform continous variables:
vars <- colnames(trait[,c(1:10)])
trait[vars] <- lapply(trait[vars], log)

#save for later:
#write.table(trait, paste(my.wd, "trait_data_filtered.csv", sep=""))

#calculate functional dissimilarities between species using GAWDIS:
trait.dis<-gawdis(trait) #can take a while

#or calculate functional dissimilarities between species using Gower distances:
trait.dis<-cluster::daisy(trait, metric="gower") #can take a while

#convert distance object into phylotree using the Ward's algorithm (https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust):
ftree<-as.phylo(hclust(trait.dis, method="average"))

#save output:
write.tree(ftree, paste(my.wd, "ftree.scorre.gowdis.log.upgma.tre", sep=""))

#clean-up:
rm(list = ls())
