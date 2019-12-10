

######
# Script to calculate phylogenetic diversity/structure
#####

#load packages:
library(tidyverse)
library(reshape2)
library(PhyloMeasures)
library(V.PhyloMaker)

#read data:
comm<-read.table("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/CoRRE_relative_abundance_Nov2019.csv", header=T, sep=",", fill = TRUE)
spp<-read.table("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/CoRRE_TRY_species_list.csv", header=T, sep=",", fill = TRUE)

#reduce spp to original and new name:
spp<-unique(spp[c(1,2)])

#do some preliminary cleaning:
comm$genus_species <- trimws(comm$genus_species, which="right")

#merge with original dataset:
comm<-merge(comm, spp, by="genus_species", all.x=T)

#create new column with unique plot identifier:
comm <- comm %>% mutate(plot_id2 = paste(site_code, project_name, community_type,
                                         treatment_year, plot_id, sep = "_"))

#create list of sites:
sites <- unique(comm$site_code)

#Calculate Faith's PD within sites:

#we will need a loop here:
comm2<-subset(comm, site_code == sites[2]) #subset plot within each site
comm2<-comm2[c(13,12)] #reorder table
comm2 <- dcast(comm2, plot_id2 ~ species_matched) # apply dcast to turn into community matrix:
head(comm2)
head(comm2)
