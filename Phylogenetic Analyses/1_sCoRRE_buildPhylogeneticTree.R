################################################################################
##  build_phylo_tree.R: Building a phylogenetic tree for the species in the CoRRE database.
##
##  Author: Josep Padulles Cubino, Kimberly Komatsu
##  Date created: December 19, 2021
################################################################################

#load libraries:
library(devtools)
install_github("jinyizju/U.PhyloMaker") #reference paper: https://doi.org/10.1016/j.pld.2022.12.007
library(ape)
library(U.PhyloMaker)
library(Taxonstand)
library(rlist)
library(tidyverse)

#### set directory ####
my.wd <- "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/" #padu
setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data\\')  #kim's computer


#### read data ####

#spp names
names <- read.csv('trait data\\corre2trykey_2021.csv') %>% 
  select(genus_species, species_matched) %>% 
  unique()

#community data
comm <- read.table("CoRRE data\\community composition\\CoRRE_RelativeCover_Jan2023.csv", header=T, sep=",", fill = TRUE) %>% 
  left_join(names) %>% 
  left_join(read.csv('trait data\\sCoRRE categorical trait data_12142022.csv')) %>% 
  filter(growth_form!='moss' & leaf_type!='microphyll') %>% 
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block, plot_id, family, genus_species, species_matched, relcov) %>% 
  filter(!grepl("(sp.)$", species_matched)) #remove species identified only at the genus level

#species list set up for phylogenies
spp <- comm %>% 
  select(family, species_matched) %>% 
  unique() %>% 
  rename(species=species_matched) %>% 
  mutate(genus=gsub("([A-Za-z]+).*", "\\1", species),
         species.relative=rep(NA, nrow(.)),
         genus.relative=rep(NA, nrow(.))) %>% 
  select(species, genus, family, species.relative, genus.relative)

#unify family names:
spp$family[spp$family=="Compositae"]<-"Asteraceae"
spp$family[spp$family=="Leguminosae"]<-"Fabaceae"
spp$family[spp$family=="Viburnaceae"]<-"Adoxaceae"
spp$family[spp$family=="Polypodiaceae"]<-"Dryopteridaceae"
spp$family[spp$genus=="Athyrium"]<-"Athyriaceae"
spp$family[spp$genus=="Blechnum"]<-"Blechnaceae"
spp$family[spp$genus=="Gymnocarpium"]<-"Cystopteridaceae"
spp$family[spp$genus=="Matteuccia"]<-"Onocleaceae"
spp$family[spp$genus=="Onoclea"]<-"Onocleaceae"
spp$family[spp$genus=="Phegopteris"]<-"Thelypteridaceae"
spp$family[spp$genus=="Thelypteris"]<-"Thelypteridaceae"
spp$family[spp$genus=="Bassia"]<-"Amaranthaceae"
spp$family[spp$genus=="Chenopodium"]<-"Amaranthaceae"
spp$family[spp$genus=="Corispermum"]<-"Amaranthaceae"
spp$family[spp$genus=="Carex"]<-"Cyperaceae"
spp$family[spp$genus=="Juncus"]<-"Juncaceae"

# spp$species_matched <- trimws(spp$species_matched, which="right") #delete empty spaces to the right of plant names


#### create phylogeny ####

#build phylo tree based on Scenario 3
megatree <- read.tree("Phylogenies\\plant_megatree.tre")

gen.list <- read.csv("Phylogenies\\plant_genus_list.csv")

# generate a phylogeny for the sample species list
result <- phylo.maker(spp, megatree, gen.list, nodes.type = 1, scenario = 3)

#save tree:
write.tree(result$phylo, "Phylogenies\\scorre.phylo.tree.S3_20230427.tre")

plot.phylo(result$phylo)


# # Question 2) Produce a X number of trees randomly placing tips on the phylogeny. Scenario 2.
# 
# #This operation might be time consuming. In my laptop, creating 1 tree takes about 76 seconds.
# #This means that 100 trees would take 2 hours. I'm setting it to 10 trees, but we'll have to change it accordingly.
# 
# #note that some taxa might not be included because their families are not included in the tree.
# scorre.trees<-list()
# for (i in 1:10) #replace 10 by 100 to produce 100 trees.
# {
#   print(i)
#   scorre.trees[[i]] <- phylo.maker(sp.list = spp1, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S2")
# }
# list.save(scorre.trees, paste(my.wd, 'scorre.tree.S2.rdata', sep=""))
