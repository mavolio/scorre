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
# library(Taxonstand)
library(rlist)
library(tidyverse)


#### read data ####

spp <- readRDS("PhyDiv_FuncDiv/sppList.rds") %>% 
  mutate(genus=gsub("([A-Za-z]+).*", "\\1", species),
         species.relative=rep(NA, nrow(.)),
         genus.relative=rep(NA, nrow(.))) %>% 
  select(species, genus, family, species.relative, genus.relative)

#unify family names:
spp$family[spp$family=="Polypodiaceae"]<-"Dryopteridaceae"


#### create phylogeny ####

#build phylo tree based on Scenario 3
megatree <- read.tree('https://raw.githubusercontent.com/megatrees/plant_20221117/refs/heads/main/plant_megatree.tre')

gen.list <- read.csv('https://raw.githubusercontent.com/megatrees/plant_20221214/refs/heads/main/plant_genus_list.csv')

# generate a phylogeny for the sample species list
result <- phylo.maker(spp, megatree, gen.list, nodes.type = 1, scenario = 3)

#save tree:
# write.tree(result$phylo, "PhyDiv_FuncDiv/corre.phylo.tree.S3_20250521.tre")

plot.phylo(result$phylo)