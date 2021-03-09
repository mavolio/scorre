library(ggtree)
library(ape)
library(googledrive)
library(readxl)
library(tidyverse)
library(phytools)
library(pez)
vcapply <- function(X, FUN, ...) {
  vapply(X, FUN, character(1), ...)
}

split_genus <- function(str) {
  str_split <- strsplit(str, "[_ ]+")
  vcapply(str_split, "[[", 1L)
}

#read a tree
phylo <- read.tree("data/zanne1.1.tre")

#read in categoriacl data
if(!file.exists("data/sCoRRE_categorical_trait_data.xlsx")){
  drive_download("sCoRRE categorical trait data",path="data/sCoRRE_categorical_trait_data.xlsx")
}



read_xlsx("data/sCoRRE_categorical_trait_data.xlsx", sheet = 1) %>%
  filter(!is.na(species_matched)) %>%
  mutate(species_matched = gsub("\\s", replacement = "_",x = species_matched),
         genus = split_genus(species_matched)) -> categorical_data

categorical_data %>%
  select(genus, photosynthetic_pathway) %>%
  filter(!is.na(genus)) %>%
  count(genus) %>%
  mutate(phot_path = ifelse(n>0, 1,0)) %>%
  select(-n)->trait_data

data_spp <- data.frame(species=phylo$tip.label)
data_spp$genus <- split_genus(phylo$tip.label)

spp_selected <- data_spp[!duplicated(data_spp$genus),]
spp_delete <- data_spp[!(data_spp$species %in% spp_selected$species),]
spp_delete <- as.character(spp_delete$species)
phyl_collapse <- pez::drop_tip(phylo,spp = spp_delete)
phyl_collapse$tip.label <- spp_selected$genus

ancs <-mrca(phyl_collapse, full = FALSE)
#mrca(phyl_collapse)
#MRCA(phyl_collapse,c("Alismataceae","Zingiberales"))
#MRCA(phyl_collapse,c("Ceratophyllum","Adoxa"))

trait_data <- as.data.frame(trait_data)
rownames(trait_data) <- trait_data$genus
trait_data <- trait_data[-1,]


pdf("plant_phylo.pdf",width=11,height=8.5)
ggg <- ggtree(phyl_collapse,layout="radial")
g <-gheatmap(ggg, trait_data, offset = 0.1, width=0.3,low = "blue",high = "red", colnames = FALSE)+
  scale_fill_viridis_d()+
  geom_cladelabel(node=10061,label="Magno",color = "purple",fontsize = 10)+
  geom_cladelabel(node=8339,label="Mono",color = "red",fontsize = 10)+
  geom_cladelabel(node=10259,label="Eudicots",color = "blue",fontsize = 10)+
  theme(legend.position = "none")

print(g)
dev.off()
