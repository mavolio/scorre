###
### Check shared species between CoRRE and GEx
### Author: Wilcox (kevin.wilcox@uwyo.edu); created Dec 2019
###

library(tidyverse)
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\CoRRE data\\") # Kevin's laptop

### Read in species lists
corre_species_list_full <- read.csv("CoRRE_TRY_species_list.csv") 
gex_species_vector <- read.csv("GEx_species names from TRY-BIEN merged_Aug 13 2019.csv") %>%
  rename(gex_present=x, spname=sp_names)
gex_species_full <- read.csv("GEx_All_Species_Trait.csv")

gex_species_vector <- data.frame(spname = unique(gex_species_full$clean_ejf),
                                 gex_present = 1)
corre_species_vector <- data.frame(spname = unique(corre_species_list_full$species_matched), 
                                   corre_present = 1)

### Match GEx and TRY species lists
corre_gex_matched_full <- corre_species_vector %>%
  full_join(gex_species_vector, by="spname") %>%
  mutate(gex_present = replace_na(gex_present, 0)) %>%
  mutate(corre_present = replace_na(corre_present, 0)) %>%
  mutate(both_present = ifelse(gex_present==1&corre_present==1, 1, 0)) %>%
  

colSums(corre_gex_matched_full[2:4])


corre_gex_matched_left <- corre_species_vector %>%
  left_join(gex_species_vector, by="spname") %>%
  mutate(gex_present = replace_na(gex_present, 0))


colSums(corre_gex_matched_left[2:4])



