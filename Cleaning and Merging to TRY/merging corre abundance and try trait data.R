###
### Combine CoRRE abundance with traits
### Author: Wilcox (kevin.wilcox@uwyo.edu)
### Created Dec 2019

### Set up workspace
# Kevin's laptop
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\CoRRE data\\")

library(tidyverse)
library(ggthemes)

### Read in data
raw_abundance_df <- read.csv("CoRRE_raw_abundance_Nov2019.csv") %>%
  dplyr::select(-X)
trait_df <- read.csv("TRY_trait_data_continuous.csv")
species_name_key <- read.csv("CoRRE_TRY_species_list.csv") %>%
  dplyr::select(-X)

### Combine trait data with abundance data
abundance_traits_df <- raw_abundance_df %>%
  full_join(species_name_key, by="genus_species") %>%
  full_join(trait_df, by="species_matched")





