### Cleaning continuous trait data
### 
### Authors: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Created: December 14, 2021; last updated: December 14, 2021
### Created with R version 4.1.1
###

#### NOTES
###
### Continuous traits to include: LDMC, SLA, Vegetative_height, seed dry mass, seed number, rooting density, rooting depth 

### Set up working space
# rm(list=ls()) clean up workspace
#library(FD)
library(tidyverse)

# setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\CoRRE data\\") # Kevin's laptop wd

### Read in trait data -- for now I'm just cleaning the traits identified above
traits_cont_raw <- read.csv("CoRRE data\\trait data\\Final Cleaned Traits\\Continuous_Traits\\Backtrans_GapFilled_sCorre.csv") %>%
  dplyr::select(X:family, LDMC, SLA, plant_height_vegetative, seed_dry_mass, seed_number, rooting_depth, SRL) %>%
  mutate(across(everything(), ~replace(., .<0, NA)))


# look at historgrams (hist code is working but I want to source this code from other scripts and don't want plotting code to run)
# par(mfrow=c(3,3))
# for(COL in 6:ncol(traits_cont_raw)){
#   hist(traits_cont_raw[,COL], xlab=colnames(traits_cont_raw[COL]))
# }

### Read in categorical data to get moss column
moss_key <- read.csv("CoRRE data\\trait data\\Final Cleaned Traits\\sCoRRE categorical trait data_final_20211209.csv") %>%
  dplyr::select(species_matched, leaf_type) %>%
  mutate(moss = ifelse(leaf_type=="moss", "moss","non-moss")) %>%
  dplyr::select(-leaf_type)

### Merge moss column into continuous trait data

### Removing unrealistically high values for all traits -- based on trait histograms
traits_cont_clean <- traits_cont_raw %>%
  mutate(LDMC = replace(LDMC, LDMC > 0.8, NA), 
         SLA = replace(SLA, SLA > 100, NA),
         plant_height_vegetative = replace(plant_height_vegetative, plant_height_vegetative > 5, NA),
         seed_dry_mass = replace(seed_dry_mass, seed_dry_mass > 100, NA),
         seed_number = replace(seed_number, seed_number > 1e5, NA),
         rooting_depth = replace(rooting_depth, rooting_depth > 4.5, NA),
         SRL = replace(SRL, SRL > 60000, NA)) %>%
  group_by(genus, family, species_matched) %>%
 # summarize(LDMC=mean(LDMC,na.rm=T)) %>%
  summarize_at(vars(LDMC:SRL), list(mean=mean, sd=sd), na.rm=T) %>%
  ungroup() %>%
  left_join(moss_key, by="species_matched") %>%
  filter(moss!="moss") %>%
  dplyr::select(-moss)
  

# hist(traits_cont_clean$LDMC_mean)
# hist(traits_cont_clean$SLA_mean)
# hist(traits_cont_clean$plant_height_vegetative_mean)
# hist(traits_cont_clean$seed_dry_mass_mean)
# hist(traits_cont_clean$seed_number_mean)
# hist(traits_cont_clean$rooting_depth_mean)

rm(traits_cont_raw)
