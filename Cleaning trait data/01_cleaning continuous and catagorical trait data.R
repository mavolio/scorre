### Cleaning continuous trait data
### 
### Authors: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Created: December 14, 2021; last updated: December 14, 2021
### Created with R version 4.1.1
###

###################################
### Set up working space
###################################

# rm(list=ls()) clean up workspace
#library(FD)
library(tidyverse)

# setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\CoRRE data\\") # Kevin's laptop wd


##################################
### Continous traits
##################################

#### NOTES
###
### Continuous traits to include: LDMC, SLA, Vegetative_height, seed dry mass, seed number, rooting density, rooting depth 


### Read in trait data -- for now I'm just cleaning the traits identified above
traits_cont_raw <- read.csv("CoRRE data\\trait data\\Final Cleaned Traits\\Continuous_Traits\\Backtrans_GapFilled_sCorre.csv") %>%
  dplyr::select(X:family, LDMC, SLA, plant_height_vegetative, seed_dry_mass, seed_number, rooting_depth) %>% # NEED TO REMOVE SRL
  mutate(across(everything(), ~replace(., .<0, NA)))


### Look at histograms to determine cutoffs for outlier trait removal

# hist(traits_cont_raw$LDMC)
# hist(traits_cont_raw$SLA)
# hist(traits_cont_raw$plant_height_vegetative)
# hist(traits_cont_raw$seed_dry_mass)
# hist(traits_cont_raw$seed_number)
# hist(traits_cont_raw$rooting_depth)
# filter(traits_cont_raw, SLA == max(SLA))
# filter(traits_cont_raw, SLA > 200)

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
         rooting_depth = replace(rooting_depth, rooting_depth > 4.5, NA)) %>%
  group_by(genus, family, species_matched) %>%
 # summarize(LDMC=mean(LDMC,na.rm=T)) %>%
  summarize_at(vars(LDMC:rooting_depth), list(mean=mean, sd=sd), na.rm=T) %>%
  ungroup() %>%
  left_join(moss_key, by="species_matched") %>%
  mutate(moss=ifelse(moss=="moss","moss","non-moss")) %>%
  filter(moss!="moss") %>%
  dplyr::select(-moss)


# hist(traits_cont_clean$LDMC_mean)
# hist(traits_cont_clean$SLA_mean)
# hist(traits_cont_clean$plant_height_vegetative_mean)
# hist(traits_cont_clean$seed_dry_mass_mean)
# hist(traits_cont_clean$seed_number_mean)
# hist(traits_cont_clean$rooting_depth_mean)

rm(traits_cont_raw)

#########################################################################
### Categorical trait data 
#########################################################################

#### NOTES:
###
### Categorical traits to include: growth_form, life_span, mycorrhizal_type, n_fixation, clonal, photosynthetic_pathway 
###
#########


#all_traits_raw <- openxlsx::read.xlsx("sCoRRE categorical trait data_final_20211209.xlsx", sheet=1) %>%
traits_catag_clean <- read.csv("CoRRE data\\trait data\\Final Cleaned Traits\\sCoRRE categorical trait data_final_20211209.csv") %>%
  dplyr::select(species_matched, growth_form, photosynthetic_pathway, lifespan,  clonal, mycorrhizal_type, n_fixation) %>%
  mutate(photosynthetic_pathway = replace(photosynthetic_pathway, grep("possible", photosynthetic_pathway), NA)) %>%
  mutate(clonal = replace(clonal, clonal=="uncertain", NA)) %>%
  mutate(mycorrhizal_type = replace(mycorrhizal_type, mycorrhizal_type=="uncertain", NA)) %>%
  mutate(lifespan = replace(lifespan, lifespan=="uncertain", NA)) %>%
  filter(lifespan != "moss")


###########################################################
### Combine continuous and catagorical traits
###########################################################

traits_all <- dplyr::select(traits_cont_clean, -LDMC_sd:-rooting_depth_sd) %>%
  full_join(traits_catag_clean, by="species_matched")

# Find species that have continuous data but no categorical data -- 3 of these are mosses, so I will remove them from continuous data in that cleaning script
# TO DO: We will want to populate categorical data for three species: "Galium mollugo" "Heracleum sphondylium" "Trachypogon spicatus"


