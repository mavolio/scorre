### 
### Calculate Functional dispersion for CoRRE data
### 
### Last updated Dec 13, 2021

### Set up workspace
rm(list=ls())

setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\") #kevin's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\') #kim's laptop

library(FD)
library(tidyverse)
memory.limit(size=50000)

### Trait data
# Read in data
contTraits <- read.csv('CoRRE data\\trait data\\Final TRY Traits\\Imputed Continuous_Traits\\data to play with\\imputed_continuous_20220620.csv')%>%
  select(-X.1, -X, -family, -genus, -observation)%>%
  group_by(species_matched)%>%
  summarise_all(funs(mean))%>%
  ungroup()

traits <- read.csv('CoRRE data\\trait data\\sCoRRE categorical trait data - traits_complete_pre spot check_03102021.csv')%>%
  full_join(contTraits) %>%
  drop_na()

traitsOutliersRemoved <- traits %>%
  filter(!leaf_type %in% c("microphyll","frond")) %>%
  filter(!species_matched %in% c("Centrolepis aristata", "Centrolepis strigosa", "Acorus calamus"))

traitsScaled <- traitsOutliersRemoved %>% ## only scales continuous traits
  mutate_at(vars(seed_dry_mass:seed_number), scale)

# Read in relative abundance data
sp_name_key <- read.csv("CoRRE data\\trait data\\corre2trykey.csv") %>%
  dplyr::select(genus_species, species_matched) %>%
  unique(.)

rel_abun_df_raw <- read.csv("CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeCover_Dec2021.csv") %>%
  group_by(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block,
           plot_id, data_type, version, genus_species) %>%
  summarize(relcov=max(relcov))

rel_abun_df <- rel_abun_df_raw %>%
  left_join(sp_name_key, by="genus_species") %>%
  filter(!is.na(species_matched)) %>%
  dplyr::select(-genus_species) %>%
  group_by(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block,
           plot_id, data_type, version, species_matched) %>%
  summarize(relcov=sum(relcov))
  


# filter(rel_abun_df, site_code=="SIU" & plot_id==48 & treatment_year==5 & genus_species=="acalypha virginica")
# abund_species_vector <- unique(rel_abun_df$species_matched)
rel_abun_df[c(13080, 13081),]
rel_abun_env_df <- rel_abun_df %>%
  spread(key=species_matched, value=relcov) %>%
  dplyr::select(site_code:version)

rel_abun_sp_wide <- rel_abun_df %>%
  spread(key=species_matched, value=relcov) %>%
  dplyr::select(-site_code:-version)

# pplot_sp_wide <- pplot_abun_df %>%
#   dplyr::select(-genus_species) %>%
#   spread(key=species_matched, value=relcov) %>%
#   dplyr::select(-site_code:-version)
# 
# pplot_env_df <- pplot_abun_df %>%
#   dplyr::select(-genus_species) %>%
#   spread(key=species_matched, value=relcov) %>%
#   dplyr::select(site_code:version)
# 
# write.csv(pplot_sp_wide, file="paper 2_PD and FD responses\\data\\pplots species abundance wide.csv")
# write.csv(pplot_env_df, file="paper 2_PD and FD responses\\data\\pplots environmental data.csv")

dups <- rel_abun_df %>% 
  group_by(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block,
                                 plot_id, data_type, version, species_matched) %>%
  summarize(length=length(relcov)) %>%
  filter(length==2)
  
# write.csv(dups, file="paper 2_PD and FD responses\\data\\duplicates.csv") #no duplicates anymore
  
  
  #Experimental information
# C:\Users\wilco\Dropbox\shared working groups\sDiv_sCoRRE_shared\CoRRE data\CoRRE data\community composition\CoRRE_ExperimentInfoMar2021.csv

# for matched_names (key for matching CORRE and TRY species names)
#C:\Users\wilco\Dropbox\shared working groups\sDiv_sCoRRE_shared\CoRRE data\CoRRE data\trait data\corre2trykey.csv

### 

pplot_abun_df <- read.csv("CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeAbundanceMar2021.csv") %>%
  filter(project_name=="pplots") %>%
  left_join(sp_name_key, by="genus_species") %>%
  drop_na(species_matched)



