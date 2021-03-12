### 
### Calculate Functional dispersion for CoRRE data
### 

### Set up workspace
rm(list=ls())

setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\")
library(FD)
library(tidyverse)
memory.limit(size=50000)

### Trait data
#read in data
contTraits <- read.csv('Trait Data\\TRY Data\\Gap_Filled\\TRY_new.csv')%>%
  rename(species_matched=Species)%>%
  select(-X.1, -X, -Family, -Genus, -ObservationID)%>%
  group_by(species_matched)%>%
  summarise_all(funs(mean))%>%
  ungroup()

contTraitsSubset <- contTraits%>%
  rename(ssd=X4, rooting_depth=X6, SLA=X11, leaf_C_mass=X13, leaf_N_mass=X14, leaf_P_mass=X15, stem_diameter=X21, seed_mass=X26, seed_length=X27, leaf_thickness=X46, LDMC=X47, leaf_dry_mass=X55, germination_rate=X95, leaf_length=X144, leaf_width=X145, leaf_CN=X146, stem_conduit_density=X169, stem_conduit_diameter=X281, seed_number=X138, SRL=X1080)%>%
  select(-X18, -X50, -X78, -X163, -X223, -X224, -X237, -X282, -X289, -X3112, -X3113, -X3114, -X3120)

traits <- read.csv('CoRRE data\\CoRRE data\\trait data\\sCoRRE categorical trait data - traits_complete_pre spot check_03102021.csv')%>%
  full_join(contTraitsSubset) %>%
  drop_na()

traitsOutliersRemoved <- traits %>%
  filter(!leaf_type %in% c("microphyll","frond")) %>%
  filter(!species_matched %in% c("Centrolepis aristata", "Centrolepis strigosa", "Acorus calamus"))

traitsScaled <- traitsOutliersRemoved %>% ## only scales continuous traits
  mutate_at(vars(ssd:SRL), scale)

# Read in relative abundance data
sp_name_key <- read.csv("CoRRE data\\CoRRE data\\trait data\\corre2trykey.csv") %>%
  dplyr::select(genus_species, species_matched) %>%
  unique(.)

rel_abun_df_raw <- read.csv("CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeAbundanceMar2021.csv") %>%
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

pplot_sp_wide <- pplot_abun_df %>%
  dplyr::select(-genus_species) %>%
  spread(key=species_matched, value=relcov) %>%
  dplyr::select(-site_code:-version)

pplot_env_df <- pplot_abun_df %>%
  dplyr::select(-genus_species) %>%
  spread(key=species_matched, value=relcov) %>%
  dplyr::select(site_code:version)

write.csv(pplot_sp_wide, file="paper 2_PD and FD responses\\data\\pplots species abundance wide.csv")
write.csv(pplot_env_df, file="paper 2_PD and FD responses\\data\\pplots environmental data.csv")

dups <- rel_abun_df %>% 
  group_by(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block,
                                 plot_id, data_type, version, genus_species) %>%
  summarize(length=length(relcov)) %>%
  filter(length==2)
  
write.csv(dups, file="paper 2_PD and FD responses\\data\\duplicates.csv")
  
  
  #Experimental information
# C:\Users\wilco\Dropbox\shared working groups\sDiv_sCoRRE_shared\CoRRE data\CoRRE data\community composition\CoRRE_ExperimentInfoMar2021.csv

# for matched_names (key for matching CORRE and TRY species names)
#C:\Users\wilco\Dropbox\shared working groups\sDiv_sCoRRE_shared\CoRRE data\CoRRE data\trait data\corre2trykey.csv

### 

pplot_abun_df <- read.csv("CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeAbundanceMar2021.csv") %>%
  filter(project_name=="pplots") %>%
  left_join(sp_name_key, by="genus_species") %>%
  drop_na(species_matched)



