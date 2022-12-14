################################################################################
##  sCoRRE_functionalDiversityMetrics.R: Calculating functional diversity metrics.
##
##  Authors: Magda Garbowski, Kevin Wilcox, Kimberly Komatsu
##  Date created: April 7, 2021
################################################################################

library(FD)
library(car)
library(tidyverse)

##### set working directory #####
setwd("~/Dropbox/sDiv_sCoRRE_shared/")
setwd("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/") #Padu's wd
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\CoRRE data\\") # Kevin's laptop wd
setwd("C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\") # Kim's laptop/desktop


##### defining functions #####
## Standard Error Function:
se <- function(x, na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}


##### data import and cleaning #####

# continuous trait data
contTraits <- read.csv('CoRRE data\\trait data\\Final TRY Traits\\Imputed Continuous_Traits\\data to play with\\imputed_continuous_20220620.csv') %>%
  select(species_matched, leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass) %>%
  group_by(species_matched) %>%
  summarise_all(funs(mean)) %>%
  ungroup() %>%
  filter(seed_dry_mass<30, plant_height_vegetative<10, rooting_depth<3, SLA<75) 

#testing normality
hist(contTraits$leaf_C.N)
qqPlot(contTraits$leaf_C.N)
shapiro.test(contTraits$leaf_C.N)

hist(contTraits$LDMC)
qqPlot(contTraits$LDMC)
shapiro.test(contTraits$LDMC)

hist(contTraits$SLA)
qqPlot(contTraits$SLA)
shapiro.test(contTraits$SLA)

hist(contTraits$plant_height_vegetative)
qqPlot(contTraits$plant_height_vegetative)
shapiro.test(contTraits$plant_height_vegetative)

hist(contTraits$rooting_depth)
qqPlot(contTraits$rooting_depth)
shapiro.test(contTraits$rooting_depth)

hist(contTraits$seed_dry_mass)
qqPlot(contTraits$seed_dry_mass)
shapiro.test(contTraits$seed_dry_mass)


##### log transform and scale continuous traits #####
contTraitsScaled <- contTraits %>%
  mutate_at(vars(leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass), log) %>% 
  mutate_at(vars(leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass), scale) #scale continuous traits
colnames(contTraits) <- c('species_matched', 'leaf_C.N', 'LDMC', 'SLA', 'plant_height_vegetative', 'rooting_depth', 'seed_dry_mass')

#testing normality
hist(contTraitsScaled$leaf_C.N)
qqPlot(contTraitsScaled$leaf_C.N)
shapiro.test(contTraitsScaled$leaf_C.N)
#sqrt W = 0.66313, p-value < 2.2e-16
#log W = 0.87271, p-value < 2.2e-16


hist(contTraitsScaled$LDMC)
qqPlot(contTraitsScaled$LDMC)
shapiro.test(contTraitsScaled$LDMC)
#sqrt W = 0.99413, p-value = 2.207e-06
#log W = 0.9974, p-value = 0.005625

hist(contTraitsScaled$SLA)
qqPlot(contTraitsScaled$SLA)
shapiro.test(contTraitsScaled$SLA)
#sqrt W = 0.97877, p-value = 2.138e-15
#log W = 0.98822, p-value = 1.033e-10

hist(contTraitsScaled$plant_height_vegetative)
qqPlot(contTraitsScaled$plant_height_vegetative)
shapiro.test(contTraitsScaled$plant_height_vegetative)
#sqrt W = 0.91475, p-value < 2.2e-16
#log W = 0.99106, p-value = 7.733e-09

hist(contTraitsScaled$rooting_depth)
qqPlot(contTraitsScaled$rooting_depth)
shapiro.test(contTraitsScaled$rooting_depth)
#sqrt W = 0.96263, p-value < 2.2e-16
#log W = 0.98769, p-value = 5.035e-11

hist(contTraitsScaled$seed_dry_mass)
qqPlot(contTraitsScaled$seed_dry_mass)
shapiro.test(contTraitsScaled$seed_dry_mass)
#sqrt W = 0.84399, p-value < 2.2e-16
#log W = 0.98691, p-value = 1.769e-11


##### combining continuous and categorical trait datasets #####
traitsAll <- read.csv('CoRRE data\\trait data\\sCoRRE categorical trait data_12142022.csv') %>% #categorical trait data
  select(family, species_matched, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation, leaf_type) %>% 
  filter(growth_form!="moss",species_matched!="", growth_form!="lycophyte") %>% 
  mutate(mycorrhizal=ifelse(mycorrhizal_type=="none", 'no', ifelse(mycorrhizal_type=="uncertain", "unk", "yes"))) %>% 
  select(-mycorrhizal_type) %>% 
  rename(mycorrhizal_type=mycorrhizal) %>% 
  mutate(photo_path=ifelse(photosynthetic_pathway=="possible C4"|photosynthetic_pathway=="possible C4/CAM", "C4", ifelse(photosynthetic_pathway=="possible CAM", "CAM",photosynthetic_pathway))) %>% 
  select(-photosynthetic_pathway) %>%
  rename(photosynthetic_pathway=photo_path) %>%
  full_join(contTraitsScaled) %>% #merge continuous and categorical traits
  drop_na()

# remove non-target species (ferns, lycophytes, mosses) that somehow snuck into the trait data
traitsClean <- traitsAll %>%
  filter(!leaf_type %in% c("microphyll","frond")) %>%
  filter(!species_matched %in% c("Centrolepis aristata", "Centrolepis strigosa", "Acorus calamus"))

# species relative cover data
relcov_full_raw <- read.csv("CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeCover_Dec2021.csv") %>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>%
  dplyr::select(site_code:community_type, site_proj_comm, calendar_year:relcov)

# corre to try species names key
corre_to_try <- read.csv("CoRRE data\\trait data\\corre2trykey_2021.csv") %>%
  dplyr::select(genus_species, species_matched) %>%
  unique(.)

# merge species names and remove all mosses -- moss key to remove mosses from species comp data
moss_sp_vec <- read.csv("CoRRE data\\trait data\\sCoRRE categorical trait data_11302021.csv") %>%
  dplyr::select(species_matched, leaf_type) %>%
  mutate(moss = ifelse(leaf_type=="moss", "moss","non-moss")) %>%
  filter(moss=="moss") %>%
  pull(species_matched)

relcov_full_clean <- relcov_full_raw %>%
  dplyr::left_join(corre_to_try, by="genus_species") %>%
  filter(!species_matched  %in% moss_sp_vec) %>%
  filter(!(site_code == "EGN" & project_name == "Nmow" & community_type == "0" & 
             plot_id == "19" & calendar_year==2015 & treatment=="Control")) ## Stipa species just has an incorrect treatment I think, needs to be fixed in absolute abundance data frame but I'm just removing it for now

rm(moss_key, traits_catag_clean, traits_cont_clean, sp_without_catag_data)


##### calculate functional dispersion - loop through sites #####
FD_df_master <- {}
site_proj_comm_vector <- unique(relcov_full_raw$site_proj_comm)

# for(PROJ in 1:4){
for(PROJ in 1:length(site_proj_comm_vector)){
  relcov_df_temp <-relcov_full_clean %>%
    filter(site_proj_comm==site_proj_comm_vector[PROJ])
  
  #species vector for pulling traits from relative cover
  sp_df_temp <- data.frame(genus_species = unique(relcov_df_temp$genus_species), dummy=1) %>%
    left_join(corre_to_try, by="genus_species") %>%
    unique(.) 
  
  sp_vec_temp <- sp_df_temp %>%
    pull(species_matched)
  
  #subset trait data to just include the species present subset relative cover data
  traits_df_raw_temp <- traitsClean %>%
    filter(species_matched %in% sp_vec_temp)
  
  #dataframe with species in trait database and in relative cover data base
  species_in_trait_data_temp <- data.frame(species_matched = unique(traits_df_raw_temp$species_matched),
                                      dummy_traits=2) %>% #there are fewer species in the unique trait dataset than in the species comp data because there are thing like "unknown forb"
    arrange(species_matched)
  
  #vector of species not in trait database (but in relative abundance data) to remove from species abundance data
  sp_to_remove_temp <- sp_df_temp %>%
    full_join(species_in_trait_data_temp, by="species_matched") %>%
    filter(is.na(dummy_traits)) %>%
    pull(genus_species)
  
  #abundance dataset with species removed that do not have trait information
  relcov_unkn_sp_rm_temp <- relcov_df_temp %>%
    filter(!genus_species %in% sp_to_remove_temp) #removing species without trait information
 
  #abundance data into wide format
  relcov_wide_temp <- relcov_unkn_sp_rm_temp %>%
    dplyr::select(-genus_species) %>%
    group_by(site_code, project_name, community_type, site_proj_comm, calendar_year, treatment_year, treatment,
             block, plot_id, data_type, version, species_matched) %>%
    summarize(relcov=sum(relcov, na.rm=T)) %>%
    ungroup() %>%
    spread(key=species_matched, value=relcov) %>%
    replace(is.na(.), 0)
  
  plot_info_temp <- relcov_wide_temp %>%
    dplyr::select(site_code:version)
  
  relcov_only_temp <- relcov_wide_temp %>%
    dplyr::select(-site_code:-version) 
  
  # change to dataframe from tibble 
  relcov_only_temp <- as.data.frame(relcov_only_temp)
  
  # add rownames 
  row.names(relcov_only_temp) <- paste(plot_info_temp$calendar_year, plot_info_temp$plot_id, sep="::")

  #dbFD function requires species names in trait data frame be arranged A-Z and identical order to the abundance data 
  traits_df_temp <- traits_df_raw_temp %>%
    arrange(species_matched) %>%
    column_to_rownames("species_matched") %>%
    dplyr::select(-family) %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
    select(growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation, leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass)

  # change to dataframe from tibble 
  traits_df_temp <- as.data.frame(traits_df_temp)
  
  #changing all categorical traits to factors
  traits_df_temp[,c(1:6)] <- lapply(traits_df_temp[,c(1:6)], as.factor)
  
  #changing all continuous to numerical
  traits_df_temp[,c(7:12)] <- lapply(traits_df_temp[,c(7:12)], as.numeric)
  
  ### Calculate functional diversity metrics -- had to use Cailliez correlations becuase Euclidean distances could be calculated
  FD_temp <- dbFD(x=traits_df_temp, a=relcov_only_temp, cor="cailliez", calc.FRic=F) # FRich is causing problems with most datasets (I think because of missing data?) so I'm removing it for now
  
  FD_df_temp <- do.call(cbind.data.frame, FD_temp) %>%
    mutate(year_plotid = row.names(.)) %>%
    separate(year_plotid, into=c("calendar_year","plot_id"), sep="::") %>%
    mutate(calendar_year = as.numeric(calendar_year)) %>%
    full_join(plot_info_temp, by=c("calendar_year","plot_id"))
  
  FD_df_master <- rbind(FD_df_master, FD_df_temp)
  
  rm(list=ls()[grep("temp", ls())])
}


# write.csv(FD_df_master, 'C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\paper 2_PD and FD responses\\CoRRE_functionalDiversity_2022-12-13.csv',row.names=F)


