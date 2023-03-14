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
setwd("C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\") # Kim's laptop/desktop
setwd("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/") #Padu's wd
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Working groups\\sDiv\\") # Kevin's laptop wd


##### defining functions #####
## Standard Error Function:
se <- function(x, na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}


##### data import and cleaning #####

# trait data
traits <- read.csv('CoRRE data\\trait data\\AllTraits\\CoRRE_allTraitData_March2023.csv') %>% 
  select(family, species_matched, leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation) %>% 
  filter(growth_form!="moss", species_matched!="") %>% #keeping lycophytes
  mutate(mycorrhizal=ifelse(mycorrhizal_type=="none", 'no', ifelse(mycorrhizal_type=="uncertain", "unk", "yes"))) %>% 
  select(-mycorrhizal_type) %>% 
  rename(mycorrhizal_type=mycorrhizal) %>% 
  mutate(photo_path=ifelse(photosynthetic_pathway=="possible C4"|photosynthetic_pathway=="possible C4/CAM", "C4", ifelse(photosynthetic_pathway=="possible CAM", "CAM",photosynthetic_pathway))) %>% 
  select(-photosynthetic_pathway) %>%
  rename(photosynthetic_pathway=photo_path) %>%
  drop_na() #only keep trait data that is complete for all traits (drops 600 species with only categorical trait data)


##### testing normality #####
hist(traits$leaf_C.N)
qqPlot(traits$leaf_C.N)
shapiro.test(traits$leaf_C.N)

hist(traits$LDMC)
qqPlot(traits$LDMC)
shapiro.test(traits$LDMC)

hist(traits$SLA)
qqPlot(traits$SLA)
shapiro.test(traits$SLA)

hist(traits$plant_height_vegetative)
qqPlot(traits$plant_height_vegetative)
shapiro.test(traits$plant_height_vegetative)

hist(traits$rooting_depth)
qqPlot(traits$rooting_depth)
shapiro.test(traits$rooting_depth)

hist(traits$seed_dry_mass)
qqPlot(traits$seed_dry_mass)
shapiro.test(traits$seed_dry_mass)


##### log transform and scale continuous traits #####
traitsScaled <- traits %>%
  mutate_at(vars(leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass), log) %>% 
  mutate_at(vars(leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass), scale) #scale continuous traits
colnames(contScaled) <- c('family', 'species_matched', 'leaf_C.N', 'LDMC', 'SLA', 'plant_height_vegetative', 'rooting_depth', 'seed_dry_mass')

#testing normality
hist(traitsScaled$leaf_C.N)
qqPlot(traitsScaled$leaf_C.N)
shapiro.test(traitsScaled$leaf_C.N)
#log W = 0.97422, p-value < 2.2e-16
#sqrt W = 0.89825, p-value < 2.2e-16

hist(traitsScaled$LDMC)
qqPlot(traitsScaled$LDMC)
shapiro.test(traitsScaled$LDMC)
#log W = 0.94689, p-value < 2.2e-16
#sqrt W = 0.97428, p-value < 2.2e-16

hist(traitsScaled$SLA)
qqPlot(traitsScaled$SLA)
shapiro.test(traitsScaled$SLA)
#log W = 0.96329, p-value < 2.2e-16
#sqrt W = 0.92173, p-value < 2.2e-16

hist(traitsScaled$plant_height_vegetative)
qqPlot(traitsScaled$plant_height_vegetative)
shapiro.test(traitsScaled$plant_height_vegetative)
#log W = 0.99396, p-value = 1.029e-06
#sqrt W = 0.82253, p-value < 2.2e-16

hist(traitsScaled$rooting_depth)
qqPlot(traitsScaled$rooting_depth)
shapiro.test(traitsScaled$rooting_depth)
#log W = 0.98391, p-value = 2.411e-13
#sqrt W = 0.8344, p-value < 2.2e-16

hist(traitsScaled$seed_dry_mass)
qqPlot(traitsScaled$seed_dry_mass)
shapiro.test(traitsScaled$seed_dry_mass)
#log W = 0.9931, p-value = 1.846e-07
#sqrt W = 0.7173, p-value < 2.2e-16


##### relative cover datasets #####

# species relative cover data
relcov_full_raw <- read.csv("CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeCover_Jan2023.csv") %>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>%
  dplyr::select(site_code:community_type, site_proj_comm, calendar_year:relcov)

# corre to try species names key
corre_to_try <- read.csv("CoRRE data\\trait data\\corre2trykey_2021.csv") %>%
  dplyr::select(genus_species, species_matched) %>%
  unique(.)

# merge species names and remove all mosses -- moss key to remove mosses from species comp data
moss_sp_vec <- read.csv("CoRRE data\\trait data\\sCoRRE categorical trait data_12142022.csv") %>%
  dplyr::select(species_matched, leaf_type) %>%
  mutate(moss = ifelse(leaf_type=="moss", "moss","non-moss")) %>%
  filter(moss=="moss") %>%
  pull(species_matched)

relcov_full_clean <- relcov_full_raw %>%
  dplyr::left_join(corre_to_try, by="genus_species") %>%
  filter(!species_matched  %in% moss_sp_vec) %>%
  mutate(plot_id=ifelse(site_proj_comm=='DL_NSFC_0', paste(plot_id, treatment, sep='__'), plot_id))

rm(moss_sp_vec, traits_catag_clean, traits_cont_clean, sp_without_catag_data)


##### calculate functional dispersion - loop through sites #####
distance_df_master <- {}
site_proj_comm_vector <- unique(relcov_full_raw$site_proj_comm)
site_proj_comm_vector[-103] #PIE TIDE is causing errors with FRich for some reason

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
  traits_df_raw_temp <- traitsScaled %>%
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
  
  ### Calculate MNTD and functional diversity metrics
  FD_temp <- dbFD(x=traits_df_temp, # matrix of traits
                  a=relcov_only_temp, # matrix of species abundances
                  w.abun=F, # don't weight by abundance
                  cor="cailliez", # use Cailliez correlations because Euclidean distances could be calculated
                  # calc.CWM=T, CWM.type='all', # calculate CWM across all spp
                  calc.FRic=T) # FRich is causing problems with most datasets (I think because of missing data?) so I'm removing it for now
  
  FD_df_temp <- do.call(cbind.data.frame, FD_temp) %>%
    mutate(year_plotid = row.names(.)) %>%
    separate(year_plotid, into=c("calendar_year","plot_id"), sep="::") %>%
    mutate(calendar_year = as.numeric(calendar_year)) %>%
    full_join(plot_info_temp, by=c("calendar_year","plot_id"))
  
  comp_matrix_temp <- as.matrix(relcov_only_temp)
  trait_dist_temp <- as.matrix(gowdis(traits_df_temp))

  MNTD_df_temp <- data.frame(
    plot_info_temp[,c("calendar_year", "plot_id")],
    MNTD_traits = picante::mntd(comp_matrix_temp, trait_dist_temp)
    )
  
 distance_df_temp <- FD_df_temp %>%
   full_join(MNTD_df_temp, by=c("calendar_year","plot_id"))
 
 distance_df_master <- rbind(distance_df_master, distance_df_temp)
 
  rm(list=ls()[grep("temp", ls())])
}


# write.csv(distance_df_master, 'C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\paper 2_PD and FD responses\\data\\CoRRE_functionalDiversity_2023-03-10.csv',row.names=F)


