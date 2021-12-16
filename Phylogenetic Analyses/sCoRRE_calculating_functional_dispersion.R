### Calculating functional diversity and dispersion using categorical and continuous traits
### Last updated: Dec 13, 2021
### R version: 4.1.1

library(FD)
library(tidyverse)
library(beepr)
# install.packages("mFD")
# library("mFD")
# rm(list=ls())

## testing github
setwd("~/Dropbox/sDiv_sCoRRE_shared/")
setwd("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/") #Padu's wd
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\CoRRE data\\") # Kevin's laptop wd

# Standard Error Function:
se <- function(x, na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}

###
### Read in and prep all data
###

### Read in a clean continuous trait data
source("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Git projects\\scorre\\Cleaning trait data\\01_cleaning continuous trait data.R")

### Categorical trait data 
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

traits_all <- dplyr::select(traits_cont_clean, -LDMC_sd:-rooting_depth_sd) %>%
  full_join(traits_catag_clean, by="species_matched")

# Find species that have continuous data but no categorical data -- 3 of these are mosses, so I will remove them from continuous data in that cleaning script
# TO DO: We will want to populate categorical data for three species: "Galium mollugo" "Heracleum sphondylium" "Trachypogon spicatus"

sp_without_catag_data <- filter(traits_all, is.na(n_fixation))

### full species relative cover data
relcov_full_raw <- read.csv("CoRRE data\\community composition\\CoRRE_RelativeCover_Dec2021.csv") %>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>%
  dplyr::select(site_code:community_type, site_proj_comm, calendar_year:relcov)

### corre to try species names key
corre_to_try <- read.csv("CoRRE data\\trait data\\corre2trykey_2021.csv") %>%
  dplyr::select(genus_species, species_matched) %>%
  unique(.)

### merge species names and remove all mosses

# moss key to remove mosses from species comp data
moss_sp_vec <- read.csv("CoRRE data\\trait data\\Final Cleaned Traits\\sCoRRE categorical trait data_final_20211209.csv") %>%
  dplyr::select(species_matched, leaf_type) %>%
  mutate(moss = ifelse(leaf_type=="moss", "moss","non-moss")) %>%
  filter(moss=="moss") %>%
  pull(species_matched)

#Heracleum sphondylium

relcov_full_clean <- relcov_full_raw %>%
  dplyr::left_join(corre_to_try, by="genus_species") %>%
  filter(!species_matched  %in% moss_sp_vec) %>%
  filter(!(site_code == "EGN" & project_name == "Nmow" & community_type == "0" & 
             plot_id == "19" & calendar_year==2015 & treatment=="Control")) ## Stipa species just has an incorrect treatment I think, needs to be fixed in absolute abundance data frame but I'm just removing it for now
#  filter(relcov_full_clean, species_matched=="Homalothecium pinnatifidum")

# 
# test5 <- relcov_full_raw %>%
#   dplyr::left_join(corre_to_try, by="genus_species") %>%
#   filter(!species_matched  %in% moss_sp_vec) %>%
#   filter(site_code == "EGN" & project_name == "Nmow" & community_type == "0" & 
#              plot_id == "19" & calendar_year==2015 )

rm(moss_key, traits_catag_clean, traits_cont_clean, sp_without_catag_data)

###
### looping through site_proj_comm and calculating FD for all
###

FD_df_master <- {}
site_proj_comm_vector <- unique(relcov_full_raw$site_proj_comm)

for(PROJ in 1:length(site_proj_comm_vector)){
  relcov_df_temp <-relcov_full_clean %>%
    filter(site_proj_comm==site_proj_comm_vector[PROJ])
  
  ### Get species vector for pulling traits from pplots relative cover
  sp_df_temp <- data.frame(genus_species = unique(relcov_df_temp$genus_species), dummy=1) %>%
    left_join(corre_to_try, by="genus_species") %>%
    unique(.) 
  
  sp_vec_temp <- sp_df_temp %>%
    pull(species_matched)
  
  ### Subset trait data to just include the species present in the pplots relative cover data
  traits_df_raw_temp <- traits_all %>%
    filter(species_matched %in% sp_vec_temp)
  
  # test1 <- sp_df_temp %>%
  #   full_join(traits_df_raw_temp, by="species_matched")
  
  ### Get data frame with species in trait data base and in pplots abundance data base
  species_in_trait_data_temp <- data.frame(species_matched = unique(traits_df_raw_temp$species_matched),
                                      dummy_traits=2) %>% ## there are less species in the unique trait dataset than in the species comp data because they're things like "unknown forb"
    arrange(species_matched)
  ### Get vector of species not in trait database (but in relative abundance data) to remove from species abundance data
  sp_to_remove_temp <- sp_df_temp %>%
    full_join(species_in_trait_data_temp, by="species_matched") %>%
    filter(is.na(dummy_traits)) %>%
    pull(genus_species)
  
  ### Abundance data set with species removed that do not have trait information
  relcov_unkn_sp_rm_temp <- relcov_df_temp %>%
    filter(!genus_species %in% sp_to_remove_temp) # removing species without trait information
 
  # relcov_unkn_sp_rm_temp <- relcov_df_temp %>% ### this one is for testing Alberta data set.. checking whether having NAs for continuous data is breaking things
  #   filter(!genus_species %in% sp_to_remove_temp) %>% # removing species without trait information
  #   filter(!species_matched %in% c("Potentilla bipinnatifida", "Festuca hallii"))

  ### get abundance data wide format
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
  
  row.names(relcov_only_temp) <- paste(plot_info_temp$calendar_year, plot_info_temp$plot_id, sep="::")
  
  ### dbFD function requires species names in trait data frame be arranged A-Z and identical order to the abundance data 
  traits_df_temp <- traits_df_raw_temp %>%
    arrange(species_matched) %>%
    column_to_rownames("species_matched") %>%
    dplyr::select(-genus, -family) %>%
    mutate_all(~ifelse(is.nan(.), NA, .))

   ### Changing all categorical traits to factors
  traits_df_temp$growth_form <- as.factor(traits_df_temp$growth_form)
  traits_df_temp$photosynthetic_pathway <- as.factor(traits_df_temp$photosynthetic_pathway)
  traits_df_temp$lifespan <- as.factor(traits_df_temp$lifespan)
  traits_df_temp$clonal <- as.factor(traits_df_temp$clonal)
  traits_df_temp$mycorrhizal_type <- as.factor(traits_df_temp$mycorrhizal_type)
  traits_df_temp$n_fixation <- as.factor(traits_df_temp$n_fixation)
  
  ### create distance matrix for incorporation into dbFD function
  gowdis_temp <- gowdis(traits_df_temp)
  
  ### Calculate functional diversity metrics -- had to use Cailliez correlations becuase Euclidean distances could be calculated
  FD_temp <- dbFD(x=gowdis_temp, a=relcov_only_temp, cor="cailliez", calc.FRic=F) # FRich is causing problems with most datasets (I think because of missing data?) so I'm removing it for now
  
  FD_df_temp <- do.call(cbind.data.frame, FD_temp) %>%
    mutate(year_plotid = row.names(.)) %>%
    separate(year_plotid, into=c("calendar_year","plot_id"), sep="::") %>%
    mutate(calendar_year = as.numeric(calendar_year)) %>%
    full_join(plot_info_temp, by=c("calendar_year","plot_id"))
  
  FD_df_master <- rbind(FD_df_master, FD_df_temp)
  
  rm(list=ls()[grep("temp", ls())])
  #beep(sound = 2, expr = NULL)
}

beep(sound=4)
write.csv(FD_df_master, file=paste0("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Working groups\\sDiv\\Dec2021\\Functional diversity metrics_", Sys.Date(),".csv"),
                                    row.names=F)
  xz
ggplot(filter(FD_df_master, site_proj_comm == site_proj_comm_vector[3]), aes(x=treatment_year, y=FDis, col=treatment)) +
  geom_point()
  
FD_df_means <- FD_df_master %>%
  group_by(site_proj_comm, site_code, project_name, community_type, treatment_year, calendar_year, treatment) %>%
  summarize_at(vars(FEve, FDis, RaoQ), list(mean=mean, se=se), na.rm=T)

### fun plotting
# F div
for(PROJ in 1:length(site_proj_comm_vector)){
  ggplot(filter(FD_df_means, site_proj_comm == site_proj_comm_vector[PROJ]), aes(x=treatment_year, y=FDis_mean, col=treatment,
                                                                              ymin=FDis_mean-FDis_se, ymax=FDis_mean+FDis_se)) +
    geom_errorbar(width=0.1) +
    geom_point() +
    geom_path() +
    ggtitle(site_proj_comm_vector[PROJ]) +
    theme_bw()
  ggsave(filename=paste0("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Working groups\\sDiv\\Dec2021\\figures\\FDis\\",
                         site_proj_comm_vector[PROJ], "FDis.jpg"))
}

# F evenness
for(PROJ in 1:length(site_proj_comm_vector)){
  ggplot(filter(FD_df_means, site_proj_comm == site_proj_comm_vector[PROJ]), aes(x=treatment_year, y=FEve_mean, col=treatment,
                                                                                 ymin=FEve_mean-FEve_se, ymax=FEve_mean+FEve_se)) +
    geom_errorbar(width=0.1) +
    geom_point() +
    geom_path() +
    ggtitle(site_proj_comm_vector[PROJ]) +
    theme_bw()
  ggsave(filename=paste0("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Working groups\\sDiv\\Dec2021\\figures\\FEve\\",
                         site_proj_comm_vector[PROJ], "FEve.jpg"))
}

# F div versus species number
for(PROJ in 1:length(site_proj_comm_vector)){
  ggplot(filter(FD_df_master, site_proj_comm == site_proj_comm_vector[PROJ]), aes(x=nbsp, y=FDis, col=treatment)) +
    geom_point() +
    ggtitle(site_proj_comm_vector[PROJ]) +
    theme_bw()
  ggsave(filename=paste0("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Working groups\\sDiv\\Dec2021\\figures\\FDis\\",
                         site_proj_comm_vector[PROJ], "FDis.jpg"))
}

ggplot(FD_df_master, aes(x=nbsp, y=FDis, col=site_proj_comm)) +
#  geom_point() +
  geom_smooth(method="lm", se=T) +
  theme_bw() +
  theme(legend.position = "none")

################## End of current script ##########################




