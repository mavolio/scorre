# functional dissimilarity - winners and losers 
# code up until line 50 is from 3_sCoRRE_functionalDiversityMetrics.R" 

setwd("/Users/MagdaGarbowski 1/Dropbox/sDiv_sCoRRE_shared/")
library(tidyverse)
library(funrar)

# continuous trait data
contTraits <- read.csv("CoRRE data/trait data/Final TRY Traits/Imputed Continuous_Traits/data to play with/imputed_continuous_20220620.csv")%>%
  select(-X.1, -X, -family, -genus, -observation)%>%
  group_by(species_matched)%>%
  summarise_all(funs(mean))%>%
  ungroup()%>%
  mutate_at(vars(seed_dry_mass:seed_number), scale) #scale continuous traits

traitsAll <- read.csv("CoRRE data/trait data/sCoRRE categorical trait data_11302021.csv")%>% #categorical trait data
  full_join(contTraits) %>% #merge continuous and categorical traits
  drop_na()

#remove non-target species (ferns, lycophytes, mosses) that somehow snuck into the trait data
traitsClean <- traitsAll %>%
  filter(!leaf_type %in% c("microphyll","frond")) %>%
  filter(!species_matched %in% c("Centrolepis aristata", "Centrolepis strigosa", "Acorus calamus"))

#test - are we missing any categorical data?
sp_without_catag_data <- filter(traitsAll, is.na(clonal)) #none, we did a great job collecting categorical data!

# species relative cover data
relcov_full_raw <- read.csv("CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Dec2021.csv") %>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>%
  dplyr::select(site_code:community_type, site_proj_comm, calendar_year:relcov)

# corre to try species names key
corre_to_try <- read.csv("CoRRE data/trait data/corre2trykey_2021.csv") %>%
  dplyr::select(genus_species, species_matched) %>%
  unique(.)

### merge species names and remove all mosses
# moss key to remove mosses from species comp data
moss_sp_vec <- read.csv("CoRRE data/trait data/sCoRRE categorical trait data_11302021.csv") %>%
  dplyr::select(species_matched, leaf_type) %>%
  mutate(moss = ifelse(leaf_type=="moss", "moss","non-moss")) %>%
  filter(moss=="moss") %>%
  pull(species_matched)

relcov_full_clean <- relcov_full_raw %>%
  dplyr::left_join(corre_to_try, by="genus_species") %>%
  filter(!species_matched  %in% moss_sp_vec) %>%
  filter(!(site_code == "EGN" & project_name == "Nmow" & community_type == "0" & 
             plot_id == "19" & calendar_year==2015 & treatment=="Control")) ## Stipa species just has an incorrect treatment I think, needs to be fixed in absolute abundance data frame but I'm just removing it for now


#----------------------------- calculating distances ---------------------------

# split into site_proj_comm-level datasets 
site_splits <- split(relcov_full_clean, list(relcov_full_clean$site_proj_comm))

distances_function <- function(site_df){
  # split to plot level 
  plot_splits <- split(site_df, list(site_df$plot_id, site_df$calendar_year), drop = TRUE)
  distances_function_2 <- function(plot_df){
    # species vector from relative cover for pulling traits 
    sp_df_temp <- unique(plot_df$species_matched)[!is.na(unique(plot_df$species_matched))] 
    
    #subset trait data to just include the species present subset relative cover data
    traits_df_raw_temp <- traitsClean[traitsClean$species_matched %in% c(sp_df_temp),]
    
    #vector of species for which trait data are available 
    traits_df_raw_temp_vec <- unique(traits_df_raw_temp$species_matched)
    
    #vector of species not in trait database (but in relative abundance data) to remove from species abundance data
    sp_to_remove_temp <- plot_df[is.na(plot_df$species_matched),]$genus_species
    
    #abundance dataset with species removed that do not have trait information
    relcov_unkn_sp_rm_temp <- plot_df[!plot_df$genus_species %in% c(sp_to_remove_temp),]
    
    # dataframe for categorical and continuous traits
    traits_df_temp_cat_cont <- traits_df_raw_temp[c("species_matched", "photosynthetic_pathway","lifespan",  "mycorrhizal_type", "n_fixation","clonal",
                                                    "seed_dry_mass","LDMC","plant_height_vegetative", "SLA", "rooting_depth")]
    
    # make values numeric (continuous traits) and factors (categorical traits)
    # will need to be careful if selected traits change 
    traits_df_temp_cat_cont[,c(7:length(traits_df_temp_cat_cont))] <- lapply(traits_df_temp_cat_cont[,c(7:length(traits_df_temp_cat_cont))],
                                                                              as.numeric)
    traits_df_temp_cat_cont[,c(2:6)] <- lapply(traits_df_temp_cat_cont[,c(2:6)], as.factor)
    
    rownames(traits_df_temp_cat_cont) <- traits_df_temp_cat_cont$species_matched
    
    traits_df_temp_cat_cont$species_matched <- NULL 

    # calculate gower distance 
    gower_dist<- as.data.frame(compute_dist_matrix(traits_df_temp_cat_cont, "gower"))

    dist_mat_ls <- list(gower_dist = gower_dist)
    
    distance_avg_function <- function(dis_mat){
      avg_fun <- function(df){
        out_avg <- data.frame(species = rownames(df),
                              allsps = colSums(df, na.rm = TRUE)/(nrow(df) -1),
                              NN = apply(df,2,min, na.rm = TRUE))
        return(out_avg)
      }
      
      plot_info_function <- function(plot_df){
        plot_df = out
        out_plot_info <- cbind(relcov_unkn_sp_rm_temp[1,(1:11)], plot_df)
        return(out_plot_info)
      }
      
      # make all 0s NA 
      dis_mat[dis_mat == 0] <- NA 
      df = dis_mat 
      
      out <- avg_fun(df)
      out_2 <- plot_info_function(out)
      return(out_2)
    } 
    out_all <- do.call(rbind, lapply(dist_mat_ls, distance_avg_function))
    out_all$metric <- gsub("\\..*", "", rownames(out_all))
    return(out_all)
  }
  plot_splits_out <- do.call(rbind, lapply(plot_splits, distances_function_2))
  rownames(plot_splits_out) <- NULL
  return(plot_splits_out)
}

# drop a few sites from main list because they are causing issues: 21,41,88,93
# need to figure out what is causing issues 

site_splits_2 <- site_splits[-c(21,41,88,93)]

# get distances for all species at all sites 
# this takes about 10 minutes 
sites_distances_out <- do.call(rbind, lapply(site_splits_2, distances_function))
write.csv(sites_distances_out, "/Users/MagdaGarbowski 1/scorre/community assembly winnerslosers/distances.csv")

