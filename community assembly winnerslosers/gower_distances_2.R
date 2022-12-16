# sCoRRe - Winners losers community assembly 
# December 16, 2022 
# Authors: Magda Garbowski

# goals: datasets: 

# (1) control_distances_all: dataset with every species compared to every other species in control plots 
# include distance of colonizers to control plots 
# do this for presence/absence distances and abundance weighted distances 

# (2) treatment_distances_all: dataset with every species compared to every other species in treatment plots 
# colonizers are just part of the community
# do this for presence/absence distances and abundance weighted distances

# (3) comparisons_dataset: dataset with w_w, w_c, w_l, w_n, c_c, c_l, c_n, c_l, l_n, n_n
# (3a) do this with distances from control plots, both presence/abs and abundance weighted 
# (3b) do this with ditances from treatment plots,both presence/abs and abundance weighted 

setwd("/Users/MagdaGarbowski 1/scorre/community assembly winnerslosers/")
# setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Working groups\\sDiv\\")
library(tidyverse)
library(funrar)

#----------------------------- datasets  ---------------------------
rel_cov <- read.csv("/Users/MagdaGarbowski 1/scorre/community assembly winnerslosers/generated_data/relcov.csv")
traitsClean <- read.csv("/Users/MagdaGarbowski 1/scorre/community assembly winnerslosers/generated_data/traitsClean.csv")
diff_quantiles <- read.csv("/Users/MagdaGarbowski 1/scorre/community assembly winnerslosers/generated_data/diff_quantile.csv")
CT_diff <- read.csv("/Users/MagdaGarbowski 1/scorre/community assembly winnerslosers/generated_data/CT_diff.csv")

# merge datasets to get winners and losers 
CT_diff_quantiles <- merge(CT_diff, diff_quantiles, by = c("site_proj_comm", "treatment"), all.x= TRUE)

# get species status based on 90 and 10 quantiles 

win_lose_function <- function(df, quant_low, quant_high, quantile){
  df$status <- ifelse(df$DCi==0 & df$treatDCi>0, "colonizer",
                                     ifelse(df$diff < df[[quant_low]], "loser",
                                            ifelse(df$diff > df[[quant_high]], "winner","neutral")))
  df$quant <- quantile
  return(df)
}

CT_diff_quantiles_out <- do.call(rbind, list(win_lose_function(CT_diff_quantiles, "p05", "p95", "p05_p95"),
                                             win_lose_function(CT_diff_quantiles, "p10", "p90", "p10_p90"),
                                             win_lose_function(CT_diff_quantiles, "p20", "p80", "p20_p80")))


# select out quantiles and columns of interest 
CT_diff_quantiles_ss <- 
  CT_diff_quantiles_out[CT_diff_quantiles_out$quant == "p05_p95",][c("site_proj_comm", "treatment", "site_code", "project_name", "species_matched", "status", "quant")]
                          
# add "_" to site_proj_comm for matching 
CT_diff_quantiles_ss$site_proj_comm <- gsub(" ", "_", CT_diff_quantiles_ss$site_proj_comm)

# ---------------- control distances all along with colonizers  ---------------- 
# need plot_level species list (controls)
# need trt_site_comm level species list for colonizers 

relcov_full_clean_control <- relcov_full_clean[relcov_full_clean$treatment == "control",]

# split by site_proj_comm, plot_id, and year 
relcov_full_clean_control_splits <- split(relcov_full_clean_control , 
                                          list(relcov_full_clean_control$site_proj_comm,
                                               relcov_full_clean_control$plot_id, 
                                               relcov_full_clean_control$calendar_year),
                                          drop = TRUE)

# keep plots with >1 species  
relcov_full_clean_control_splits <- relcov_full_clean_control_splits[sapply(relcov_full_clean_control_splits, nrow)>1] 

#----------------------------- calculating distances ---------------------------

# all distances of species to one another in given plot (so far in controls)

all_distances_function <- function(plot_df){
  
  # species vector from relative cover for pulling traits 
  sp_df_temp <- unique(plot_df$species_matched)[!is.na(unique(plot_df$species_matched))] 
  
  # subset trait data to just include the species present subset relative cover data
  traits_df_raw_temp <- traitsClean[traitsClean$species_matched %in% c(sp_df_temp),]
  
  if(nrow(traits_df_raw_temp) > 2){
  
  # vector of species for which trait data are available 
  traits_df_raw_temp_vec <- unique(traits_df_raw_temp$species_matched)
  
  # vector of species not in trait database (but in relative abundance data) to remove from species abundance data
  sp_to_remove_temp <- plot_df[is.na(plot_df$species_matched),]$genus_species
  
  # abundance dataset with species removed that do not have trait information
  relcov_unkn_sp_rm_temp <- plot_df[!plot_df$genus_species %in% c(sp_to_remove_temp),]
  
  # dataframe for categorical and continuous traits
  traits_df_temp_cat_cont <- traits_df_raw_temp[c("species_matched", "photosynthetic_pathway","lifespan",  "mycorrhizal_type", "n_fixation","clonal",
                                                  "seed_dry_mass","LDMC","plant_height_vegetative", "SLA", "rooting_depth")]
  # make categorical traits factors
  traits_df_temp_cat_cont[,c(2:6)] <- lapply(traits_df_temp_cat_cont[,c(2:6)], as.factor)
  
  # get data read for compute_dist_matrix 
  rownames(traits_df_temp_cat_cont) <- traits_df_temp_cat_cont$species_matched
  traits_df_temp_cat_cont$species_matched <- NULL
  
  # calculate gower distance 
  gower_dist <- as.matrix(compute_dist_matrix(traits_df_temp_cat_cont, "gower"))
  
  # create output dataset 
  out <- as.data.frame(matrix(gower_dist, dimnames=list(t(outer(colnames(gower_dist), rownames(gower_dist), FUN=paste, sep = "_")), NULL)))
  out$species_1 <- sub("_.*", "", rownames(out)) 
  out$species_2 <- sub(".*_", "", rownames(out)) 
  
  colnames(out)[1] <- "distance"
  rownames(out) <- NULL
  
  out <- out[c("species_1", "species_2", "distance")]
  
  # cbind with plot info and export 
  info_dat <- plot_df[1,c(1:11)][rep(seq_len(nrow(plot_df[1,c(1:11)])), nrow(out)),]
  out_3 <- cbind(info_dat, out)
  return(out_3)
  }
}

out <- lapply(relcov_full_clean_control_splits, all_distances_function)

# These datasets can now be used to get specific averages? 


  
  # create subset datasets for comparisons 
  # full dataset
  focal_all <- gower_w_winners_losers[,c(3:length(gower_w_winners_losers))]
  focal_all <- focal_all[c(rownames(focal_all))]
  
  if(any(gower_w_winners_losers == "winner")){
    win_vec <- gower_w_winners_losers[gower_w_winners_losers$species_status_90_10 == "winner",]$species_matched
    win_cols <- gower_w_winners_losers[c("species_matched", "species_status_90_10", win_vec)]
    
    # winner to winner
    win_winner <- win_cols[win_cols$species_status_90_10 %in% "winner",]
    winner_winner <- as.data.frame(win_winner[,c(3:length(win_winner))])
    colnames(winner_winner) <- win_vec
    
    # winner to neutral 
    win_neutral <- win_cols[win_cols$species_status_90_10 %in% "neutral",]
    winner_neutral <- as.data.frame(win_neutral[,c(3:length(win_neutral))])
    colnames(winner_neutral) <- win_vec
    
    # winner to loser 
    win_loser <- win_cols[win_cols$species_status_90_10 %in% "loser",]
    winner_loser <- as.data.frame(win_loser[,c(3:length(win_loser))])
    colnames(winner_loser) <- win_vec
    
    # put datasets into a list 
    dist_mat_ls <- list(focal_all = focal_all, 
                        winner_winner = winner_winner,
                        winner_neutral = winner_neutral,
                        winner_loser = winner_loser)
  }
  
  else{
    dist_mat_ls <- list(focal_all = focal_all)
  }
  
  # keep list elements with length > 0 
  dist_mat_ls_2 <- dist_mat_ls[sapply(dist_mat_ls, nrow)>1]
  
  # function to get distances of different datasets 
  distance_avg_function <- function(dis_mat){
    # make all 0s NA 
    dis_mat[dis_mat == 0] <- NA 
    df = dis_mat 
    
    avg_fun <- function(df){
      
      out_avg <- data.frame(species = colnames(df),
                            distance = colSums(df, na.rm = TRUE)/(colSums(!is.na(df[1]))),
                            NN = apply(df,2,min, na.rm = TRUE))
      return(out_avg)
    }
    
    plot_info_function <- function(plot_df){
      out_plot_info <- cbind(relcov_unkn_sp_rm_temp[1,c(1:4,6:11)], plot_df)
      return(out_plot_info)
    }
    
    out <- avg_fun(df)
    out_2 <- plot_info_function(out)
    return(out_2)
  } 
  out_all <- do.call(rbind, lapply(dist_mat_ls_2, distance_avg_function))
  
  out_all$comparison <- gsub("\\..*", "", rownames(out_all))
  return(out_all)
}


# get distances for all species at all sites 
# this takes about 10 minutes 
plots_out_splits_2 <- plots_out_splits[-c(1388, 4752,7205,7248,7542,7715,7847,7879,7885,8264,8750,8753,10002)]
sites_distances_out_comparisons <- do.call(rbind, lapply(plots_out_splits_2, distances_function)) 

#----------------------------- calculating averages ---------------------------

# select out "focal_all" dataset - can later work with winner_winner or other subset 
dat_focal_all <- sites_distances_out_comparisons[sites_distances_out_comparisons$comparison == "focal_all",]
# split by site_proj_comm, species, and treatment for averages 
dat_focal_all_sps_splits <- split(dat_focal_all, 
                                  list(dat_focal_all$site_proj_comm, dat_focal_all$species, dat_focal_all$treatment), drop = TRUE)

# winner winner 
dat_winner_winner <- sites_distances_out_comparisons[sites_distances_out_comparisons$comparison == "winner_winner",]
# split by site_proj_comm, species, and treatment for averages 
dat_winner_winner_sps_splits <- split(dat_focal_winner_winner, 
                                            list(dat_winner_winner$site_proj_comm, dat_winner_winner$species, dat_winner_winner$treatment), drop = TRUE)

# winner loser
dat_winner_loser <- sites_distances_out_comparisons[sites_distances_out_comparisons$comparison == "winner_loser",]
# split by site_proj_comm, species, and treatment for averages 
dat_winner_loser_sps_splits <- split(dat_winner_loser, 
                                            list(dat_winner_loser$site_proj_comm, dat_winner_loser$species, dat_winner_loser$treatment), drop = TRUE)

# winner neutral
dat_winner_neutral <- sites_distances_out_comparisons[sites_distances_out_comparisons$comparison == "winner_neutral",]
# split by site_proj_comm, species, and treatment for averages 
dat_winner_neutral_sps_splits <- split(dat_winner_neutral, 
                                     list(dat_winner_neutral$site_proj_comm, dat_winner_neutral$species, dat_winner_neutral$treatment), drop = TRUE)

# get averages
avg_function <- function(df, datset){
  out <-  data.frame(site_proj_comm = df$site_proj_comm[1],
                     treatment = df$treatment[1], 
                     site_code = df$site_code[1],
                     project_name = df$project_name[1],
                     data_type = df$data_type[1],
                     species_matched = df$species[1],
                     avg_distance_all = mean(df$distance, na.rm = TRUE), 
                     avg_distance_NN = mean(df$NN, na.rm = TRUE),
                     sd_distance_all = sd(df$distance, na.rm = TRUE), 
                     sd_distance_NN = sd(df$NN, na.rm = TRUE), 
                     n_obs = nrow(df),
                     data_set = datset)
  return(out)
}

dat_focal_all_sps_df <- do.call(rbind, lapply(dat_focal_all_sps_splits, avg_function, "full_all"))
dat_winner_winner_sps_df <- do.call(rbind, lapply(dat_winner_winner_sps_splits, avg_function, "winner_winner"))
dat_winner_loser_sps_df <- do.call(rbind, lapply(dat_winner_loser_sps_splits, avg_function, "winner_loser"))
dat_winner_neutral_sps_df <- do.call(rbind, lapply(dat_winner_neutral_sps_splits, avg_function, "winner_neutral"))

dat_avgs <- do.call(rbind, list(dat_focal_all_sps_df, dat_winner_winner_sps_df, dat_winner_loser_sps_df, dat_winner_neutral_sps_df ))


write.csv(sites_distances_out, "/Users/MagdaGarbowski 1/scorre/community assembly winnerslosers/distances.csv")
write.csv(sites_distances_out_comparisons, "/Users/MagdaGarbowski 1/scorre/community assembly winnerslosers/generated_data/distances_winner_comparisons.csv")
write.csv(dat_avgs, "/Users/MagdaGarbowski 1/scorre/community assembly winnerslosers/generated_data/distances_avgs.csv")



# ---------------------------------
# to delete 
relcov_full_clean_win_los <- merge(relcov_full_clean, CT_diff_quantiles_ss, 
                                   by = c("site_proj_comm", "treatment", "site_code", "project_name", "species_matched"),
                                   all.x = TRUE)

# pull out only sites/treatments with winners, losers, colonizers, neutrals identified 
relcov_full_clean_win_los_splits <- split(relcov_full_clean_win_los, 
                                          list(relcov_full_clean_win_los$site_proj_comm, relcov_full_clean_win_los$treatment),
                                          drop = TRUE)

# keep only datasets with some winners and losers in them 
return_df_function <- function(x){
  if(sum(is.na(x$species_status_90_10)) != nrow(x))
    return(x)
}


#----------------- merge winners and losers with relcov_full_clean CONTROL only  -------------------------

relcov_full_clean_control <- relcov_full_clean[relcov_full_clean$treatment == "control",]

relcov_full_clean_win_los_control <- merge(relcov_full_clean_control, CT_diff_quantiles_ss[], 
                                           by = c("site_proj_comm", "site_code", "project_name", "species_matched"),
                                           all.x = TRUE)

# pull out only sites/treatments with winners, losers, colonizers, neutrals identified 
relcov_full_clean_win_los_control_splits <- split(relcov_full_clean_win_los_control , 
                                                  list(relcov_full_clean_win_los_control$site_proj_comm, relcov_full_clean_win_los_control$treatment),
                                                  drop = TRUE)
prep_function <- function(plot_df){
  # split to plot level 
  plot_splits <- split(plot_df, list(plot_df$plot_id, plot_df$calendar_year), drop = TRUE)
  
  # keep communities with more than one species 
  plot_splits_2 <- plot_splits[sapply(plot_splits, nrow)>1] 
  
  # drop NA from winner_loser column 
  plot_splits_2 <- lapply(plot_splits_2, function(x) {xx = x[!is.na(x$species_status_90_10),]; return(xx)})
  
  # keep plots with winners 
  winners_keep_function <- function(df){
    if(any(df == "winner")){
      return(df)
    }
  }
  
  # keep plots with more than one species 
  plot_splits_w_winners <- lapply(plot_splits_2, winners_keep_function)
  plot_splits_w_winners[sapply(plot_splits_w_winners, is.null)] <- NULL
  plot_splits_w_winners <- plot_splits_w_winners[sapply(plot_splits_w_winners, nrow)>1]
  
  if(length(plot_splits_w_winners) > 0){
    plots_out <- do.call(rbind, plot_splits_w_winners)
    return(plots_out)
  }
}

plots_out_ls <- lapply(site_splits, prep_function)
plots_out_ls[sapply(plots_out_ls, is.null)] <- NULL

plots_out_df <- do.call(rbind, plots_out_ls)

plots_out_splits <- split(plots_out_df, list(plots_out_df$site_proj_comm, plots_out_df$treatment, plots_out_df$calendar_year, plots_out_df$plot_id), drop = TRUE)









# winner loser dataframe 
win_los_df <- plot_df[c("species_matched", "species_status_90_10")]
win_los_df <- win_los_df[win_los_df$species_matched %in% traits_df_raw_temp_vec,]

# set 1 to get gower distances out
set_1 <- traits_df_temp_cat_cont
rownames(set_1) <- set_1$species_matched
set_1$species_matched <- NULL 

# calculate gower distance 
gower_dist <- as.data.frame(compute_dist_matrix(set_1, "gower"))
gower_dist$species_matched <- rownames(gower_dist)

# get into same order for cbind 
gower_dist <- gower_dist[order(gower_dist$species_matched),]
win_los_df <- win_los_df[order(win_los_df$species_matched),]
win_los_df <- unique(win_los_df)

gower_w_winners_losers <- cbind(win_los_df,gower_dist[,!names(gower_dist) %in% c("species_matched")])
row.names(gower_w_winners_losers) <- gower_w_winners_losers$species_matched

