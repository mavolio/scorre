# sCoRRe - Winners losers community assembly 
# December 14, 2022 
# Authors: Magda Garbowski 

# Goals: 
# get averages of gower distances by species x site_proj_com_trt
# plot distances ~ diff 

library(ggplot2)

my.wd <- "C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Working groups\\sDiv\\"
  
dat_focal_all_sps_df <- read.csv("/Users/MagdaGarbowski 1/scorre/community assembly winnerslosers/generated_data/distances_avgs.csv")
CT_diff <- read.csv("/Users/MagdaGarbowski 1/scorre/community assembly winnerslosers/generated_data/CT_diff.csv")
# reading Kevin's files in
dat_focal_all_sps_df <- read.csv(paste0(my.wd, "distances_avgs.csv"))
CT_diff <- read.csv(paste0(my.wd, "CT_diff_2022Dec15.csv")) %>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_"))

# get plot information from CT diff
CT_diff_ss <- CT_diff[c("site_code", "project_name", "treatment", "treatDCi", "DCi", "diff", "trt_type", "species_matched", "site_proj_comm")]
#CT_diff_ss$site_proj_comm <- gsub(" ", "_", CT_diff_ss$site_proj_comm)


dat_focal_all_sps_df_merged <- merge(dat_focal_all_sps_df, CT_diff_ss, 
                                     by = c("site_code", "project_name", "treatment", "species_matched", "site_proj_comm"), all.x = TRUE)

# selected out "full" dataset 
dat_focal_all <- dat_focal_all_sps_df_merged[dat_focal_all_sps_df_merged$data_set == "full_all",]

# get different treatment datasets  
dat_focal_all_N <- dat_focal_all[dat_focal_all$trt_type %in% c("N"),] # nitrogen 
dat_focal_all_mult <- dat_focal_all[dat_focal_all$trt_type %in% c("mult_nutrient"),] # multiple nutrients
dat_focal_all_drought <- dat_focal_all[dat_focal_all$trt_type %in% c("drought"),] # drought 

dat_prep <- function(df){
  dat_ls <- split(df, list(df$site_proj_comm))
  # keep list elements with length > 4 
  dat_ls <- dat_ls[sapply(dat_ls, nrow)>4]
  dat_out <- do.call(rbind, dat_ls)
  return(dat_out)
}

ls_out <- lapply(list(dat_focal_all_N = dat_focal_all_N, 
                      dat_focal_all_mult = dat_focal_all_mult, 
                      dat_focal_all_drought = dat_focal_all_drought), dat_prep)

list2env(ls_out, .GlobalEnv)

plot_function <- function(df, distance, y_var, y_var_title, title){
  plot_out <- ggplot(df, aes(x = distance, y = y_var))+
    geom_point() + 
    geom_smooth(method = "lm") + 
    facet_wrap(~site_code) + labs(title = title, y = y_var_title)
  return(plot_out)
}

nitrogen_all_diff <- plot_function(dat_focal_all_N, dat_focal_all_N$avg_distance_all, dat_focal_all_N$diff, "diff", "Nitrogen TRT - All Comparisons")
nitrogen_NN_diff <- plot_function(dat_focal_all_N, dat_focal_all_N$avg_distance_NN,dat_focal_all_N$diff, "diff", "Nitrogen TRT - Nearest neighbor")

mult_all_diff <- plot_function(dat_focal_all_mult, dat_focal_all_mult$avg_distance_all, dat_focal_all_mult$diff, "diff", "Mult Nutrient TRT - All Comparisons")
mult_NN_diff <- plot_function(dat_focal_all_mult, dat_focal_all_mult$avg_distance_NN, dat_focal_all_mult$diff, "diff", "Mult Nutrient TRT - Nearest neighbor")

drought_all_diff <- plot_function(dat_focal_all_drought, dat_focal_all_drought$avg_distance_all, dat_focal_all_drought$diff, "diff", "Drought TRT - All Comparisons")
drought_NN_diff <- plot_function(dat_focal_all_drought, dat_focal_all_drought$avg_distance_NN, dat_focal_all_drought$diff, "diff", "Drought  TRT - Nearest neighbor")

# plot absolute values of delta Dci and color by species catagory
nitrogen_all_absdiff <- plot_function(dat_focal_all_N, dat_focal_all_N$avg_distance_all, abs(dat_focal_all_N$diff), "|diff|", "Nitrogen TRT - All Comparisons")
nitrogen_NN_absdiff <- plot_function(dat_focal_all_N, dat_focal_all_N$avg_distance_NN,abs(dat_focal_all_N$diff), "|diff|", "Nitrogen TRT - Nearest neighbor")

mult_all_absdiff <- plot_function(dat_focal_all_mult, dat_focal_all_mult$avg_distance_all, abs(dat_focal_all_mult$diff), "|diff|", "Mult Nutrient TRT - All Comparisons")
mult_NN_absdiff <- plot_function(dat_focal_all_mult, dat_focal_all_mult$avg_distance_NN, abs(dat_focal_all_mult$diff), "|diff|", "Mult Nutrient TRT - Nearest neighbor")

drought_all_absdiff <- plot_function(dat_focal_all_drought, dat_focal_all_drought$avg_distance_all, abs(dat_focal_all_drought$diff), "|diff|", "Drought TRT - All Comparisons")
drought_NN_absdiff <- plot_function(dat_focal_all_drought, dat_focal_all_drought$avg_distance_NN, abs(dat_focal_all_drought$diff), "|diff|", "Drought  TRT - Nearest neighbor")


# check a few plots with DCi 

nitrogen_all_DCi <- plot_function(dat_focal_all_N, dat_focal_all_N$avg_distance_all, dat_focal_all_N$DCi, "DCi", "Nitrogen TRT - All Comparisons")
nitrogen_NN_DCi <- plot_function(dat_focal_all_N, dat_focal_all_N$avg_distance_NN,dat_focal_all_N$DCi, "DCi", "Nitrogen TRT - Nearest neighbor")

