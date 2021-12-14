# sCoRRE - Covariance in community-level traits 
# Key question: Do environmental change drivers increase covariance in community traits?
# Based on Dwyer and Laughlin 2017, Ecology Letters 
# 

library(ggplot2)
library(gridExtra)

ab_data <- read.csv("/Users/MagdaGarbowski 1/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_RawAbundance_Dec2021.csv")
trait_data_cont <- read.csv("/Users/MagdaGarbowski 1/Dropbox/sDiv_sCoRRE_shared/Trait Data/TRY Data/TRY Continuous data/TRY_trait_data_continuous_Nov2021.csv")
trait_data_cont$species_matched <- tolower(trait_data_cont$species_matched)

ab_data_nutnet <- ab_data[(ab_data$project_name %in% c("NutNet") & ab_data$treatment %in% c("Control", "NPK")),]
nutnet_sp_list <- unique(ab_data_nutnet[c("genus_species")])
nutnet_traits <- trait_data_cont[trait_data_cont$species_matched %in% c(nutnet_sp_list$genus_species),]    

ab_data_edge <- ab_data[(ab_data$project_name %in% c("EDGE") & ab_data$treatment %in% c("con", "int")),]
edge_sp_list <- unique(ab_data_edge[c("genus_species")])
edge_traits <- trait_data_cont[trait_data_cont$species_matched %in% c(edge_sp_list$genus_species),]
  
# select out traits to work with 
# 4 - seed dry mass
# 50 - leaf nitrogen content 
# 3106 - vegetative height
# 3107 - generative height 
# 3115, 3116, 3117 - SLA 


########### Traits ####################
nutnet_traits_ss <- nutnet_traits[c("species_matched", "leaf_C", "leaf_N", "SLA", "seed_dry_mass", "plant_height_generative", "plant_height_vegetative")]
edge_traits_ss <- edge_traits[c("species_matched", "leaf_C", "leaf_N", "SLA", "seed_dry_mass", "plant_height_generative", "plant_height_vegetative")]

df_traits_avg <- function(df_traits_ss){
  aggregate(list(leaf_C = df_traits_ss$leaf_C, 
                 leaf_N = df_traits_ss$leaf_N, 
                 SLA = df_traits_ss$SLA, 
                 seed_dry_mass = df_traits_ss$seed_dry_mass, 
                 plant_height_generative = df_traits_ss$plant_height_generative, 
                 plant_height_vegetative = df_traits_ss$plant_height_vegetative), 
            by = list(species_matched = df_traits_ss$species_matched), FUN = mean, na.rm = TRUE)
  }

nutnet_traits_avg <- df_traits_avg(nutnet_traits_ss)
edge_traits_avg <- df_traits_avg(edge_traits_ss)

# get single column for height 
nutnet_traits_avg$height_comb <- ifelse(is.na(nutnet_traits_avg$plant_height_vegetative), nutnet_traits_avg$plant_height_generative, nutnet_traits_avg$plant_height_vegetative)
edge_traits_avg$height_comb <- ifelse(is.na(edge_traits_avg$plant_height_vegetative), edge_traits_avg$plant_height_generative, edge_traits_avg$plant_height_vegetative)

# nutnet 
nutnet_comm_traits <- merge(ab_data_nutnet, nutnet_traits_avg, by.x = "genus_species", by.y = "species_matched", all.x = TRUE)
nutnet_splits <- split(nutnet_comm_traits, nutnet_comm_traits$site_code)

nutnet_comm_traits_lastyr <-  lapply(nutnet_splits, function (x) {
  x = x; x$yr_group = ifelse(x$treatment_year == max(x$treatment_year),"last", "other");
  y = x[x$yr_group =="last",]; return (y)})

#edge 
edge_comm_traits <- merge(ab_data_edge, edge_traits_avg, by.x = "genus_species", by.y = "species_matched", all.x = TRUE)
edge_splits <- split(edge_comm_traits, edge_comm_traits$site_code)

edge_comm_traits_lastyr <-  lapply(edge_splits, function (x) {
  x = x; x$yr_group = ifelse(x$treatment_year == max(x$treatment_year),"last", "other");
  y = x[x$yr_group =="last",]; return (y)})


sps_function <- function(df){
  df_comm_traits_lastyr_sp <- unique(df[c("site_code", "genus_species", "treatment","leaf_C", "leaf_N", "SLA", "seed_dry_mass", "height_comb")])
  return(df_comm_traits_lastyr_sp)
}

site_sps_out_nutnet <- lapply(nutnet_comm_traits_lastyr, sps_function)
site_sps_out_nutnet_all <- do.call(rbind, site_sps_out_nutnet)
site_sps_out_edge <- lapply(edge_comm_traits_lastyr, sps_function)

plot_function <- function (df, x_var, y_var, control, trt) {
  df_ss = df[c("site_code", "genus_species", "treatment", x_var, y_var)]
  df_ss = df_ss[complete.cases(df_ss),]
  df_ss_control = df_ss[df_ss$treatment == control,]
  df_ss_trt = df_ss[df_ss$treatment == trt,]
  cov_cont = round(cov(df_ss_control[[x_var]], df_ss_control[[y_var]]), digits = 2)
  cov_trt = round(cov(df_ss_trt[[x_var]], df_ss_trt[[y_var]]), digits = 2)
  corr_cont = round(cor(df_ss_control[[x_var]], df_ss_control[[y_var]]), digits = 2)
  corr_trt = round(cor(df_ss_trt[[x_var]], df_ss_trt[[y_var]]), digits = 2)
  
  text_data <- data.frame(
    label = c((paste("cov:", cov_cont, "corr:", corr_cont)), c(paste("cov:",cov_trt, "corr:", corr_trt))),
    treatment = c(control, trt))
  
  plot <- ggplot(df_ss, aes_string(x = x_var, y = y_var)) +
    geom_point() +
    geom_smooth(method = "lm") + 
    labs(title = df$site_code[1])+
    facet_wrap(~ treatment)
  
  plot2 <- plot + geom_text(
    data = text_data, 
    mapping = aes(x = -Inf, y = -Inf, label = label), 
    size = 5,
    hjust = -0.1, 
    vjust = -7)
  return(plot2)
}


nutnet_all <- plot_function(site_sps_out_nutnet_all,"height_comb",  "seed_dry_mass", "Control", "NPK")



# Plots and covariances 
do.call(grid.arrange, c(lapply(site_sps_out_nutnet, plot_function,"height_comb",  "seed_dry_mass", "Control", "NPK")))
do.call(grid.arrange, c(lapply(site_sps_out_nutnet, plot_function,"height_comb",  "leaf_N", "Control", "NPK")))
do.call(grid.arrange, c(lapply(site_sps_out_nutnet, plot_function,"height_comb",  "SLA", "Control", "NPK")))
do.call(grid.arrange, c(lapply(site_sps_out_nutnet, plot_function,"seed_dry_mass",  "SLA", "Control", "NPK")))
do.call(grid.arrange, c(lapply(site_sps_out_nutnet, plot_function,"seed_dry_mass",  "leaf_N", "Control", "NPK")))
do.call(grid.arrange, c(lapply(site_sps_out_nutnet, plot_function,"SLA",  "leaf_N", "Control", "NPK")))

# CDR plots 
# Within the LES - SLA & leaf_N 

grid.arrange(plot_function(site_sps_out_nutnet$CDR, "SLA",  "leaf_N","Control", "NPK"), 
             plot_function(site_sps_out_nutnet$CDR, "height_comb", "seed_dry_mass", "Control", "NPK"),
             plot_function(site_sps_out_nutnet$CDR, "height_comb", "SLA", "Control", "NPK"))

grid.arrange(plot_function(site_sps_out_edge$KNZ, "SLA",  "leaf_N","con", "int"), 
             plot_function(site_sps_out_edge$KNZ, "height_comb", "seed_dry_mass","con", "int"),
             plot_function(site_sps_out_edge$KNZ, "height_comb", "SLA","con", "int"))

# To do : 
# (1) Work with CDR site to show how covariance might change 
# (2) Pull out covariances for the trait combinations and plot to see if there are consistent differences

cov_var_function <- function(df, x_var, y_var){
  df_ss = df[c("site_code", "genus_species", "treatment", x_var, y_var)]
  df_ss = df_ss[complete.cases(df_ss),]
  df_ss_control = df_ss[df_ss$treatment == "Control",]
  df_ss_NPK = df_ss[df_ss$treatment == "NPK",]
  cov_cont = round(cov(df_ss_control[[x_var]], df_ss_control[[y_var]]), digits = 2)
  cov_NPK = round(cov(df_ss_NPK[[x_var]], df_ss_NPK[[y_var]]), digits = 2)
  df_out = data.frame(cov = c(cov_cont, cov_NPK))
  df_out$cov_abs = abs(df_out$cov)
  df_out$treatment = c("Control", "NPK")
  df_out$site_code = df_ss$site_code[1]
  df_out$combo = paste(x_var, y_var, sep = "_vs_")
  return(df_out)
}

cov_height_seed <-do.call(rbind, lapply(site_sps_out, cov_var_function, "height_comb",  "seed_dry_mass"))
cov_height_leafN <-do.call(rbind, lapply(site_sps_out, cov_var_function, "height_comb",  "leaf_N"))
cov_height_SLA <-do.call(rbind, lapply(site_sps_out, cov_var_function, "height_comb",  "SLA"))


ggplot(cov_height_seed, aes(x = treatment, y = cov)) + 
  geom_boxplot() + geom_point(aes(color = site_code), position = position_dodge(0.1))

ggplot(cov_height_seed, aes(x = treatment, y = cov_abs)) + 
  geom_boxplot()+ geom_point(aes(color = site_code), position = position_dodge(0.1))

ggplot(cov_height_leafN, aes(x = treatment, y = cov)) + 
  geom_boxplot() + geom_point(aes(color = site_code), position = position_dodge(0.1))

ggplot(cov_height_leafN, aes(x = treatment, y = cov_abs)) + 
  geom_boxplot()+ geom_point(aes(color = site_code), position = position_dodge(0.1))

ggplot(cov_height_SLA, aes(x = treatment, y = cov)) + 
  geom_boxplot() + geom_point(aes(color = site_code), position = position_dodge(0.1))

ggplot(cov_height_SLA, aes(x = treatment, y = cov_abs)) + 
  geom_boxplot()+ geom_point(aes(color = site_code), position = position_dodge(0.1))
