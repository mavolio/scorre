
my.wd <- "~/Dropbox/sCoRRE/sDiv_sCoRRE_shared"



# only experimennts that are at least 5 years

##############################################
# identify winners and losers per experiment #
##############################################
# DCi though time for control and treatment and then slopes in Adam's framework + extinct and immigration species
# important!remove too rare species; maybe identify via CSi values, do histograms to identify threshlds but it could be around 0.1


# Emily developed code to generate DCi through time: “emily ambient v2.R”







##############################################
# calculate mdns, dnns, wmdns
##############################################

# first step in one experiment
ab_data_nutnet <- ab_data[(ab_data$project_name %in% c("NutNet") & ab_data$treatment %in% c("Control", "NPK")),]
nutnet_sp_list <- unique(ab_data_nutnet[c("genus_species")])
nutnet_traits <- trait_data_cont[trait_data_cont$species_matched %in% c(nutnet_sp_list$genus_species),]    

nutnet_traits_ss <- nutnet_traits[c("species_matched", "LDMC", "leaf_N", "SLA", "seed_dry_mass", "plant_height_generative",
                                    "plant_height_vegetative", "root_density", "rooting_depth")]

# take average trait values 
df_traits_avg <- function(df_traits_ss){
  aggregate(list(LDMC = df_traits_ss$LDMC, 
                 leaf_N = df_traits_ss$leaf_N, 
                 SLA = df_traits_ss$SLA, 
                 seed_dry_mass = df_traits_ss$seed_dry_mass, 
                 plant_height_generative = df_traits_ss$plant_height_generative, 
                 plant_height_vegetative = df_traits_ss$plant_height_vegetative,
                 root_density = df_traits_ss$root_density, 
                 rooting_depth = df_traits_ss$rooting_depth), 
            by = list(species_matched = df_traits_ss$species_matched), FUN = mean, na.rm = TRUE)
}

nutnet_traits_avg <- df_traits_avg(nutnet_traits_ss)

# merge with abundance values 
nutnet_comm_traits <- merge(ab_data_nutnet, nutnet_traits_avg, by.x = "genus_species", by.y = "species_matched", all.x = TRUE)













#