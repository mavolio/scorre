#library(tidyverse)
library(tidyr)
library(ggplot2)
library(dplyr)
library(plyr)
library(reshape2)
library(lmerTest)
library(ggpubr)
library(vegan)
library(FD)
library(visreg)
library(ggthemes)
library(codyn)

#Read in data
traits_cat <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/sCoRRE categorical trait data_11302021.csv") #categorical trait data

traits1 <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/Final TRY Traits/Imputed Continuous_Traits/data to play with/imputed_continuous_20220620.csv")
corre2trykey <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/corre2trykey_2021.csv") #contrinuous trait data

cover <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Dec2021.csv") %>% #community comp relative cover data
    mutate(drop=ifelse(site_code=="CDR"&treatment==2|site_code=="CDR"&treatment==3|site_code=="CDR"&treatment==4|site_code=="CDR"&treatment==5|site_code=="CDR"&treatment==7, 1,0))%>%
  filter(drop==0) #remove some Cedar Creek treatments since that site is somewhat overrepresented


corre2trykey <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/corre2trykey_2021.csv") #matched species names between trait data and relative cover data
corre2trykey <- corre2trykey[,c("genus_species","species_matched")]
corre2trykey <- unique(corre2trykey)
cover <- left_join(cover, corre2trykey, by = "genus_species", all.x = TRUE)

experimentinfo <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfo_Dec2021.csv")#Information about the treatments which gets used to test how treatment magnitude explains efect sizes




#create average trait values per species and clean outliers
traits<-traits1%>%
  select(-X.1, -X, -family, -genus, -observation)%>%
  select(species_matched, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass, leaf_C.N) %>%
  group_by(species_matched)%>%
  summarise_all(funs(mean))%>%
  ungroup()  %>%
  filter(seed_dry_mass<30, plant_height_vegetative<10, rooting_depth<3, SLA<75)


#standardize the scale of continuous traits
cols <- c( "seed_dry_mass", 
           # "stem_spec_density",
           #"leaf_N",
           #"leaf_P",
           "LDMC",
           #"leaf_C",
           #"leaf_dry_mass",
           "plant_height_vegetative",
           "leaf_C.N",
           "SLA",
           #"water_content",
           "rooting_depth"#,
           #"SRL",
           #"seed_number"
)

traits[cols] <- scale(traits[cols])

traits <- left_join(traits, traits_cat, by = "species_matched")#merge categorical traits



#Reduce cover data to focal data using a series of merges.

#minimum number of replicates
repnum <- cover%>%
  dplyr::select(site_code, project_name, community_type, treatment, plot_id)%>%
  unique()%>%
  #dplyr::mutate(present = 1)%>%
  dplyr::group_by(site_code, project_name, community_type, treatment)%>%
  dplyr::summarise(rep_num = length(plot_id))%>%
  dplyr::ungroup()

#Last year of experiment
lastyear <- ddply(experimentinfo, .(site_code, project_name, community_type),
                  function(x)data.frame(
                    last_trt_yr = max(x$calendar_year)
                  ))

#minimum treatment length
nyear <- experimentinfo[c("site_code", "project_name", "community_type", "treatment_year")] %>%
  unique()%>%
  subset(treatment_year != 0)%>%
  ddply(.(site_code, project_name, community_type),
        function(x)data.frame(
          n.trt.yrs = max(x$treatment_year)
        ))

#Merge all the datasets above to create columns to subset by
crest <- cover %>%
  left_join(nyear,by = c("site_code", "project_name", "community_type"))%>%
  left_join(lastyear, by = c("site_code", "project_name", "community_type"))%>%
  left_join( experimentinfo, by = c("site_code", "project_name", "community_type", "treatment", "calendar_year", "treatment_year"))%>%
  left_join(repnum, by = c("site_code", "project_name", "community_type", "treatment"))

#subset by criteria
test <- crest %>%
  #subset(treats_wanted != "NA")%>%
  subset( rep_num >=5)%>%
  subset(last_trt_yr == calendar_year)%>%
  subset(treatment_year == n.trt.yrs)

test$trt_type <-  revalue(test$trt_type, c("N*P" = "mult_nutrient","CO2*temp" = "mult_GCD", "drought*CO2*temp" = "mult_GCD","irr*CO2" = "mult_GCD","irr*CO2*temp" = "mult_GCD","N*CO2*temp" = "mult_GCD","N*irr*CO2" = "mult_GCD", "mult_nutrient*irr" = "mult_GCD","N*irr*CO2*temp" = "mult_GCD", "N*CO2" = "mult_GCD","N*drought" = "mult_GCD","N*irr" = "mult_GCD","N*irr*temp" = "mult_GCD","N*temp" = "mult_GCD","mult_nutrient*temp" = "mult_GCD","N*P*temp" = "mult_GCD","drought*temp" = "mult_GCD","irr*temp" = "mult_GCD") ) #all expect for the first term are used for mult_GCD category which is no longer being used

test <- test%>%
  subset( trt_type == "control" | trt_type == "N" | trt_type == "P" | trt_type == "irr" | 
            trt_type == "drought"  | trt_type == "temp"| trt_type == "mult_nutrient" #|trt_type == "mult_GCD"| trt_type == "CO2"
  )#%>%  #keep only the focal treatments

#Set minimum treatment years. Note that criteria is relaxed for drought experiments for: reasons
N <-  subset(test[test$trt_type %in% "N",], n.trt.yrs >= 6)
P <-  subset(test[test$trt_type %in% "P",], n.trt.yrs >= 6)
irr <-  subset(test[test$trt_type %in% "irr",], n.trt.yrs >= 6)
#CO2 <-  subset(test[test$trt_type %in% "CO2",], n.trt.yrs >= 6)
temp <-  subset(test[test$trt_type %in% "temp",], n.trt.yrs >= 6)
mult_nutrient <-  subset(test[test$trt_type %in% "mult_nutrient",], n.trt.yrs >= 6)
#mult_GCD <-  subset(test[test$trt_type %in% "mult_GCD",], n.trt.yrs >= 6)
drought <-  subset(test[test$trt_type %in% "drought",], n.trt.yrs >= 4)
control <-  subset(test[test$trt_type %in% "control",], n.trt.yrs >= 4)

test <- bind_rows(N, P, irr, temp, mult_nutrient, drought, control
                  #, mult_GCD, CO2
                  )

test <- test[c("site_code", "project_name", "community_type", "treatment_year", "plot_id", "species_matched", "relcov", "trt_type", "plot_mani", "treatment")]%>%
  unique()

plot.treatment <- test[c("site_code", "project_name", "community_type", "plot_id", "trt_type", "treatment")]%>%
  unique()
plot.treatment <- tidyr::unite(plot.treatment, "rep", c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = FALSE)

df <- left_join(test, traits, by = "species_matched", all.x = TRUE)

df <- unite(df, rep, c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = FALSE)

df <- unite(df, expgroup, c("site_code", "project_name", "community_type"), sep = "::")

#a few lines to remove NAs from the continuous trait data
df$ok <- complete.cases(df[,c(#"stem_spec_density", #"leaf_N",#"leaf_P",
  "LDMC",
  "plant_height_vegetative",
  "SLA",
  "rooting_depth",
  "seed_dry_mass",
  "leaf_C.N"
)])
df <- subset(df, ok == TRUE)
df <- subset(df, species_matched != "NA")


########################
##Summarize sites being used

sites <- test%>%
  dplyr::select(site_code, project_name, community_type, treatment_year, trt_type, treatment)%>%
  unique()%>%
  subset(trt_type != "control")

sites <- unite(sites, temp, c("project_name", "community_type"), sep = "::", remove = FALSE)

##############################
####CREATING AND TESTING BETA DIVERSITY RESULTS

#'test' dataframe has all the cover data but only with focal sites and treatments and such

#For each treatment at each site, pull the treatment and control data, spread, calculate distance matrix, then betadisper 

kevin <- unite(test, expgroup, c("site_code", "project_name", "community_type"), remove = FALSE, sep = "::" ) #named Kevin because Kevin Wilcox helped make this loop
kevin <- subset(kevin, species_matched != "NA")%>%
  ddply(.(expgroup, site_code, project_name, community_type, treatment_year, plot_id, species_matched, trt_type, plot_mani, treatment),
        function(x)data.frame(
          relcov = sum(x$relcov)
        )) #with species names matched to trait data, separate observations in the cover data can become multiple observations of the same species, therefore, must sum cover values

expgroup_vector <- unique(kevin$expgroup)

distances_master <- {}

for(i in 1:length(expgroup_vector)) {
  temp.df <- subset(kevin, expgroup == expgroup_vector[i])
  
  temp.wide <- temp.df%>%
    pivot_wider(names_from = species_matched, values_from = relcov, values_fill = 0)
  temp.distances <- vegdist(temp.wide[10:ncol(temp.wide)], method = "bray")
  temp.mod <- betadisper(temp.distances, group = temp.wide$trt_type, type = "centroid")
  distances_temp <- data.frame(expgroup = expgroup_vector[i], trt_type = temp.wide$trt_type, treatment = temp.wide$treatment, plot_mani = temp.wide$plot_mani, dist = temp.mod$dist)
  #distances_temp <- subset(distances_temp, dist > 0.00000000001) #not necessary when cO2 treatment excluded
  #distances_temp$dist <- ifelse(distances_temp$dist > 0.00000000001, distances_temp$dist, 0.001) #changes value for single serc experiment where distance equals essentially 0 which doesn't work with response ratios
  distances_master <- rbind(distances_master, distances_temp )
  rm(temp.df, temp.wide, temp.distances, temp.mod, distances_temp)
}



mean.dist.df <- ddply(distances_master,.(expgroup, trt_type, treatment, plot_mani), function(x)data.frame( mean_dist = mean(x$dist)))

trt.df <- subset(mean.dist.df, plot_mani >= 1)%>%
  dplyr::rename(dist.trt = mean_dist)
con.df <- subset(mean.dist.df, plot_mani == 0)%>%
  dplyr::rename(dist.con = mean_dist)%>%
  dplyr::select(expgroup, dist.con)

lrr.df <- merge(trt.df, con.df, by = "expgroup", all.x = TRUE)%>%
  mutate(lrr = log(dist.trt/dist.con))%>%
  mutate(con_minus_trt = dist.trt/dist.con)

lrr.df.conf <- lrr.df%>%
  ddply(.(trt_type), function(x)data.frame(
    lrr.mean = mean(x$lrr),
    lrr.error = qt(0.975, df=length(x$trt_type)-1)*sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type)-1),
    num_experiments = length(x$expgroup)
  ))

lrr.df.conf$trt_type <- factor(lrr.df.conf$trt_type, levels = c("drought", "irr", "temp", "N", "P", "mult_nutrient"#, "mult_GCD", "CO2"
                                                                ))
#visualize
ggplot(lrr.df.conf, aes(trt_type, lrr.mean, color = trt_type))+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
    geom_pointrange(aes(ymin = lrr.mean-lrr.error, ymax = lrr.mean+lrr.error), size = 1.5)+
  xlab("")+
  ylab("Species composition LRR distance between plots within treatment")+
  scale_color_manual(values = c("#df0000","#0099f6", "orange", "#00b844","#f2c300","#6305dc", "black"))+
  theme_base()

ggplot(lrr.df, aes(trt_type, lrr, color = trt_type))+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  
      geom_beeswarm(cex = 2)+
  xlab("")+
  ylab("Species composition LRR distance between plots within treatment")+
  scale_color_manual(values = c("#df0000","#0099f6", "orange", "#00b844","#f2c300","#6305dc"))+
  theme_base()


#models to test results
distances_master.1 <- tidyr::separate(distances_master, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(distances_master.1, trt_type == "control" | trt_type == "drought"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(distances_master.1, trt_type == "control" | trt_type == "irr"))
summary(mod)
#mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(distances_master.1, trt_type == "control" | trt_type == "CO2"))
#summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(distances_master.1, trt_type == "control" | trt_type == "N"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(distances_master.1, trt_type == "control" | trt_type == "P"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(distances_master.1, trt_type == "control" | trt_type == "mult_nutrient"))
summary(mod)
#mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(distances_master.1, trt_type == "control" | trt_type == "mult_GCD"))
#summary(mod)
            
            
lrr.df_species <- lrr.df


###Summarize sites being used
sites <- test%>%
  dplyr::select(site_code, project_name, community_type, treatment_year, trt_type, treatment)%>%
  unique()%>%
  subset(trt_type != "control")




######
###Try the same stuff with traits but they include categorical traits
CoRRE_CWMtraits <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/paper 2_PD and FD responses/data/CoRRE_CWMtraits_12142022.csv") #for now I'll just use this for categorical traits
CoRRE_CWMtraits_cat <- CoRRE_CWMtraits[, c(   "site_code", "project_name","community_type", "plot_id", "treatment_year", "CWM.growth_form", "CWM.photosynthetic_pathway", "CWM.lifespan", "CWM.clonal", "CWM.mycorrhizal_type", "CWM.n_fixation")]

CoRRE_CWMtraits_cat <- tidyr::unite(CoRRE_CWMtraits_cat, "rep", c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = TRUE)

CoRRE_CWMtraits_cat$is.graminoid <- ifelse(CoRRE_CWMtraits_cat$CWM.growth_form == "graminoid", 1, 0)
CoRRE_CWMtraits_cat$is.C4 <- ifelse(CoRRE_CWMtraits_cat$CWM.photosynthetic_pathway == "C$", 1, 0)
CoRRE_CWMtraits_cat$is.perennial <- ifelse(CoRRE_CWMtraits_cat$CWM.lifespan == "perennial", 1, 0)
CoRRE_CWMtraits_cat$is.clonal <- ifelse(CoRRE_CWMtraits_cat$CWM.clonal == "yes", 1, 0)
CoRRE_CWMtraits_cat$is.AM <- ifelse(CoRRE_CWMtraits_cat$CWM.mycorrhizal_type == "arbuscular", 1, 0)
CoRRE_CWMtraits_cat$is.n_fixer <- ifelse(CoRRE_CWMtraits_cat$CWM.n_fixation == "yes", 1, 0)
CoRRE_CWMtraits_cat <- CoRRE_CWMtraits_cat[,c("rep", "treatment_year", "is.graminoid", "is.C4", "is.perennial", "is.clonal", "is.AM", "is.n_fixer")]


summarize.cwm <-   # New dataframe where we can inspect the result
  df %>%   # First step in the next string of statements
  dplyr::group_by(rep, expgroup, treatment_year, plot_id, trt_type, treatment, plot_mani) %>%   # Groups the summary file by Plot number
  dplyr::summarize(           # Coding for how we want our CWMs summarized
    seed_dry_mass.cwm = weighted.mean(seed_dry_mass, relcov),   # Actual calculation of CWMs
    LDMC.cwm = weighted.mean(LDMC, relcov),
    plant_height_vegetative.cwm = weighted.mean(plant_height_vegetative, relcov),
    SLA.cwm = weighted.mean(SLA, relcov),
    rooting_depth.cwm = weighted.mean(rooting_depth, relcov),
    leaf_C.N.cwm = weighted.mean(rooting_depth, relcov)
    )%>%
  left_join(CoRRE_CWMtraits_cat, by = c("rep", "treatment_year"))

summarize.traits.continuous <- traits[,c("species_matched", "seed_dry_mass", "LDMC", "plant_height_vegetative", "SLA", "rooting_depth", "leaf_C.N")]
summarize.traits.continuous <- unique(summarize.traits.continuous)
summarize.traits.categorical <- traits[,c("species_matched", "growth_form", "photosynthetic_pathway", "lifespan", "clonal", "mycorrhizal_type", "n_fixation")]
summarize.traits.categorical <- subset(summarize.traits.categorical, photosynthetic_pathway == "C3" | photosynthetic_pathway == "C4" | photosynthetic_pathway == "CAM")
                                       
summarize.traits <- left_join(summarize.traits.continuous, summarize.traits.categorical, by = "species_matched")                                       
# reassigning row names
summarize.traits <- unique(summarize.traits)



#########CALCULATE ALPHA FDIS TO COMPARE



##### calculate functional dispersion - loop through sites #####
#distance_df_master <- {}
#site_proj_comm_vector <- unique(df$expgroup)


#for(i in 1:length(site_proj_comm_vector)){

#  temp.df.comm <- subset(kevin, expgroup == site_proj_comm_vector[i])

  
  #species vector for pulling traits from relative cover
#  sp_df_temp <- data.frame(species_matched = unique(temp.df.comm$species_matched), dummy=1) 
#  site.traits.temp <- left_join(sp_df_temp, summarize.traits, by = "species_matched")
#  rownames(site.traits.temp) <- site.traits.temp$species_matched
#  site.traits.temp <- site.traits.temp[ , -which(names(site.traits.temp) %in% c("species_matched", "dummy"))]
  
#  temp.wide.comm <- temp.df.comm%>%
#    pivot_wider(names_from = species_matched, values_from = relcov, values_fill = 0)
  
#  temp.alpha <- dbFD(x = site.traits.temp, a = temp.wide.comm[10:ncol(temp.wide.comm)], m = 2)
  
  
  
#  sp_vec_temp <- sp_df_temp %>%
#    pull(species_matched)
  
  #subset trait data to just include the species present subset relative cover data
#  traits_df_raw_temp <- traitsClean %>%
#    filter(species_matched %in% sp_vec_temp)
  
  #dataframe with species in trait database and in relative cover data base
#  species_in_trait_data_temp <- data.frame(species_matched = unique(traits_df_raw_temp$species_matched),
#                                           dummy_traits=2) %>% #there are fewer species in the unique trait dataset than in the species comp data because there are thing like "unknown forb"
#    arrange(species_matched)
  
  #vector of species not in trait database (but in relative abundance data) to remove from species abundance data
#  sp_to_remove_temp <- sp_df_temp %>%
#    full_join(species_in_trait_data_temp, by="species_matched") %>%
#    filter(is.na(dummy_traits)) %>%
#    pull(genus_species)
  
  #abundance dataset with species removed that do not have trait information
#  relcov_unkn_sp_rm_temp <- relcov_df_temp %>%
#    filter(!genus_species %in% sp_to_remove_temp) #removing species without trait information
  
  #abundance data into wide format
#  relcov_wide_temp <- 
    
  # add rownames 
#  row.names(relcov_only_temp) <- paste(plot_info_temp$calendar_year, plot_info_temp$plot_id, sep="::")
  
  #dbFD function requires species names in trait data frame be arranged A-Z and identical order to the abundance data 
#  traits_df_temp <- traits_df_raw_temp %>%
#    arrange(species_matched) %>%
#    column_to_rownames("species_matched") %>%
#    dplyr::select(-family) %>%
#    mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
#    select(growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation, leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass)
  
  # change to dataframe from tibble 
#  traits_df_temp <- as.data.frame(traits_df_temp)
  
  #changing all categorical traits to factors
#  traits_df_temp[,c(1:6)] <- lapply(traits_df_temp[,c(1:6)], as.factor)
  
  #changing all continuous to numerical
#  traits_df_temp[,c(7:12)] <- lapply(traits_df_temp[,c(7:12)], as.numeric)
  
  ### Calculate MNTD and functional diversity metrics -- had to use Cailliez correlations becuase Euclidean distances could be calculated
#  FD_temp <- dbFD(x=traits_df_temp, a=relcov_only_temp, cor="cailliez", calc.FRic=F) # FRich is causing problems with most datasets (I think because of missing data?) so I'm removing it for now
  
#  FD_df_temp <- do.call(cbind.data.frame, FD_temp) %>%
#    mutate(year_plotid = row.names(.)) %>%
#    separate(year_plotid, into=c("calendar_year","plot_id"), sep="::") %>%
#    mutate(calendar_year = as.numeric(calendar_year)) %>%
#    full_join(plot_info_temp, by=c("calendar_year","plot_id"))
  
#  comp_matrix_temp <- as.matrix(relcov_only_temp)
#  trait_dist_temp <- as.matrix(gowdis(traits_df_temp))
  
#  MNTD_df_temp <- data.frame(
#    plot_info_temp[,c("calendar_year", "plot_id")],
#    MNTD_traits = picante::mntd(comp_matrix_temp, trait_dist_temp)
#  )
  
#  distance_df_temp <- FD_df_temp %>%
#    full_join(MNTD_df_temp, by=c("calendar_year","plot_id"))
  
#  distance_df_master <- rbind(distance_df_master, distance_df_temp)
  
#  rm(list=ls()[grep("temp", ls())])
#}





###########CALCULATE BETA DIVERSITY WITH TRAITS
expgroup_vector <- unique(df$expgroup)

tdistances_master <- {}

for(i in 1:length(expgroup_vector)) {
  temp.df <- subset(summarize.cwm, expgroup == expgroup_vector[i])
  temp.gow <- gowdis(temp.df[7:ncol(temp.df)])
  temp.beta <- betadisper(temp.gow, group = temp.df$trt_type, type = "centroid")
  tdistances_temp <- data.frame(expgroup = expgroup_vector[i], trt_type = temp.df$trt_type, treatment = temp.df$treatment,  dist = temp.beta$dist, plot_mani = temp.df$plot_mani)
#  tdistances_temp <- subset(tdistances_temp, dist > 0.00000000001) #not necesssary when excluding CO2 treatment
#  tdistances_temp$dist <- ifelse(tdistances_temp$dist > 0.00000000001, tdistances_temp$dist, 0.001) #changes value for single serc experiment where distance equals essentially 0 which doesn't work with response ratios
  tdistances_master <- rbind(tdistances_master, tdistances_temp )
  rm(temp.df, temp.gow, temp.beta, tdistances_temp)
  
}

mean.dist.df <- ddply(tdistances_master,.(expgroup, trt_type, treatment, plot_mani), function(x)data.frame( mean_dist = mean(x$dist)))

trt.df <- subset(mean.dist.df, plot_mani >= 1)%>%
  dplyr::rename(dist.trt = mean_dist)
con.df <- subset(mean.dist.df, plot_mani == 0)%>%
  dplyr::rename(dist.con = mean_dist)%>%
  dplyr::select(expgroup, dist.con)

lrr.df <- merge(trt.df, con.df, by = "expgroup", all.x = TRUE)%>%
  mutate(lrr = log(dist.trt/dist.con))%>%
  mutate(con_minus_trt = dist.trt/dist.con)

lrr.df.conf <- lrr.df%>%
  ddply(.(trt_type), function(x)data.frame(
    lrr.mean = mean(x$lrr),
    lrr.error = qt(0.975, df=length(x$trt_type)-1)*sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type)-1),
    num_experiments = length(x$expgroup)
  ))

#visualize
lrr.df.conf$trt_type <- factor(lrr.df.conf$trt_type, levels = c("drought", "irr", "temp", "N", "P", "mult_nutrient" 
                                                                #,"mult_GCD", "CO2"
                                                                ))
ggplot(lrr.df.conf, aes(trt_type, lrr.mean, color = trt_type))+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  geom_pointrange(aes(ymin = lrr.mean-lrr.error, ymax = lrr.mean+lrr.error), size = 1.5)+
  xlab("")+
  ylab("Trait LRR distance between plots within treatment")+
  scale_color_manual(values = c("#df0000","#0099f6", "orange", "#00b844","#f2c300","#6305dc"))+
  theme_base()

ggplot(lrr.df, aes(trt_type, lrr, color = trt_type))+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  geom_beeswarm(cex = 2)+
  xlab("")+
  ylab("Trait composition LRR distance between plots within treatment")+
  scale_color_manual(values = c("#df0000","#0099f6", "orange","#00b844","#f2c300","#6305dc"))+
  theme_base()

lrr.df_traits <- lrr.df

#models to test results
tdistances_master.1 <- tidyr::separate(tdistances_master, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)
tdistances_master.1 <- tidyr::separate(tdistances_master, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(tdistances_master.1, trt_type == "control" | trt_type == "drought"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(tdistances_master.1, trt_type == "control" | trt_type == "irr"))
summary(mod)
#mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(tdistances_master.1, trt_type == "control" | trt_type == "CO2"))
#summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(tdistances_master.1, trt_type == "control" | trt_type == "N"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(tdistances_master.1, trt_type == "control" | trt_type == "P"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(tdistances_master.1, trt_type == "control" | trt_type == "mult_nutrient"))
summary(mod)
#mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(tdistances_master.1, trt_type == "control" | trt_type == "mult_GCD"))
#summary(mod)

#####################
###Compare species and trait responses
lrr_sp.tr <- merge(lrr.df_species, lrr.df_traits, by = c("expgroup", "trt_type", "treatment", "plot_mani"), all = TRUE)%>%
  dplyr::rename(c(lrr.species = lrr.x, lrr.traits = lrr.y))

#visualize
lrr_sp.tr$trt_type <- factor(lrr_sp.tr$trt_type, levels = c("drought", "irr", "temp", "N", "P", "mult_nutrient" 
                                                                #,"mult_GCD", "CO2"
))
ggplot(lrr_sp.tr, aes(lrr.species, lrr.traits))+
  facet_wrap(~trt_type)+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
  geom_vline(xintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
  ylim(-1.1, 2.4)+
  xlim(-1.1, 2.4)+
  theme_base()

#models to test results
summary(lm(lrr.traits~lrr.species, data = subset(lrr_sp.tr, trt_type == "drought")))
summary(lm(lrr.traits~lrr.species, data = subset(lrr_sp.tr, trt_type == "irr")))
#summary(lm(lrr.traits~lrr.species, data = subset(lrr_sp.tr, trt_type == "CO2")))
summary(lm(lrr.traits~lrr.species, data = subset(lrr_sp.tr, trt_type == "N")))
summary(lm(lrr.traits~lrr.species, data = subset(lrr_sp.tr, trt_type == "P")))
summary(lm(lrr.traits~lrr.species, data = subset(lrr_sp.tr, trt_type == "mult_nutrient")))
#summary(lm(lrr.traits~lrr.species, data = subset(lrr_sp.tr, trt_type == "mult_GCD")))
summary(lm(lrr.traits~lrr.species, data = subset(lrr_sp.tr, trt_type == "temp")))


############
##BRING IN COVARIATES AND SEE IF THEY EXPLAIN BETA DIVERSITY
CoRRE_siteLocationClimate_Dec2021 <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/environmental data/CoRRE_siteLocationClimate_Dec2021.csv")

CoRRE_project_summary <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/CoRRE_project_summary.csv")
CoRRE_project_summary$project <- CoRRE_project_summary$project_name
CoRRE_project_summary$community <- CoRRE_project_summary$community_type
CoRRE_project_summary <- CoRRE_project_summary %>% dplyr::select(-c(project_name, community_type))


lrr.df_species <- tidyr::separate(lrr.df_species, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)
lrr.df_traits <- tidyr::separate(lrr.df_traits, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)

lrr_covariate <- left_join(lrr.df_species, CoRRE_project_summary, by = c("site_code", "project", "community"))
lrr_covariate_traits <- left_join(lrr.df_traits, CoRRE_project_summary, by = c("site_code", "project", "community"))

##DROUGHT
drought <- subset(lrr_covariate, trt_type == "drought")
mod <- lmer(lrr~MAP + (1|site_code), data = drought)
#mod <- lm(lrr~MAP, data = drought)
summary(mod)
visreg(mod, ylab = "lrr beta diversity", main = "Drought treatment")
mod <- lmer(lrr~MAT + (1|site_code), data = drought)
summary(mod)
visreg(mod, ylab = "lrr beta diversity", main = "Drought treatment")
mod <- lmer(lrr~rrich + (1|site_code), data = drought)
summary(mod)
visreg(mod,  ylab = "lrr beta diversity", main = "Drought treatment")
mod <- lmer(lrr~experiment_length + (1|site_code), data = drought)
summary(mod)
visreg(mod,  ylab = "lrr beta diversity", main = "Drought treatment")

##IRRIGATION
irrigation <- subset(lrr_covariate, trt_type == "irr")
mod <- lmer(lrr~MAP + (1|site_code), data = irrigation)
summary(mod)
visreg(mod, ylab = "lrr beta diversity", main = "Irrigation treatment")
mod <- lmer(lrr~MAT + (1|site_code), data = irrigation)
summary(mod)
visreg(mod)
mod <- lmer(lrr~rrich + (1|site_code), data = irrigation)
summary(mod)
visreg(mod)
mod <- lmer(lrr~experiment_length + (1|site_code), data = irrigation)
summary(mod)
visreg(mod)

##NITROGEN
N <- subset(lrr_covariate, trt_type == "N")
mod <- lmer(lrr~MAP + (1|site_code), data = N)
summary(mod)
visreg(mod)
mod <- lmer(lrr~MAT + (1|site_code), data = N)
summary(mod)
visreg(mod)
mod <- lmer(lrr~rrich + (1|site_code), data = N)
summary(mod)
visreg(mod)
mod <- lmer(lrr~experiment_length + (1|site_code), data = N)
summary(mod)
visreg(mod)

##PHOSPORUS
P <- subset(lrr_covariate, trt_type == "P")
mod <- lmer(lrr~MAP + (1|site_code), data = P)
summary(mod)
visreg(mod)
mod <- lmer(lrr~MAT + (1|site_code), data = P)
summary(mod)
visreg(mod)
mod <- lmer(lrr~rrich + (1|site_code), data = P)
summary(mod)
visreg(mod)
mod <- lmer(lrr~experiment_length + (1|site_code), data = P)
summary(mod)
visreg(mod)


#MULT_NUTRIENT
mult_nutrient <- subset(lrr_covariate, trt_type == "mult_nutrient")
mod <- lmer(lrr~MAP + (1|site_code), data = mult_nutrient)
summary(mod)
mod <- lmer(lrr~MAT + (1|site_code), data = mult_nutrient)
summary(mod)
visreg(mod)
mod <- lmer(lrr~rrich + (1|site_code), data = mult_nutrient)
summary(mod)
visreg(mod)
visreg(mod, xvar = "rrich", yvar = "lrr", ylab = "lrr beta diversity", main = "Mutliple nutrient addition", gg = TRUE)+
  theme_base()
mod <- lmer(lrr~experiment_length + (1|site_code), data = mult_nutrient)
summary(mod)
visreg(mod)

###trait models
##DROUGHT
drought <- subset(lrr_covariate_traits, trt_type == "drought")
mod <- lmer(lrr~MAP + (1|site_code), data = drought)
summary(mod)
visreg(mod)
mod <- lmer(lrr~MAT + (1|site_code), data = drought)
summary(mod)
visreg(mod)
mod <- lmer(lrr~rrich + (1|site_code), data = drought)
summary(mod)
visreg(mod)
mod <- lmer(lrr~experiment_length + (1|site_code), data = drought)
summary(mod)
visreg(mod)

##NITROGEN
N <- subset(lrr_covariate_traits, trt_type == "N")
mod <- lmer(lrr~MAP + (1|site_code), data = N)
summary(mod)
visreg(mod)
mod <- lmer(lrr~MAT + (1|site_code), data = N)
summary(mod)
visreg(mod)
mod <- lmer(lrr~rrich + (1|site_code), data = N)
summary(mod)
visreg(mod)
mod <- lmer(lrr~experiment_length + (1|site_code), data = N)
summary(mod)
visreg(mod)

##mult_nutrient
mult_nutrient <- subset(lrr_covariate_traits, trt_type == "mult_nutrient")
mod <- lmer(lrr~MAP + (1|site_code), data = mult_nutrient)
summary(mod)
visreg(mod)
mod <- lmer(lrr~MAT + (1|site_code), data = mult_nutrient)
summary(mod)
visreg(mod)
mod <- lmer(lrr~rrich + (1|site_code), data = mult_nutrient)
summary(mod)
visreg(mod)
mod <- lmer(lrr~experiment_length + (1|site_code), data = mult_nutrient)
summary(mod)
visreg(mod)


##with treatment information
treatment_info <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/basic dataset info/ExperimentInfo.csv")
treatment_info$trt_type <- revalue(treatment_info$trt_type, c("N*P" = "mult_nutrient"))
treatment_info$project <- treatment_info$project_name
treatment_info$community <- treatment_info$community_type
treatment_info <- treatment_info[,c("site_code", "project", "community", "plot_mani", "trt_type", "treatment", "n", "p",  "CO2", "precip", "temp")]
treatment_info <- unique(treatment_info)
lrr_treat_species <- left_join(lrr_covariate, treatment_info, by = c("site_code", "project", "community", "plot_mani", "trt_type", "treatment"))
lrr_treat_traits <- left_join(lrr_covariate_traits, treatment_info, by = c("site_code", "project", "community", "plot_mani", "trt_type", "treatment"))


##ANY WATER MANIPULATION
water_mani <- subset(lrr_treat_species, trt_type == "drought"| trt_type == "irr")
mod <- lmer(lrr~precip + (1|expgroup) ,data = water_mani)
summary(mod)
visreg(mod)
ggplot(water_mani, aes(precip, lrr, color = MAP))+
  geom_point(size = 2)+
  geom_smooth(method = "lm", se = FALSE)+
  theme_base()

##Nitrogen gradient
temp <- subset(lrr_treat_species, trt_type == "N")
mod <- lmer(lrr~n + (1|expgroup) ,data = temp)
mod <- lmer(lrr~n*MAP + (1|expgroup), data = temp)
summary(mod)
visreg(mod)

##traits
water_mani <- subset(lrr_treat_traits, trt_type == "drought" | trt_type == "irr")
mod <- lmer(lrr~precip + (1|expgroup) ,data = water_mani)
summary(mod)
ggplot(water_mani, aes(precip, lrr, color = MAP))+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
  geom_vline(xintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
  geom_point(size = 2)+
  ylab("LRR trait beta diversity")+
  xlab("Precip treatment as percentage of MAP")+
  geom_smooth(method = "lm", se = FALSE, color = "black")+
  theme_base()

visreg(mod, xvar = "precip", yvar = "lrr", ylab = "lrr trait beta diversity", xlab = "Precipitation treatment", gg = TRUE)+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  theme_base()
  
N <- subset(lrr_treat_traits, trt_type == "N")
mod <- lmer(lrr~n  + (1|expgroup) ,data = N)
summary(mod)
ggplot(N, aes(n, lrr, color = dist.con))+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
  geom_point(size = 2)+
  ylab("LRR trait beta diversity")+
  xlab("N addition treatment (units)")+
  #geom_smooth(method = "lm", se = FALSE)+ #makes sense to remove the geom_smooth layer as long as it's not a significant relationship
  theme_base()
visreg(mod, xvar = "n", yvar = "lrr", ylab = "lrr beta diversity", xlab = "Nitrogen application", gg = TRUE)+
  geom_hline(yintercept = 0)+
  theme_base()


###########################################
####WHAT CONTRIBUTES TO BETA DIVERSITY CHANGE??


####RANK DIFF AMONG REPLICATES (COMMUNITY)


kevin <- unite(test, expgroup, c("site_code", "project_name", "community_type"), remove = FALSE, sep = "::" ) #named Kevin because Kevin Wilcox helped make this loop
kevin <- subset(kevin, species_matched != "NA")%>%
  ddply(.(expgroup, site_code, project_name, community_type, treatment_year, plot_id, species_matched, trt_type, plot_mani, treatment),
        function(x)data.frame(
          relcov = sum(x$relcov)
        )) #with species names matched to trait data, separate observations in the cover data can become multiple observations of the same species, therefore, must sum cover values

expgroup_vector <- unique(kevin$expgroup)

rank_diff_master <- {}

for(i in 1:length(expgroup_vector)) {
  temp.df <- subset(kevin, expgroup == expgroup_vector[i])
  temp.expinfo <- temp.df[,c("expgroup", "site_code", "project_name", "community_type", "plot_id", "trt_type", "plot_mani", "treatment")]%>%
                  unique()
  
rank_diff_temp <-  RAC_difference(
    df = temp.df,
    time.var = "treatment_year",
    species.var = "species_matched",
    abundance.var = "relcov",
    replicate.var = "plot_id",
    treatment.var = "treatment",
    pool = FALSE,
    block.var = NULL,
    reference.treatment = NULL
  )%>%
   subset( treatment == treatment2)%>%
   left_join(temp.expinfo, by = c("plot_id", "treatment"))
 

  rank_diff_master <- rbind(rank_diff_master, rank_diff_temp )
  rm(temp.df, temp.expinfo, rank_diff_temp)
}


mean.diff.df <- ddply(rank_diff_master,.(expgroup, trt_type, treatment, plot_mani), function(x)data.frame( mean_diff = mean(x$rank_diff)))

trt.df <- subset(mean.diff.df, plot_mani >= 1)%>%
  dplyr::rename(diff.trt = mean_diff)
con.df <- subset(mean.diff.df, plot_mani == 0)%>%
  dplyr::rename(diff.con = mean_diff)%>%
  dplyr::select(expgroup, diff.con)

lrr.df <- merge(trt.df, con.df, by = "expgroup", all.x = TRUE)%>%
  mutate(lrr = log(diff.trt/diff.con))%>%
  mutate(con_minus_trt = diff.trt/diff.con)

lrr.df.conf <- lrr.df%>%
  ddply(.(trt_type), function(x)data.frame(
    lrr.mean = mean(x$lrr, na.rm = TRUE),
    lrr.error = qt(0.975, df=length(x$trt_type)-1)*sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type)-1),
    num_experiments = length(x$expgroup)
  ))

lrr.df.conf$trt_type <- factor(lrr.df.conf$trt_type, levels = c("drought", "irr", "temp", "N", "P", "mult_nutrient"#, "mult_GCD", "CO2"
))
          
ggplot(lrr.df.conf, aes(trt_type, lrr.mean, color = trt_type))+ #this figure doesn't match up with the model results below
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  geom_pointrange(aes(ymin = lrr.mean-lrr.error, ymax = lrr.mean+lrr.error), size = 1.5)+
  xlab("")+
  ylab("LRR rank difference between replicated")+
  scale_color_manual(values = c("#df0000","#0099f6", "orange", "#00b844","#f2c300","#6305dc", "black"))+
  theme_base()                                            
                                                                
#models to test results
rank_diff_master.1 <- tidyr::separate(rank_diff_master, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)
mod <- lmer(rank_diff~trt_type + (1|site_code/expgroup), data = subset(rank_diff_master.1, trt_type == "control" | trt_type == "drought"))
summary(mod)
mod <- lmer(rank_diff~trt_type + (1|site_code/expgroup), data = subset(rank_diff_master.1, trt_type == "control" | trt_type == "irr"))
summary(mod)
mod <- lmer(rank_diff~trt_type + (1|site_code/expgroup), data = subset(rank_diff_master.1, trt_type == "control" | trt_type == "N"))
summary(mod)
mod <- lmer(rank_diff~trt_type + (1|site_code/expgroup), data = subset(rank_diff_master.1, trt_type == "control" | trt_type == "P"))
summary(mod)
mod <- lmer(rank_diff~trt_type + (1|site_code/expgroup), data = subset(rank_diff_master.1, trt_type == "control" | trt_type == "mult_nutrient"))
summary(mod)
                                                                

####SINGLE TRAIT VARIANCE AMONG REPLICATES (TRAIT)





