##A little digression to make sure that mult nutrient effects aren't just nitrogn effects
###the next file after tim_betadiv_v2. used to calculate local distances of species composition and traits. Includes statistical analyses and figures 

library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggeffects)
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
library(emmeans)
library(fixest)
library(nlme)
library(RcmdrMisc)
library(lmtest)

#Read in trait data
traits_cat <- read.csv('https://pasta.lternet.edu/package/data/eml/edi/1533/3/5ebbc389897a6a65dd0865094a8d0ffd')%>%
  dplyr::select(-family, -source, -error_risk_overall)%>%
  pivot_wider(names_from = trait, values_from = trait_value)%>%
  dplyr::rename(species_matched = species)

traits <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/1533/3/169fc12d10ac20b0e504f8d5ca0b8ee8")%>% #continuous trait data
  mutate(species_matched = species)%>%
  dplyr::select(species_matched, trait, trait_value)%>%
  pivot_wider(names_from = trait, values_from = trait_value)

#standardize the scale of continuous traits
cols <- c( 
  "SLA", 
  "LDMC", 
  "leaf_N", 
  "plant_height_vegetative", 
  "seed_dry_mass", 
  "SRL"
)

traits[cols] <- scale(traits[cols])

traits <- left_join(traits, traits_cat, by = "species_matched")#merge w/ categorical traits

##Read in cover data
cover <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_RelativeCoverMarch2024.csv") %>% #community comp relative cover data
  mutate(drop=ifelse(site_code=="CDR"&treatment==2|site_code=="CDR"&treatment==3|site_code=="CDR"&treatment==4|site_code=="CDR"&treatment==5|site_code=="CDR"&treatment==7, 1,0))%>%
  filter(drop==0)%>% #remove some Cedar Creek treatments since that site is somewhat overrepresented
  subset(treatment_year <=10& treatment_year > 0) #only use treatment data and subset the number of years to be used

corre2trykey <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/corre2trykey_2021.csv") #matched species names among trait data and relative cover data
corre2trykey <- corre2trykey[,c("genus_species","species_matched")]
corre2trykey <- unique(corre2trykey)
cover <- left_join(cover, corre2trykey, by = "genus_species", keep = FALSE)

experimentinfo <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfo_March2024.csv")%>%#Information about the treatments which gets used to test how treatment magnitude explains efect sizes
  unique()

#siteLocationClimate <- read.csv("C:/Users/ohler/Dropbox/CoRRE_database/Data/CompiledData/siteLocationClimate.csv") #information about sites


#Reduce cover data to focal data using a series of merges.
#minimum number of replicates
repnum <- cover%>%
  dplyr::select(site_code, project_name, community_type, treatment, plot_id)%>%
  unique()%>%
  #dplyr::mutate(present = 1)%>%
  dplyr::group_by(site_code, project_name, community_type, treatment)%>%
  dplyr::summarise(rep_num = length(plot_id))%>%
  dplyr::ungroup()
#Merge all the datasets above to create columns to subset by
crest <- cover %>%
  left_join( experimentinfo, by = c("site_code", "project_name", "community_type", "treatment", "calendar_year", "treatment_year"))%>%
  left_join(repnum, by = c("site_code", "project_name", "community_type", "treatment"))

#subset by criteria
test <- crest %>%
  subset( rep_num >=5)

test$trt_type <-  revalue(test$trt_type, c("N*P" = "mult_nutrient",
                                           "CO2*temp" = "mult_GCD", 
                                           "drought*CO2*temp" = "mult_GCD",
                                           #"irr*CO2" = "mult_GCD",
                                           "irr*CO2*temp" = "mult_GCD",
                                           "N*CO2*temp" = "mult_GCD",
                                           #"N*irr*CO2" = "mult_GCD", 
                                           #"mult_nutrient*irr" = "mult_GCD",
                                           "N*irr*CO2*temp" = "mult_GCD", 
                                           #"N*CO2" = "mult_GCD",
                                           "N*drought" = "mult_GCD",
                                           #"N*irr" = "mult_GCD",
                                           "N*irr*temp" = "mult_GCD",
                                           "N*temp" = "mult_GCD",
                                           "mult_nutrient*temp" = "mult_GCD",
                                           "N*P*temp" = "mult_GCD",
                                           "drought*temp" = "mult_GCD",
                                           "irr*temp" = "mult_GCD") ) #all expect for the first term are used for mult_GCD category which is no longer being used. basically used to subset out all the mult_gcd category below

test <- test%>%
  subset( trt_type == "control" | trt_type == "N" | trt_type == "P" | trt_type == "irr" |# trt_type == "drought"  | trt_type == "temp"| 
            trt_type == "mult_nutrient" #|trt_type == "mult_GCD"
          | trt_type == "CO2" | trt_type == "irr*CO2"  |trt_type == "N*irr*CO2" | trt_type == "mult_nutrient*irr" |trt_type == "N*CO2"|trt_type == "N*irr"
  )  #keep only the focal treatments

test <- test[c("site_code", "project_name", "community_type", "treatment_year", "plot_id", "species_matched", "relcov", "trt_type", "plot_mani", "treatment")]%>%
  unique()

plot.treatment <- test[c("site_code", "project_name", "community_type", "plot_id", "trt_type", "treatment")]%>%
  unique()
plot.treatment <- tidyr::unite(plot.treatment, "rep", c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = FALSE)

df <- left_join(test, traits, by = "species_matched", keep = FALSE)

df <- unite(df, rep, c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = FALSE)

df <- unite(df, expgroup, c("site_code", "project_name", "community_type"), sep = "::")

#a few lines to remove NAs from the continuous trait data
df$ok <- complete.cases(df[,c("SLA", "LDMC", "leaf_N", "plant_height_vegetative", "seed_dry_mass", "SRL"
)])
df <- subset(df, ok == TRUE)

########################
##Summarize sites being used
sites <- test%>%
  dplyr::select(site_code, project_name, community_type, treatment_year, trt_type, treatment)%>%
  unique()%>%
  subset(trt_type != "control")
length(unique(sites$site_code))

sites <- unite(sites, temp, c("site_code", "project_name", "community_type"), sep = "::", remove = FALSE)

length(unique(sites$temp))

sites <- unite(sites, temp, c("site_code", "project_name", "community_type", "trt_type", "treatment"), sep = "::", remove = FALSE)

#write.csv(dplyr::select(sites, site_code,project_name)%>%unique(), "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/sites_table_local.csv")

length(unique(sites$temp))

length(unique(test$species_matched))

###Summarize sites being used
sites <- test%>%
  dplyr::select(site_code, project_name, community_type, trt_type, treatment)%>%
  tidyr::unite("expgroup", c("site_code", "project_name", "community_type"))%>%
  unique()%>%
  subset(trt_type != "control")

n <- sites%>%
  ddply(.(trt_type), function(x)data.frame(n = length(x$expgroup)))




##############################
####CREATING AND TESTING BETA DIVERSITY RESULTS

#'test' dataframe has all the cover data but only with focal sites and treatments and such
#For each treatment at each site, pull the treatment and control data, spread, calculate distance matrix, then betadisper 

kevin <- unite(test, expgroup, c("site_code", "project_name", "community_type"), remove = FALSE, sep = "::" ) #named Kevin because Kevin Wilcox helped make this loop
kevin <- subset(kevin, species_matched != "NA")%>%
  ddply(.(expgroup, site_code, project_name, community_type, treatment_year, plot_id, species_matched, trt_type, plot_mani, treatment, treatment_year),
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
  temp.mod <- betadisper(temp.distances, group = temp.wide$treatment, type = "centroid")
  distances_temp <- data.frame(expgroup = expgroup_vector[i], trt_type = temp.wide$trt_type, treatment = temp.wide$treatment, plot_mani = temp.wide$plot_mani, dist = temp.mod$dist, treatment_year = temp.wide$treatment_year)
  distances_master <- rbind(distances_master, distances_temp )
  rm(temp.df, temp.wide, temp.distances, temp.mod, distances_temp)
}


mean.dist.df <- ddply(distances_master,.(expgroup, trt_type, treatment, plot_mani, treatment_year), function(x)data.frame( mean_dist = mean(x$dist)))%>%
  separate(expgroup, into = c("site", "project", "community"), sep = "::", remove = FALSE)

mean.dist.comp <- mean.dist.df


mult_nut <- subset(mean.dist.comp, trt_type == "mult_nutrient")
control <- subset(mean.dist.comp, trt_type == "control")
n <- subset(mean.dist.comp, trt_type == "N")


##treatment info
treatment_info <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/basic dataset info/ExperimentInfo.csv")
treatment_info$trt_type <- revalue(treatment_info$trt_type, c("N*P" = "mult_nutrient"))
treatment_info$project <- treatment_info$project_name
treatment_info$community <- treatment_info$community_type
treatment_info <- treatment_info[,c("site_code", "project", "community", "plot_mani", "trt_type", "treatment", "n", "p",  "CO2", "precip", "temp")]
treatment_info <- unique(treatment_info)%>%
  unite( expgroup, c("site_code", "project", "community"), sep = "::")%>%
  group_by(expgroup, plot_mani, trt_type, treatment)%>%
  dplyr::summarize(n = mean(n), p = mean(p), CO2 = mean(CO2), precip = mean(precip), temp = mean(temp))

n_comb <- left_join(n, treatment_info, by = c("expgroup", "treatment", "trt_type"))%>%
  mutate(set = "n_only")%>%
  dplyr::select(mean_dist, set, site, expgroup,treatment_year)

comb <- left_join(mult_nut, treatment_info, by = c("expgroup", "treatment", "trt_type"))


mult_full <- subset(comb, trt_type == "mult_nutrient")%>%
  mutate(set = "full")%>%
  dplyr::select(mean_dist, set, site, expgroup,treatment_year)

mult_n_only <- subset(comb, trt_type == "mult_nutrient"& n >0)%>%
  mutate(set = "n_and_friends")%>%
  dplyr::select(mean_dist, set, site, expgroup,treatment_year)

mult_minus_n <- subset(comb, trt_type == "mult_nutrient"& n ==0)%>%
  mutate(set = "no_n")%>%
  dplyr::select(mean_dist, set, site, expgroup,treatment_year)

control2 <- control%>%
  mutate(set = "control")%>%
  dplyr::select(mean_dist, set, site, expgroup,treatment_year)

comb2 <- rbind(mult_full, mult_n_only)%>%
          rbind(mult_minus_n)%>%
          rbind(n_comb)%>%
          rbind(control2)

#check with numbers
mod <- feols(mean_dist ~ set | site+expgroup+as.character(treatment_year), data = comb2)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))

x <- ggpredict(mod, "set")
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  ylim(0.3,0.5)+
  ylab("Taxonomic dispersion within sites")+
  theme_base()


###########
#
#
#
#
#
#

######
###The same stuff with traits but they include categorical traits

#turn categorical traits into binary
df$growth_form <- ifelse(df$growth_form == "graminoid", 1, 0)
df$photosynthetic_pathway <- ifelse(df$photosynthetic_pathway == "C4", 1, 0)
df$lifespan <- ifelse(df$lifespan == "perennial", 1, 0)
df$clonal <- ifelse(df$clonal == "yes", 1, 0)
df$mycorrhizal_type <- ifelse(df$mycorrhizal_type == "AM", 1, 0)
df$n_fixation_type <- ifelse(df$n_fixation_type == "rhizobial", 1, 0)


#the below chunk collates the raw trait data into CWMs and adds the categorical data
summarize.cwm <-   
  df %>%   # First step in the next string of statements
  dplyr::group_by(rep, expgroup, treatment_year, plot_id, trt_type, treatment, plot_mani) %>%   # Groups the summary file by Plot number
  dplyr::summarize(           # Coding for how we want our CWMs summarized
    SLA.cwm = weighted.mean(SLA, relcov),
    LDMC.cwm = weighted.mean(LDMC, relcov),
    leaf_N.cwm = weighted.mean(leaf_N, relcov),
    plant_height_vegetative.cwm = weighted.mean(plant_height_vegetative, relcov),
    seed_dry_mass.cwm = weighted.mean(seed_dry_mass, relcov),   # Actual calculation of CWMs
    SRL.cwm = weighted.mean(SRL, relcov),
    growth_form.cwm = weighted.mean(growth_form, relcov),
    photosynthetic_pathway.cwm = weighted.mean(photosynthetic_pathway, relcov),
    lifespan.cwm = weighted.mean(lifespan, relcov),
    clonal.cwm = weighted.mean(clonal, relcov),
    mycorrhizal_type.cwm = weighted.mean(mycorrhizal_type, relcov),
    n_fixation_type.cwm = weighted.mean(n_fixation_type, relcov)
    
    
  )


###########CALCULATE BETA DIVERSITY WITH TRAITS
expgroup_vector <- unique(summarize.cwm$expgroup)

tdistances_master <- {}

for(i in 1:length(expgroup_vector)) {
  temp.df <- subset(summarize.cwm, expgroup == expgroup_vector[i])
  temp.gow <- gowdis(temp.df[8:ncol(temp.df)
  ])
  temp.beta <- betadisper(temp.gow, group = temp.df$trt_type, type = "centroid")
  tdistances_temp <- data.frame(expgroup = expgroup_vector[i], trt_type = temp.df$trt_type, treatment = temp.df$treatment,  dist = temp.beta$dist, plot_mani = temp.df$plot_mani, treatment_year = temp.df$treatment_year)
  #  tdistances_temp <- subset(tdistances_temp, dist > 0.00000000001) #not necesssary when excluding CO2 treatment
  #  tdistances_temp$dist <- ifelse(tdistances_temp$dist > 0.00000000001, tdistances_temp$dist, 0.001) #changes value for single serc experiment where distance equals essentially 0 which doesn't work with response ratios
  tdistances_master <- rbind(tdistances_master, tdistances_temp )
  rm(temp.df, temp.gow, temp.beta, tdistances_temp)
  
}

mean.dist.df <- ddply(tdistances_master,.(expgroup, trt_type, treatment, plot_mani, treatment_year), function(x)data.frame( mean_dist = mean(x$dist)))%>%
  separate(expgroup, into = c("site", "project", "community"), sep = "::", remove = FALSE)

mean.dist.trait <- mean.dist.df


mult_nut <- subset(mean.dist.trait, trt_type == "mult_nutrient")
control <- subset(mean.dist.trait, trt_type == "control")
n <- subset(mean.dist.trait, trt_type == "N")


##treatment info
treatment_info <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/basic dataset info/ExperimentInfo.csv")
treatment_info$trt_type <- revalue(treatment_info$trt_type, c("N*P" = "mult_nutrient"))
treatment_info$project <- treatment_info$project_name
treatment_info$community <- treatment_info$community_type
treatment_info <- treatment_info[,c("site_code", "project", "community", "plot_mani", "trt_type", "treatment", "n", "p",  "CO2", "precip", "temp")]
treatment_info <- unique(treatment_info)%>%
  unite( expgroup, c("site_code", "project", "community"), sep = "::")%>%
  group_by(expgroup, plot_mani, trt_type, treatment)%>%
  dplyr::summarize(n = mean(n), p = mean(p), CO2 = mean(CO2), precip = mean(precip), temp = mean(temp))

n_comb <- left_join(n, treatment_info, by = c("expgroup", "treatment", "trt_type"))%>%
  mutate(set = "n_only")%>%
  dplyr::select(mean_dist, set, site, expgroup,treatment_year)

comb <- left_join(mult_nut, treatment_info, by = c("expgroup", "treatment", "trt_type"))


mult_full <- subset(comb, trt_type == "mult_nutrient")%>%
  mutate(set = "full")%>%
  dplyr::select(mean_dist, set, site, expgroup,treatment_year)

mult_n_only <- subset(comb, trt_type == "mult_nutrient"& n >0)%>%
  mutate(set = "n_and_friends")%>%
  dplyr::select(mean_dist, set, site, expgroup,treatment_year)

mult_minus_n <- subset(comb, trt_type == "mult_nutrient"& n ==0)%>%
  mutate(set = "no_n")%>%
  dplyr::select(mean_dist, set, site, expgroup,treatment_year)

control2 <- control%>%
  mutate(set = "control")%>%
  dplyr::select(mean_dist, set, site, expgroup,treatment_year)

comb2 <- rbind(mult_full, mult_n_only)%>%
  rbind(mult_minus_n)%>%
  rbind(n_comb)%>%
  rbind(control2)

#check with numbers
mod <- feols(mean_dist ~ set | site+expgroup+as.character(treatment_year), data = comb2)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))

x <- ggpredict(mod, "set")
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  ylim(0.15,0.2)+
  ylab("Functional dispersion within sites")+
  theme_base()





#
#
#
#
#
#
#
#