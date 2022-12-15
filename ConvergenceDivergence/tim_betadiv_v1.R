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


#Read in data
#traits <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/trait data/Final Cleaned Traits/Continuous_Traits/Backtrans_GapFilled_sCorre.csv")

traits_cat <- read.csv("C:/Users/Timothy/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/sCoRRE categorical trait data_11302021.csv")

traits1 <- read.csv("C:/Users/Timothy/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/Final TRY Traits/Imputed Continuous_Traits/data to play with/imputed_continuous_20220620.csv")
corre2trykey <- read.csv("C:/Users/Timothy/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/corre2trykey_2021.csv")

cover <- read.csv("C:/Users/Timothy/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Dec2021.csv") %>% 
    mutate(drop=ifelse(site_code=="CDR"&treatment==2|site_code=="CDR"&treatment==3|site_code=="CDR"&treatment==4|site_code=="CDR"&treatment==5|site_code=="CDR"&treatment==7, 1,0))%>%
  filter(drop==0)


corre2trykey <- read.csv("C:/Users/Timothy/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/corre2trykey_2021.csv")
corre2trykey <- corre2trykey[,c("genus_species","species_matched")]
corre2trykey <- unique(corre2trykey)
cover <- left_join(cover, corre2trykey, by = "genus_species", all.x = TRUE)

experimentinfo <- read.csv("C:/Users/Timothy/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfo_Dec2021.csv")

treatment_info <- read.csv("C:/Users/Timothy/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/basic dataset info/ExperimentInfo.csv")


#create average trait values per species and clean outliers
traits<-traits1%>%
  select(-X.1, -X, -family, -genus, -observation)%>%
  select(species_matched, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass) %>%
  group_by(species_matched)%>%
  summarise_all(funs(mean))%>%
  ungroup()  %>%
  filter(seed_dry_mass<30, plant_height_vegetative<10, rooting_depth<3, SLA<75)


#standardize the scale of all the traits
cols <- c( "seed_dry_mass", 
           # "stem_spec_density",
           #"leaf_N",
           #"leaf_P",
           "LDMC",
           #"leaf_C",
           #"leaf_dry_mass",
           "plant_height_vegetative",
           #"leaf_C.N",
           "SLA",
           #"water_content",
           "rooting_depth"#,
           #"SRL",
           #"seed_number"
)

traits[cols] <- scale(traits[cols])

#merge categorical traits

traits <- left_join(traits, traits_cat, by = "species_matched")


#Reduce cover data to focal data using a series of merges


#minimum number of replicates
repnum <- cover%>%
  dplyr::select(site_code, project_name, community_type, treatment, plot_id)%>%
  unique()%>%
  #dplyr::mutate(present = 1)%>%
  dplyr::group_by(site_code, project_name, community_type, treatment)%>%
  dplyr::summarise(rep_num = length(plot_id))%>%
  dplyr::ungroup()


#only certain manipulations
foctrt <- experimentinfo[c("site_code", "project_name", "community_type", "trt_type")]%>%
  unique()%>%
  subset( trt_type == "N" | trt_type == "P" | trt_type == "irr" | trt_type == "drought"| trt_type == "N*P" | trt_type == "mult_nutrient")

foctrt$treats_wanted <- foctrt$trt_type
foctrt <- foctrt[c("site_code", "project_name", "community_type", "treats_wanted")]



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
#  subset(last_trt_yr == calendar_year)%>%
  subset(treatment_year == n.trt.yrs)%>%

  subset( trt_type == "control" | trt_type == "N" | trt_type == "P" | trt_type == "irr" | 
            trt_type == "drought" | trt_type == "N*P" | trt_type == "mult_nutrient"
  )#%>%

test$trt_type <- revalue(test$trt_type, c("N*P" = "mult_nutrient"))

N <-  subset(test[test$trt_type %in% "N",], n.trt.yrs >= 5)
P <-  subset(test[test$trt_type %in% "P",], n.trt.yrs >= 5)
irr <-  subset(test[test$trt_type %in% "irr",], n.trt.yrs >= 5)
mult_nutrient <-  subset(test[test$trt_type %in% "mult_nutrient",], n.trt.yrs >= 5)
drought <-  subset(test[test$trt_type %in% "drought",], n.trt.yrs >= 4)
control <-  subset(test[test$trt_type %in% "control",], n.trt.yrs >= 4)

test <- bind_rows(N, P, irr, mult_nutrient, drought, control)

test <- test[c("site_code", "project_name", "community_type", "treatment_year", "plot_id", "species_matched", "relcov", "trt_type", "plot_mani", "treatment")]%>%
  unique()


plot.treatment <- test[c("site_code", "project_name", "community_type", "plot_id", "trt_type", "treatment")]%>%
  unique()
plot.treatment <- tidyr::unite(plot.treatment, "rep", c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = FALSE)

df <- left_join(test, traits, by = "species_matched", all.x = TRUE)

df <- unite(df, rep, c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = FALSE)

df <- unite(df, expgroup, c("site_code", "project_name", "community_type"), sep = "::")

df$ok <- complete.cases(df[,c(#"seed_dry_mass",
  #"stem_spec_density",
  #"leaf_N",
  #"leaf_P",
  "LDMC",
  "plant_height_vegetative",
  "SLA",
  "rooting_depth",
  "seed_dry_mass"
  #"leaf_C.N"
)])
df <- subset(df, ok == TRUE)
df <- subset(df, species_matched != "NA")

length(unique(df$species_matched))

########################
##Summarize sites being used

sites <- test%>%
  dplyr::select(site_code, project_name, community_type, treatment_year, trt_type, treatment)%>%
  unique()%>%
  subset(trt_type != "control")

sites <- unite(sites, temp, c("project_name", "community_type"), sep = "::", remove = FALSE)



##############################
####



#'test' dataframe has all the cover data but only with focal sites and treatments and such

#For each treatment at each site, pull the treatment and control data, spread, calculate distance matrix, then betadisper 

kevin <- unite(test, expgroup, c("site_code", "project_name", "community_type"), remove = FALSE, sep = "::" )
kevin <- subset(kevin, species_matched != "NA")%>%
  ddply(.(expgroup, site_code, project_name, community_type, treatment_year, plot_id, species_matched, trt_type, plot_mani, treatment),
        function(x)data.frame(
          relcov = sum(x$relcov)
        ))

expgroup_vector <- unique(kevin$expgroup)

distances_master <- {}

for(i in 1:length(expgroup_vector)) {
  temp.df <- subset(kevin, expgroup == expgroup_vector[i])
  
  temp.wide <- temp.df%>%
    pivot_wider(names_from = species_matched, values_from = relcov, values_fill = 0)
  temp.distances <- vegdist(temp.wide[10:ncol(temp.wide)], method = "bray")
  temp.mod <- betadisper(temp.distances, group = temp.wide$trt_type, type = "centroid")
  distances_temp <- data.frame(expgroup = expgroup_vector[i], trt_type = temp.wide$trt_type, treatment = temp.wide$treatment, plot_mani = temp.wide$plot_mani, dist = temp.mod$dist)
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

lrr.df.conf$trt_type <- factor(lrr.df.conf$trt_type, levels = c("drought", "irr", "N", "P", "mult_nutrient"))

ggplot(lrr.df.conf, aes(trt_type, lrr.mean, color = trt_type))+
  #geom_point()+
  #geom_errorbar(aes(ymin = lrr.mean-lrr.error, ymax = lrr.mean+lrr.error))+
  geom_pointrange(aes(ymin = lrr.mean-lrr.error, ymax = lrr.mean+lrr.error), size = 1.5)+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  xlab("")+
  ylab("Species composition LRR distance between plots within treatment")+
  scale_color_manual(values = c("#df0000","#0099f6","#00b844","#f2c300","#6305dc"))+
  theme_bw()

ggplot(lrr.df_species, aes(trt_type, lrr, color = trt_type))+
    geom_beeswarm(cex = 2)+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  xlab("")+
  ylab("Species composition LRR distance between plots within treatment")+
  scale_color_manual(values = c("#df0000","#0099f6","#00b844","#f2c300","#6305dc"))+
  theme_bw()

distances_master.1 <- tidyr::separate(distances_master, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(distances_master.1, trt_type == "control" | trt_type == "drought"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(distances_master.1, trt_type == "control" | trt_type == "irr"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(distances_master.1, trt_type == "control" | trt_type == "N"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(distances_master.1, trt_type == "control" | trt_type == "P"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(distances_master.1, trt_type == "control" | trt_type == "mult_nutrient"))
summary(mod)
            
            
lrr.df_species <- lrr.df


###Summarize sites being used
sites <- test%>%
  dplyr::select(site_code, project_name, community_type, treatment_year, trt_type, treatment)%>%
  unique()%>%
  subset(trt_type != "control")




######
###Try the same stuff with traits but they include categorical traits
CoRRE_CWMtraits <- read.csv("C:/Users/Timothy/Dropbox/sDiv_sCoRRE_shared/paper 2_PD and FD responses/data/CoRRE_CWMtraits_12142022.csv") #for now I'll just use this for categorical traits
CoRRE_CWMtraits_cat <- CoRRE_CWMtraits[, c(   "site_code", "project_name","community_type", "plot_id", "treatment_year", "CWM.growth_form", "CWM.photosynthetic_pathway", "CWM.lifespan", "CWM.clonal", "CWM.mycorrhizal_type", "CWM.n_fixation")]

CoRRE_CWMtraits_cat <- tidyr::unite(CoRRE_CWMtraits_cat, "rep", c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = TRUE)

CoRRE_CWMtraits_cat$is.graminoid <- ifelse(CoRRE_CWMtraits_cat$CWM.growth_form == "graminoid", 1, 0)
CoRRE_CWMtraits_cat$is.C4 <- ifelse(CoRRE_CWMtraits_cat$CWM.photosynthetic_pathway == "C$", 1, 0)
CoRRE_CWMtraits_cat$is.perennial <- ifelse(CoRRE_CWMtraits_cat$CWM.growth_form == "perennial", 1, 0)
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
    plant_vegetative_height.cwm = weighted.mean(plant_height_vegetative, relcov),
    SLA.cwm = weighted.mean(SLA, relcov),
    rooting_depth.cwm = weighted.mean(rooting_depth, relcov)
    )%>%
  left_join(CoRRE_CWMtraits_cat, by = c("rep", "treatment_year"))



expgroup_vector <- unique(df$expgroup)

tdistances_master <- {}

for(i in 1:length(expgroup_vector)) {
  temp.df <- subset(summarize.cwm, expgroup == expgroup_vector[i])

  
  temp.gow <- gowdis(temp.df[7:ncol(temp.df)])
  temp.beta <- betadisper(temp.gow, group = temp.df$trt_type, type = "centroid")
  
  tdistances_temp <- data.frame(expgroup = expgroup_vector[i], trt_type = temp.df$trt_type, treatment = temp.df$treatment,  dist = temp.beta$dist, plot_mani = temp.df$plot_mani)
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

lrr.df.conf$trt_type <- factor(lrr.df.conf$trt_type, levels = c("drought", "irr", "N", "P", "mult_nutrient"))

ggplot(lrr.df.conf, aes(trt_type, lrr.mean, color = trt_type))+
  geom_pointrange(aes(ymin = lrr.mean-lrr.error, ymax = lrr.mean+lrr.error), size = 1.5)+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  xlab("")+
  ylab("Trait LRR distance between plots within treatment")+
  scale_color_manual(values = c("#df0000","#0099f6","#00b844","#f2c300","#6305dc"))+
  theme_bw()

ggplot(lrr.df_traits, aes(trt_type, lrr, color = trt_type))+
  geom_beeswarm(cex = 2)+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  xlab("")+
  ylab("Trait composition LRR distance between plots within treatment")+
  scale_color_manual(values = c("#df0000","#0099f6","#00b844","#f2c300","#6305dc"))+
  theme_bw()

lrr.df_traits <- lrr.df

tdistances_master.1 <- tidyr::separate(tdistances_master, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(tdistances_master.1, trt_type == "control" | trt_type == "drought"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(tdistances_master.1, trt_type == "control" | trt_type == "irr"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(tdistances_master.1, trt_type == "control" | trt_type == "N"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(tdistances_master.1, trt_type == "control" | trt_type == "P"))
summary(mod)
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(tdistances_master.1, trt_type == "control" | trt_type == "mult_nutrient"))
summary(mod)







############
##BRING IN COVARIATES AND SEE IF THEY EXPLAIN BETA DIVERSITY


CoRRE_siteLocationClimate_Dec2021 <- read.csv("C:/Users/Timothy/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/environmental data/CoRRE_siteLocationClimate_Dec2021.csv")

CoRRE_project_summary <- read.csv("C:/Users/Timothy/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/CoRRE_project_summary.csv")
CoRRE_project_summary$project <- CoRRE_project_summary$project_name
CoRRE_project_summary$community <- CoRRE_project_summary$community_type
CoRRE_project_summary <- CoRRE_project_summary %>% dplyr::select(-c(project_name, community_type))


lrr.df_species <- tidyr::separate(lrr.df_species, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)
lrr.df_traits <- tidyr::separate(lrr.df_traits, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)

lrr_covariate <- left_join(lrr.df_species, CoRRE_project_summary, by = c("site_code", "project", "community"))
lrr_covariate_traits <- left_join(lrr.df_traits, CoRRE_project_summary, by = c("site_code", "project", "community"))

mod <- lmer(lrr~MAP + (1|site_code), data = subset(lrr_covariate, trt_type == "drought"))
summary(mod)
mod <- lmer(lrr~MAT + (1|site_code), data = subset(lrr_covariate, trt_type == "drought"))
summary(mod)
mod <- lmer(lrr~rrich + (1|site_code), data = subset(lrr_covariate, trt_type == "drought"))
summary(mod)
mod <- lmer(lrr~experiment_length + (1|site_code), data = subset(lrr_covariate, trt_type == "drought"))
summary(mod)

mod <- lmer(lrr~MAP + (1|site_code), data = subset(lrr_covariate, trt_type == "N"))
summary(mod)
mod <- lmer(lrr~MAT + (1|site_code), data = subset(lrr_covariate, trt_type == "N"))
summary(mod)
mod <- lmer(lrr~rrich + (1|site_code), data = subset(lrr_covariate, trt_type == "N"))
summary(mod)
mod <- lmer(lrr~experiment_length + (1|site_code), data = subset(lrr_covariate, trt_type == "N"))
summary(mod)

mod <- lmer(lrr~MAP + (1|site_code), data = subset(lrr_covariate, trt_type == "mult_nutrient"))
summary(mod)
mod <- lmer(lrr~MAT + (1|site_code), data = subset(lrr_covariate, trt_type == "mult_nutrient"))
summary(mod)
mod <- lmer(lrr~rrich + (1|site_code), data = subset(lrr_covariate, trt_type == "mult_nutrient"))
summary(mod)
mod <- lmer(lrr~experiment_length + (1|site_code), data = subset(lrr_covariate, trt_type == "mult_nutrient"))
summary(mod)




###trait models
mod <- lmer(lrr~MAP + (1|site_code), data = subset(lrr_covariate_traits, trt_type == "drought"))
summary(mod)
mod <- lmer(lrr~MAT + (1|site_code), data = subset(lrr_covariate_traits, trt_type == "drought"))
summary(mod)
mod <- lmer(lrr~rrich + (1|site_code), data = subset(lrr_covariate_traits, trt_type == "drought"))
summary(mod)
mod <- lmer(lrr~experiment_length + (1|site_code), data = subset(lrr_covariate_traits, trt_type == "drought"))
summary(mod)

mod <- lmer(lrr~MAP + (1|site_code), data = subset(lrr_covariate_traits, trt_type == "N"))
summary(mod)
mod <- lmer(lrr~MAT + (1|site_code), data = subset(lrr_covariate_traits, trt_type == "N"))
summary(mod)
mod <- lmer(lrr~rrich + (1|site_code), data = subset(lrr_covariate_traits, trt_type == "N"))
summary(mod)
mod <- lmer(lrr~experiment_length + (1|site_code), data = subset(lrr_covariate_traits, trt_type == "N"))
summary(mod)

mod <- lmer(lrr~MAP + (1|site_code), data = subset(lrr_covariate_traits, trt_type == "mult_nutrient"))
summary(mod)
mod <- lmer(lrr~MAT + (1|site_code), data = subset(lrr_covariate_traits, trt_type == "mult_nutrient"))
summary(mod)
mod <- lmer(lrr~rrich + (1|site_code), data = subset(lrr_covariate_traits, trt_type == "mult_nutrient"))
summary(mod)
mod <- lmer(lrr~experiment_length + (1|site_code), data = subset(lrr_covariate_traits, trt_type == "mult_nutrient"))
summary(mod)


##with treatment information
treatment_info$trt_type <- revalue(treatment_info$trt_type, c("N*P" = "mult_nutrient"))
treatment_info$project <- treatment_info$project_name
treatment_info$community <- treatment_info$community_type
treatment_info <- treatment_info[,c("site_code", "project", "community", "plot_mani", "trt_type", "treatment", "n", "p", "k", "CO2", "precip", "temp")]
treatment_info <- unique(treatment_info)
lrr_treat_species <- left_join(lrr_covariate, treatment_info, by = c("site_code", "project", "community", "plot_mani", "trt_type", "treatment"))


mod <- lmer(lrr~precip + (1|site_code) ,data = subset(lrr_treat_species, trt_type == "drought"))
summary(mod)
mod <- lmer(lrr~n + (1|site_code) ,data = subset(lrr_treat_species, trt_type == "N"))
summary(mod)
