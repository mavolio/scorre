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
traits_cat <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/sCoRRE categorical trait data_12142022.csv") #categorical trait data

#traits1 <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/old files/Final TRY Traits/Imputed Continuous_Traits/data to play with/imputed_continuous_20220620.csv")
traits <- read.csv("C:/Users/ohler/Downloads/CoRRE_allTraitData_wide_June2023.csv")%>%
  dplyr::select(species_matched, trait, trait_value)%>%
  pivot_wider(names_from = trait, values_from = trait_value)



corre2trykey <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/old files/corre2trykey_2021.csv") #contrinuous trait data

cover <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Jan2023.csv") %>% #community comp relative cover data
    mutate(drop=ifelse(site_code=="CDR"&treatment==2|site_code=="CDR"&treatment==3|site_code=="CDR"&treatment==4|site_code=="CDR"&treatment==5|site_code=="CDR"&treatment==7, 1,0))%>%
  filter(drop==0) #remove some Cedar Creek treatments since that site is somewhat overrepresented


corre2trykey <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/old files/corre2trykey_2021.csv") #matched species names between trait data and relative cover data
corre2trykey <- corre2trykey[,c("genus_species","species_matched")]
corre2trykey <- unique(corre2trykey)
cover <- left_join(cover, corre2trykey, by = "genus_species", all.x = TRUE)

experimentinfo <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfo_Dec2021.csv")#Information about the treatments which gets used to test how treatment magnitude explains efect sizes




#create average trait values per species and clean outliers
#traits<-traits1%>%
  #select(-X.1, -X, -family, -genus, -observation)%>%
#  select(species_matched, SLA, LDMC, leaf_N, plant_height_vegetative, seed_dry_mass, SRL) %>% #SLA, LDMC, leaf N, plant height, seed mass, and SRL
#  group_by(species_matched)%>%
#  summarise_all(funs(mean))%>%
#  ungroup() # %>%
  #filter(seed_dry_mass<30, plant_height_vegetative<10, rooting_depth<3, SLA<75)


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
df$ok <- complete.cases(df[,c("SLA", "LDMC", "leaf_N", "plant_height_vegetative", "seed_dry_mass", "SRL"
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
  mutate(con_minus_trt = dist.con-dist.trt)

lrr.df.conf <- lrr.df%>%
  ddply(.(trt_type), function(x)data.frame(
    lrr.mean = mean(x$lrr),
    lrr.error = qt(0.975, df=length(x$trt_type)-1)*sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type)-1),
    lrr.se = sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type)),
    num_experiments = length(x$expgroup)
  ))

lrr.df.conf$trt_type <- factor(lrr.df.conf$trt_type, levels = c("drought", "irr", "temp", "N", "P", "mult_nutrient"#, "mult_GCD", "CO2"
                                                                ))
lrr.df.conf$min <- lrr.df.conf$lrr.mean-lrr.df.conf$lrr.error
lrr.df.conf$max <- lrr.df.conf$lrr.mean+lrr.df.conf$lrr.error


#visualize
ggplot(lrr.df.conf, aes(trt_type, lrr.mean, color = trt_type))+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
    geom_pointrange(aes(ymin = lrr.mean-lrr.error, ymax = lrr.mean+lrr.error), size = 1.5)+
  xlab("")+
  ylab("Species composition LRR distance between plots within treatment")+
  scale_color_manual(values = c("#df0000","#0099f6", "orange", "#00b844","#f2c300","#6305dc", "black"))+
  theme_base()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


ggsave(
  "C:/Users/ohler/Documents/converge-diverge/fig1_comp.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 7,
  height = 6,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)

#models to test results
#mod <- lmer(lrr~0+ (1|expgroup), data = subset(lrr.df, trt_type == "drought" ))
#summary(mod)

mod <- lmer(lrr~0+ trt_type + (1|expgroup/trt_type), data = lrr.df)
summary(mod)




distances_master.1 <- tidyr::separate(distances_master, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)

expgroup_drought <- c(distances_master.1%>%
                        dplyr::select(expgroup, trt_type)%>%
                        subset(trt_type == "drought")%>%
                        unique())$expgroup
tempdf <- subset(distances_master.1, expgroup %in% expgroup_drought)%>%
  subset(trt_type == "control" | trt_type == "drought")
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)

tempdf <- subset(lrr.df, trt_type == "drought")
mod <- lmer(lrr~0 + (1|expgroup/trt_type), data = tempdf)
summary(mod)


expgroup_irr <- c(distances_master.1%>%
  dplyr::select(expgroup, trt_type)%>%
  subset(trt_type == "irr"))$expgroup
tempdf <- subset(distances_master.1, expgroup %in% expgroup_irr)%>%
          subset(trt_type == "control" | trt_type == "irr")
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)

expgroup_temp <- c(distances_master.1%>%
                    dplyr::select(expgroup, trt_type)%>%
                    subset(trt_type == "temp"))$expgroup
tempdf <- subset(distances_master.1, expgroup %in% expgroup_temp)%>%
  subset(trt_type == "control" | trt_type == "temp")
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)

expgroup_N <- c(distances_master.1%>%
                     dplyr::select(expgroup, trt_type)%>%
                     subset(trt_type == "N"))$expgroup
tempdf <- subset(distances_master.1, expgroup %in% expgroup_N)%>%
  subset(trt_type == "control" | trt_type == "N")
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)

expgroup_P <- c(distances_master.1%>%
                     dplyr::select(expgroup, trt_type)%>%
                     subset(trt_type == "P"))$expgroup
tempdf <- subset(distances_master.1, expgroup %in% expgroup_P)%>%
  subset(trt_type == "control" | trt_type == "P")
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)

expgroup_mult_nutrient <- c(distances_master.1%>%
                     dplyr::select(expgroup, trt_type)%>%
                     subset(trt_type == "mult_nutrient"))$expgroup
tempdf <- subset(distances_master.1, expgroup %in% expgroup_mult_nutrient)%>%
  subset(trt_type == "control" | trt_type == "mult_nutrient")
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)
#mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = subset(distances_master.1, trt_type == "control" | trt_type == "mult_GCD"))
#summary(mod)
            
            
lrr.df_species <- lrr.df


###Summarize sites being used
sites <- test%>%
  dplyr::select(site_code, project_name, community_type, treatment_year, trt_type, treatment)%>%
  unique()%>%
  subset(trt_type != "control")

n <- sites%>%
#  tidyr::unite("expgroup", c("site_code", "project_name", "community_type"))%>%
      ddply(.(trt_type), function(x)data.frame(n = length(x$site_code)))


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
    SLA.cwm = weighted.mean(SLA, relcov),
    LDMC.cwm = weighted.mean(LDMC, relcov),
    leaf_N.cwm = weighted.mean(leaf_N, relcov),
    plant_height_vegetative.cwm = weighted.mean(plant_height_vegetative, relcov),
    seed_dry_mass.cwm = weighted.mean(seed_dry_mass, relcov),   # Actual calculation of CWMs
    SRL.cwm = weighted.mean(SRL, relcov)
    )%>%
  left_join(CoRRE_CWMtraits_cat, by = c("rep", "treatment_year"))


summarize.traits.continuous <- traits[,c("species_matched", "SLA", "LDMC", "leaf_N","plant_height_vegetative", "seed_dry_mass", "SRL")]
summarize.traits.continuous <- unique(summarize.traits.continuous)
summarize.traits.categorical <- traits[,c("species_matched", "growth_form", "photosynthetic_pathway", "lifespan", "clonal", "mycorrhizal_type", "n_fixation")]
summarize.traits.categorical <- subset(summarize.traits.categorical, photosynthetic_pathway == "C3" | photosynthetic_pathway == "C4" | photosynthetic_pathway == "CAM")
                                       
summarize.traits <- left_join(summarize.traits.continuous, summarize.traits.categorical, by = "species_matched")                                       
# reassigning row names
summarize.traits <- unique(summarize.traits)




###########CALCULATE BETA DIVERSITY WITH TRAITS
expgroup_vector <- unique(summarize.cwm$expgroup)

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
    lrr.se = sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type)),
    num_experiments = length(x$expgroup)
  ))


lrr.df.conf$min <- lrr.df.conf$lrr.mean-lrr.df.conf$lrr.error
lrr.df.conf$max <- lrr.df.conf$lrr.mean+lrr.df.conf$lrr.error



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
  theme_base()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


ggsave(
  "C:/Users/ohler/Documents/converge-diverge/fig1_trait.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 7,
  height = 6,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)




lrr.df_traits <- lrr.df

#models to test results
mod <- lmer(lrr~0+ trt_type + (1|expgroup/trt_type), data = lrr.df_traits)
summary(mod)

mod <- lmer(lrr~0 + (1|expgroup), data = subset(lrr.df_traits, trt_type == "drought"))
mod <- lmer(lrr~0 + (1|expgroup), data = subset(lrr.df_traits, trt_type == "mult_nutrient"))



tdistances_master.1 <- tidyr::separate(tdistances_master, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)
tdistances_master.1 <- tidyr::separate(tdistances_master, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)

expgroup_drought <- c(tdistances_master.1%>%
                              dplyr::select(expgroup, trt_type)%>%
                              subset(trt_type == "drought"))$expgroup
tempdf <- subset(tdistances_master.1, expgroup %in% expgroup_drought)%>%
  subset(trt_type == "control" | trt_type == "drought")
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)

expgroup_irr <- c(tdistances_master.1%>%
                        dplyr::select(expgroup, trt_type)%>%
                        subset(trt_type == "irr"))$expgroup
tempdf <- subset(tdistances_master.1, expgroup %in% expgroup_irr)%>%
  subset(trt_type == "control" | trt_type == "irr")
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)

expgroup_temp <- c(tdistances_master.1%>%
                        dplyr::select(expgroup, trt_type)%>%
                        subset(trt_type == "temp"))$expgroup
tempdf <- subset(tdistances_master.1, expgroup %in% expgroup_temp)%>%
  subset(trt_type == "control" | trt_type == "temp")
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)

expgroup_N <- c(tdistances_master.1%>%
                        dplyr::select(expgroup, trt_type)%>%
                        subset(trt_type == "N"))$expgroup
tempdf <- subset(tdistances_master.1, expgroup %in% expgroup_N)%>%
  subset(trt_type == "control" | trt_type == "N")
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)

expgroup_P <- c(tdistances_master.1%>%
                        dplyr::select(expgroup, trt_type)%>%
                        subset(trt_type == "P"))$expgroup
tempdf <- subset(tdistances_master.1, expgroup %in% expgroup_P)%>%
  subset(trt_type == "control" | trt_type == "P")
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)

expgroup_mult_nutrient <- c(tdistances_master.1%>%
                        dplyr::select(expgroup, trt_type)%>%
                        subset(trt_type == "mult_nutrient"))$expgroup
tempdf <- subset(tdistances_master.1, expgroup %in% expgroup_mult_nutrient)%>%
  subset(trt_type == "control" | trt_type == "mult_nutrient")
mod <- lmer(dist~trt_type + (1|site_code/expgroup), data = tempdf)
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

library(ggpmisc)
ggplot(lrr_sp.tr, aes(lrr.species, lrr.traits))+
  facet_wrap(~trt_type)+
  geom_point()+
  geom_smooth(aes(color = trt_type),method = "lm", se = FALSE)+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
  geom_vline(xintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
  ylim(-1.1, 2.4)+
  xlim(-1.1, 2.4)+
  #stat_fit_glance(method = 'lm',
                  #method.args = list(formula = formula),
  #                geom = 'text',
  #                aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")))+
  xlab("LRR species composition beta diversity")+
  ylab("LRR trait composition beta diversity")+
  scale_color_manual(values = c("#df0000","#0099f6", "orange", "#00b844","#f2c300","#6305dc"))+
  theme_base()+
  theme(legend.position = "none")

ggsave(
  "C:/Users/ohler/Documents/converge-diverge/comp-trait_correlation.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 6,
  height = 4.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)




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
water_mani <- subset(lrr_treat_species, trt_type == "drought"|
                     trt_type == "irr")
mod <- lmer(lrr~precip + (1|expgroup) ,data = water_mani)
summary(mod)
ggplot(water_mani, aes(precip, lrr))+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
  geom_vline(xintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
  geom_point(aes(color = trt_type), size = 4)+
  scale_color_manual(values = c("#df0000", "#0099f6"))+
  ylab("LRR community composition beta diversity")+
  xlab("Precip treatment as percentage of MAP")+
  ylim(-1, 1)+
  geom_smooth(method = "lm", se = FALSE, color = "black")+
  theme_base()+
  theme(legend.position = "none")

ggsave(
  "C:/Users/ohler/Documents/converge-diverge/wate_mani.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 5,
  height = 5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)


##Nitrogen gradient
temp <- subset(lrr_treat_species, trt_type == "N")
mod <- lmer(lrr~n + (1|expgroup) ,data = temp)
summary(mod)
ggplot(temp, aes(n, lrr))+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
  geom_point(size = 4, color = "#00b844")+
  ylab("LRR species composition beta diversity")+
  xlab("N addition treatment (grams/m2)")+
  #geom_smooth(method = "lm", se = FALSE)+ #makes sense to remove the geom_smooth layer as long as it's not a significant relationship
  theme_base()+
  theme(legend.position = "none")


ggsave(
  "C:/Users/ohler/Documents/converge-diverge/N_gradient.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 5,
  height = 5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)


##traits
water_mani <- subset(lrr_treat_traits, trt_type == "drought" | trt_type == "irr")
mod <- lmer(lrr~precip + (1|expgroup) ,data = water_mani)
summary(mod)


visreg(mod, xvar = "precip", yvar = "lrr", ylab = "lrr trait beta diversity", xlab = "Precipitation treatment", gg = TRUE)+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
  theme_base()
  
N <- subset(lrr_treat_traits, trt_type == "N")
mod <- lmer(lrr~n  + (1|expgroup) ,data = N)
summary(mod)
ggplot(N, aes(n, lrr))+
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
  geom_point(size = 4, color = "#00b844")+
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


mean.diff.df <- ddply(rank_diff_master,.(expgroup, trt_type, treatment, plot_mani), function(x)data.frame( 
  richness_diff = mean(x$richness_diff),
  evenness_diff = mean(x$evenness_diff),
  rank_diff = mean(x$rank_diff),
  species_diff = mean(x$species_diff)
  ))

RAC_trt.df <- subset(mean.diff.df, plot_mani >= 1)

RAC_con.df <- subset(mean.diff.df, plot_mani == 0)%>%
  dplyr::rename(rich.con = richness_diff, eve.con = evenness_diff, rank.con = rank_diff, sp.con = species_diff)%>%
  dplyr::select(expgroup, rich.con, eve.con, rank.con, sp.con)

RAC_lrr.df <- merge(RAC_trt.df, RAC_con.df, by = "expgroup", all.x = TRUE)%>%
  mutate(lrr.rich = log(richness_diff/rich.con), lrr.eve = log(evenness_diff/eve.con), lrr.rank = log(rank_diff/rank.con), lrr.sp = log(species_diff/sp.con))%>%
  mutate(sub.rich = richness_diff-rich.con, sub.eve = evenness_diff-eve.con, sub.rank = rank_diff-rank.con, sub.sp = species_diff-sp.con)


full_lrr.df <- left_join(lrr.df, RAC_lrr.df, by = c("expgroup", "trt_type", "treatment", "plot_mani"))%>%
  tidyr::separate( expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)


pairs(~ sub.rich + sub.eve + sub.rank + sub.sp, data = full_lrr.df)
c("sub.rich","sub.eve" , "sub.rank" , "sub.sp")

ggsave(
  "C:/Users/ohler/Documents/converge-diverge/RAC_var_supp.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 10,
  height = 10,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)


tempdf <- subset(full_lrr.df, trt_type == "drought" | trt_type == "N" | trt_type == "mult_nutrient")
ggplot(tempdf, aes(sub.rich, lrr))+
  facet_wrap(~trt_type)+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab("LRR distance among plots between treatments")+
  xlab("Richness difference among plots (treatment-control")+
  theme_base()

ggplot(tempdf, aes(sub.eve, lrr))+
  facet_wrap(~trt_type)+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab("LRR distance among plots between treatments")+
  xlab("Evenness difference among plots (treatment-control")+
  theme_base()

ggplot(tempdf, aes(sub.rank, lrr))+
  facet_wrap(~trt_type)+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab("LRR distance among plots between treatments")+
  xlab("Rank difference among plots (treatment-control")+
  theme_base()

ggplot(tempdf, aes(sub.sp, lrr))+
  facet_wrap(~trt_type)+
  geom_point()+
  geom_smooth(method = "lm")+
  ylab("LRR distance among plots between treatments")+
  xlab("Species difference among plots (treatment-control")+
  theme_base()


library(remotes)
remotes::install_github("mastoffel/partR2") 
library(partR2)
library(tibble)

tempdf <- subset(full_lrr.df, trt_type == "drought")
mod <- lmer(lrr~sub.rich+sub.eve + sub.rank + sub.sp + (1|expgroup), data = tempdf)
summary(mod) #sig include marginal sub.rank
r2 <- partR2(mod, data = tempdf, partvars = c("sub.rich","sub.eve" , "sub.rank" , "sub.sp"), R2_type = "marginal", nboot = 10)
r2

r2$R2%>%
  subset(term == "sub.rich" | term == "sub.eve" | term == "sub.rank" | term == "sub.sp")%>%
  mutate(term = factor(term, levels = c("sub.eve","sub.rich",  "sub.sp", "sub.rank"))) %>%
  tibble::add_column(sig = c("0","0", "0", "0"))%>%
  ggplot( aes(term, estimate, fill = sig))+
  geom_bar(color = "black",stat = "identity")+
  ylim(0,.2)+
  ylab("Partial r-squared")+
  scale_fill_manual(values = c( "white", "black"))+
  coord_flip()+
  ggtitle("DROUGHT")+
  xlab("")+
  theme_classic()+
  theme(legend.position = "none")

ggsave(
  "C:/Users/ohler/Documents/converge-diverge/partialR2_drought.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 3,
  height = 10,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)


tempdf <- tempdf <- subset(full_lrr.df, trt_type == "N")%>%
  dplyr::select(lrr, sub.rich, sub.eve, sub.rank, sub.sp, expgroup)%>%
  filter(complete.cases(.))
mod <- lmer(lrr~sub.rich+sub.eve + sub.rank + sub.sp + (1|expgroup), data = tempdf)
summary(mod) #sig include sub.rank, sub.sp
r2 <- partR2(mod, data = tempdf, partvars = c("sub.rich","sub.eve" , "sub.rank" , "sub.sp"), R2_type = "marginal", nboot = 10)
r2

r2$R2%>%
  subset(term == "sub.rich" | term == "sub.eve" | term == "sub.rank" | term == "sub.sp")%>%
  mutate(term = factor(term, levels = c("sub.eve","sub.rich",  "sub.sp", "sub.rank"))) %>%
  tibble::add_column(sig = c("0","0", "1", "0"))%>%
  ggplot( aes(term, estimate, fill = sig))+
  geom_bar(color = "black",stat = "identity")+
  ylim(0,.2)+
  ylab("Partial r-squared")+
  scale_fill_manual(values = c( "white", "black"))+
  coord_flip()+
  ggtitle("NITROGEN")+
  xlab("")+
  theme_classic()+
  theme(legend.position = "none")

ggsave(
  "C:/Users/ohler/Documents/converge-diverge/partialR2_N.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 3,
  height = 10,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)



tempdf <- subset(full_lrr.df, trt_type == "mult_nutrient")%>%
  dplyr::select(lrr, sub.rich, sub.eve, sub.rank, sub.sp, expgroup)%>%
  filter(complete.cases(.))
mod <- lmer(lrr~sub.rich+sub.eve + sub.rank + sub.sp + (1|expgroup), data = tempdf)
summary(mod) #sig include sub.rank, sub.eve
r2 <- partR2(mod, data = tempdf, partvars = c("sub.rich","sub.eve" , "sub.rank" , "sub.sp"), R2_type = "marginal", nboot = 10)
r2

r2$R2%>%
  subset(term == "sub.rich" | term == "sub.eve" | term == "sub.rank" | term == "sub.sp")%>%
  mutate(term = factor(term, levels = c("sub.sp","sub.eve",  "sub.rich", "sub.rank"))) %>%
  tibble::add_column(sig = c("0","0", "1", "0"))%>%
  ggplot( aes(term, estimate, fill = sig))+
  geom_bar(color = "black",stat = "identity")+
  ylim(0,.2)+
  ylab("Partial r-squared")+
  scale_fill_manual(values = c( "white", "black"))+
  coord_flip()+
  ggtitle("MULTIPLE NUTRIENT")+
  xlab("")+
  theme_classic()+
  theme(legend.position = "none")

ggsave(
  "C:/Users/ohler/Documents/converge-diverge/partialR2_mult_nutrient.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 3,
  height = 10,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)




#models to test results
rank_diff_master.1 <- tidyr::separate(rank_diff_master, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)


##rank_diff
tempdf <- subset(rank_diff_master.1, trt_type == "control" | trt_type == "drought")
mod <- lmer(rank_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
rank_drought <- data.frame(trt_type = "drought", metric = "rank_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "irr")
mod <- lmer(rank_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
rank_irr <- data.frame(trt_type = "irr", metric = "rank_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "temp")
mod <- lmer(rank_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
rank_temp <- data.frame(trt_type = "temp", metric = "rank_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "N")
mod <- lmer(rank_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
rank_N <- data.frame(trt_type = "N", metric = "rank_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "P")
mod <- lmer(rank_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
rank_P <- data.frame(trt_type = "P", metric = "rank_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "mult_nutrient")
mod <- lmer(rank_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
rank_mult_nutrient <- data.frame(trt_type = "mult_nutrient", metric = "rank_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

#species_diff
tempdf <- subset(rank_diff_master.1, trt_type == "control" | trt_type == "drought")
mod <- lmer(species_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
species_drought <- data.frame(trt_type = "drought", metric = "species_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "irr")
mod <- lmer(species_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
species_irr <- data.frame(trt_type = "irr", metric = "species_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "temp")
mod <- lmer(species_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
species_temp <- data.frame(trt_type = "temp", metric = "species_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "N")
mod <- lmer(species_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
species_N <- data.frame(trt_type = "N", metric = "species_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "P")
mod <- lmer(species_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
species_P <- data.frame(trt_type = "P", metric = "species_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "mult_nutrient")
mod <- lmer(species_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)                           
species_mult_nutrient <- data.frame(trt_type = "mult_nutrient", metric = "species_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

#richness_diff
tempdf <- subset(rank_diff_master.1, trt_type == "control" | trt_type == "drought")
mod <- lmer(richness_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
richness_drought <- data.frame(trt_type = "drought", metric = "richness_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "irr")
mod <- lmer(richness_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
richness_irr <- data.frame(trt_type = "irr", metric = "richness_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "temp")
mod <- lmer(richness_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
richness_temp <- data.frame(trt_type = "temp", metric = "richness_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "N")
mod <- lmer(richness_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
richness_N <- data.frame(trt_type = "N", metric = "richness_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "P")
mod <- lmer(richness_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
richness_P <- data.frame(trt_type = "P", metric = "richness_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "mult_nutrient")
mod <- lmer(richness_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)                           
richness_mult_nutrient <- data.frame(trt_type = "mult_nutrient", metric = "richness_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])



#evenness_diff
tempdf <- subset(rank_diff_master.1, trt_type == "control" | trt_type == "drought")
mod <- lmer(evenness_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
evenness_drought <- data.frame(trt_type = "drought", metric = "evenness_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "irr")
mod <- lmer(evenness_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
evenness_irr <- data.frame(trt_type = "irr", metric = "evenness_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "temp")
mod <- lmer(evenness_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
evenness_temp <- data.frame(trt_type = "temp", metric = "evenness_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "N")
mod <- lmer(evenness_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
evenness_N <- data.frame(trt_type = "N", metric = "evenness_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "P")
mod <- lmer(evenness_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)
evenness_P <- data.frame(trt_type = "P", metric = "evenness_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])

tempdf <-subset(rank_diff_master.1, trt_type == "control" | trt_type == "mult_nutrient")
mod <- lmer(evenness_diff~trt_type + (1|expgroup), data = tempdf)
summary(mod)                           
evenness_mult_nutrient <- data.frame(trt_type = "mult_nutrient", metric = "evenness_diff", Estimate = summary(mod)$coefficients[2,1], se =  summary(mod)$coefficients[2,2], pvalue = summary(mod)$coefficients[2,5])



diff.df <- dplyr::bind_rows(rank_drought, rank_irr, rank_temp, rank_N, rank_P, rank_mult_nutrient,
                 species_drought, species_irr, species_temp, species_N, species_P, species_mult_nutrient,
                 richness_drought, richness_irr, richness_temp, richness_N, richness_P, richness_mult_nutrient,
                 evenness_drought, evenness_irr, evenness_temp, evenness_N, evenness_P, evenness_mult_nutrient)

diff.df$trt_type <- factor(diff.df$trt_type, levels = c("drought", "irr", "temp", "N", "P", "mult_nutrient" ))

###Make bar graph to summarize the above information
diff.df$Estimate.sig <- ifelse(diff.df$pvalue <0.05, diff.df$metric, 0)
ggplot(diff.df, aes(x=trt_type, y=Estimate, fill = metric, color = metric))+
  geom_bar(aes(fill = ifelse(pvalue <0.05, metric, "NA"), color =  metric),stat="identity",
           position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c("#ff0000", "white", "#0000ff", "#ffa500", "#4b0082"))+
  scale_color_manual(values = c("#ff0000","#0000ff", "#ffa500", "#4b0082"))+
      #geom_errorbar(aes(ymin=ifelse(pvalue <0.05, Estimate-se, NA), ymax=ifelse(pvalue <0.05, Estimate+se, NA)), width=1
       #           , position=position_dodge(width = 0.9), color = "black"
        #        )+
  theme_base()


ggplot(diff.df, aes(x=trt_type, y=Estimate, fill = metric, color = metric))+
  facet_wrap(~metric)+
  geom_bar(aes(fill = ifelse(pvalue <0.05, metric, "NA"), color =  metric),stat="identity",
           position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c("#ff0000", "white", "#0000ff", "#ffa500", "#4b0082"))+
  scale_color_manual(values = c("#ff0000","#0000ff", "#ffa500", "#4b0082"))+
  theme_base()

ggplot(diff.df, aes(x=metric, y=Estimate, fill = metric, color = metric))+
  facet_wrap(~trt_type)+
  geom_bar(aes(fill = ifelse(pvalue <0.05, metric, "NA"), color =  metric),stat="identity",
           position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c("#ff0000", "white", "#0000ff", "#ffa500", "#4b0082"))+
  scale_color_manual(values = c("#ff0000","#0000ff", "#ffa500", "#4b0082"))+
  theme_base()


################BLEOW SECTION IS A TEST THAT YOU CAN PROBABLY DELETE
##Explain variation using partial r squareds

rank_diff_expgroup <- ddply(rank_diff_master.1, .(expgroup, trt_type, treatment, plot_mani),function(x)data.frame(
                        richness_diff = mean(x$richness_diff),
                          evenness_diff = mean(x$evenness_diff),
                          rank_diff = mean(x$rank_diff),
                          species_diff = mean(x$species_diff)
))


distances_expgroup <- ddply(distances_master.1, .(expgroup, site_code, project, community, trt_type, treatment), function(x)data.frame(
                            dist = mean(x$dist)
))


full_df <- left_join(distances_expgroup, rank_diff_expgroup, by = c("expgroup", "trt_type", "treatment"))

library(remotes)
remotes::install_github("mastoffel/partR2") 
library(partR2)

tempdf <- subset(full_df, trt_type == "control" | trt_type == "drought")%>% filter(complete.cases(.))
mod <- lmer(dist~richness_diff+evenness_diff + rank_diff + species_diff + (1|site_code), data = tempdf)
summary(mod)
partR2(mod, data = tempdf, partvars = c("richness_diff", "evenness_diff", "rank_diff", "species_diff"), R2_type = "marginal", nboot = 10)



tempdf <- subset(full_df, trt_type == "control" | trt_type == "N")%>% filter(complete.cases(.))
mod <- lmer(dist~richness_diff+evenness_diff + rank_diff + species_diff + (1|site_code), data = tempdf)
summary(mod)
partR2(mod, data = tempdf, partvars = c("richness_diff", "evenness_diff", "rank_diff", "species_diff"), R2_type = "marginal", nboot = 10)

#######################################above


####SINGLE TRAIT VARIANCE AMONG REPLICATES (TRAIT)

trait_variance <- summarize.cwm%>%
    ddply(.(expgroup, trt_type, treatment, plot_mani), function(x)data.frame(
     seed_dry_mass.var = sd(x$seed_dry_mass.cwm)/mean(x$seed_dry_mass.cwm),
    LDMC.var =  sd(x$LDMC.cwm)/mean(x$LDMC.cwm),
    plant_height_vegetative.var =  sd(x$plant_height_vegetative.cwm)/mean(x$plant_height_vegetative.cwm),
    SLA.var =  sd(x$SLA.cwm)/mean(x$plant_height_vegetative.cwm),
    rooting_depth.var = sd(x$rooting_depth.cwm)/mean(x$rooting_depth.cwm),
    leaf_C.N.var =  sd(x$leaf_C.N.cwm)/mean(x$leaf_C.N.cwm)
    ))


trt.df <- subset(trait_variance, plot_mani >= 1)%>%
  dplyr::rename(c(seed_dry_mass.var.trt = seed_dry_mass.var, 
                  LDMC.var.trt = LDMC.var, 
                  plant_height_vegetative.var.trt = plant_height_vegetative.var, 
                  SLA.var.trt = SLA.var, 
                  rooting_depth.var.trt = rooting_depth.var,
                  leaf_C.N.var.trt = leaf_C.N.var))
con.df <- subset(trait_variance, plot_mani == 0)%>%
  dplyr::rename(c(seed_dry_mass.var.con = seed_dry_mass.var, 
                  LDMC.var.con = LDMC.var, 
                  plant_height_vegetative.var.con = plant_height_vegetative.var, 
                  SLA.var.con = SLA.var, 
                  rooting_depth.var.con = rooting_depth.var,
                  leaf_C.N.var.con = leaf_C.N.var))%>%
  dplyr::select(!c(trt_type, treatment, plot_mani))

var.summary <- merge(trt.df, con.df, by = "expgroup", all.x = TRUE)

trait_variance <- tidyr::separate(trait_variance, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)

mod <- lmer(seed_dry_mass.var~trt_type + (1|site_code/expgroup), data = trait_variance)
summary(mod)

mod <- lmer(LDMC.var~trt_type + (1|expgroup), data = trait_variance)
summary(mod)

mod <- lmer(plant_height_vegetative.var~trt_type + (1|expgroup), data = trait_variance)
summary(mod)

mod <- lmer(SLA.var~trt_type + (1|expgroup), data = trait_variance)
summary(mod)

mod <- lmer(rooting_depth.var~trt_type + (1|expgroup), data = trait_variance)
summary(mod)

mod <- lmer(leaf_C.N.var~trt_type + (1|expgroup), data = trait_variance)
summary(mod)


