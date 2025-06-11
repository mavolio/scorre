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

#Read in trait data
traits_cat <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/sCoRRE categorical trait data_12142022.csv") #categorical trait data

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

corre2trykey <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait data/corre2trykey_2021.csv") #matched species names between trait data and relative cover data
corre2trykey <- corre2trykey[,c("genus_species","species_matched")]
corre2trykey <- unique(corre2trykey)
cover <- left_join(cover, corre2trykey, by = "genus_species", keep = FALSE)

experimentinfo <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfo_March2024.csv")%>%#Information about the treatments which gets used to test how treatment magnitude explains efect sizes
  unique()

siteLocationClimate <- read.csv("C:/Users/ohler/Dropbox/CoRRE_database/Data/CompiledData/siteLocationClimate.csv") #information about sites


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
  #  left_join(nyear,by = c("site_code", "project_name", "community_type"))%>%
  # left_join(lastyear, by = c("site_code", "project_name", "community_type"))%>%
  left_join( experimentinfo, by = c("site_code", "project_name", "community_type", "treatment", "calendar_year", "treatment_year"))%>%
  left_join(repnum, by = c("site_code", "project_name", "community_type", "treatment"))

#number of sites for each year
n_sites <- crest %>%
  subset( rep_num >=5)%>%
  dplyr::select(site_code, project_name, community_type, treatment_year)%>%
  unique()%>%
  group_by(treatment_year)%>%
  dplyr::summarise(n_sites = n())


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
                                           "irr*temp" = "mult_GCD") ) #all expect for the first term are used for mult_GCD category which is no longer being used

test1 <- test%>%
  subset( trt_type == "control" | trt_type == "N" | trt_type == "P" | trt_type == "irr" |# trt_type == "drought"  | trt_type == "temp"| 
            trt_type == "mult_nutrient" #|trt_type == "mult_GCD"
          | trt_type == "CO2" | trt_type == "irr*CO2"  |trt_type == "N*irr*CO2" | trt_type == "mult_nutrient*irr" |trt_type == "N*CO2"|trt_type == "N*irr"
  ) %>% #keep only the focal treatments
  dplyr::select(site_code, project_name, community_type, treatment_year, trt_type, treatment,plot_mani, block, plot_id,species_matched, relcov)%>%
  subset(species_matched != "NA")%>%
  group_by(site_code, project_name, community_type, treatment_year, trt_type, treatment,plot_mani, block, plot_id,species_matched)%>%
  dplyr::summarize(relcover = sum(relcov))%>%
  pivot_wider(names_from = species_matched, values_from = relcover, values_fill = 0)%>%
  pivot_longer(cols=-(1:9), names_to = "species_matched", values_to = "relcover")%>%
  group_by(site_code, project_name, community_type, treatment_year, trt_type, treatment,plot_mani, species_matched)%>%
  dplyr::summarize(relcov = mean(relcover))


#test <- test[c("site_code", "project_name", "community_type", "treatment_year", "species_matched", "relcov", "trt_type", "plot_mani", "treatment")]%>%
#  unique()

#plot.treatment <- test[c("site_code", "project_name", "community_type", "trt_type", "treatment")]%>%
#  unique()
#plot.treatment <- tidyr::unite(plot.treatment, "rep", c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = FALSE)



















df <- left_join(test1, traits, by = "species_matched", keep = FALSE)

#df <- unite(df, rep, c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = FALSE)

df <- unite(df, expgroup, c("site_code", "project_name", "community_type"), sep = "::")

#a few lines to remove NAs from the continuous trait data
df$ok <- complete.cases(df[,c("SLA", "LDMC", "leaf_N", "plant_height_vegetative", "seed_dry_mass", "SRL"
)])
df <- subset(df, ok == TRUE)




######
###The same stuff with traits but they include categorical traits
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


#the below chunk collates the raw trait data into CWMs and adds the categorical data
summarize.cwm <-   
  df %>%   # First step in the next string of statements
  dplyr::group_by( expgroup, treatment_year,  trt_type, treatment, plot_mani) %>%   # Groups the summary file by expgroup
  dplyr::summarize(           # Coding for how we want our CWMs summarized
    SLA.cwm = weighted.mean(SLA, relcov),
    LDMC.cwm = weighted.mean(LDMC, relcov),
    leaf_N.cwm = weighted.mean(leaf_N, relcov),
    plant_height_vegetative.cwm = weighted.mean(plant_height_vegetative, relcov),
    seed_dry_mass.cwm = weighted.mean(seed_dry_mass, relcov),   # Actual calculation of CWMs
    SRL.cwm = weighted.mean(SRL, relcov)
  )#%>%
  #left_join(CoRRE_CWMtraits_cat, by = c("rep", "treatment_year"))


summarize.traits.continuous <- traits[,c("species_matched", "SLA", "LDMC", "leaf_N","plant_height_vegetative", "seed_dry_mass", "SRL")]
summarize.traits.continuous <- unique(summarize.traits.continuous)
#summarize.traits.categorical <- traits[,c("species_matched", "growth_form", "photosynthetic_pathway", "lifespan", "clonal", "mycorrhizal_type", "n_fixation")]
#summarize.traits.categorical <- subset(summarize.traits.categorical, photosynthetic_pathway == "C3" | photosynthetic_pathway == "C4" | photosynthetic_pathway == "CAM")

#summarize.traits <- left_join(summarize.traits.continuous, summarize.traits.categorical, by = "species_matched")                                       
# reassigning row names
#summarize.traits <- unique(summarize.traits)

trt_vector <- unique(summarize.cwm$trt_type)

#tdistances_master <- {}

for(i in 1:length(trt_vector)) {
  
  
  site.list <- subset(summarize.cwm,  trt_type == trt_vector[i])%>%dplyr::select(expgroup)%>%unique()
  
temp.df <- subset(subset(summarize.cwm,  expgroup%in%site.list$expgroup), trt_type == "control"|trt_type == trt_vector[i])
  
temp.gow <- gowdis(temp.df[6:ncol(temp.df)])
temp.beta <- betadisper(temp.gow, group = temp.df$trt_type, type = "centroid")
tdistances_temp <- data.frame(site = separate(temp.df,expgroup, into = c("site", "project", "community"), sep = "::")$site,
  expgroup = temp.df$expgroup, trt_type = temp.df$trt_type, treatment = temp.df$treatment,  dist = temp.beta$dist, plot_mani = temp.df$plot_mani, treatment_year = temp.df$treatment_year)

assign(paste0("df", trt_vector[i]),tdistances_temp)
}


mod <- feols(dist~trt_type | site+expgroup +treatment_year, data = dfN)
summary(mod)

mod <- feols(dist~trt_type | site+expgroup +treatment_year, data = dfP)
summary(mod)

mod <- feols(dist~trt_type | site+expgroup +treatment_year, data = dfmult_nutrient)
summary(mod)

mod <- feols(dist~trt_type | site+expgroup +treatment_year, data = dfirr)
summary(mod)

mod <- feols(dist~trt_type | site+expgroup +treatment_year, data = dfCO2)
summary(mod)

mod <- feols(dist~trt_type | site+expgroup +treatment_year, data = `dfN*irr`)
summary(mod)


ggplot(dfmult_nutrient, aes(trt_type, dist))+
  geom_boxplot()



###All together now!
temp.gow <- gowdis(summarize.cwm[6:ncol(summarize.cwm)])
temp.beta <- betadisper(temp.gow, group = summarize.cwm$trt_type, type = "centroid")
tdistances_temp <- data.frame(site = separate(summarize.cwm,expgroup, into = c("site", "project", "community"), sep = "::")$site,
                              expgroup = summarize.cwm$expgroup, trt_type = summarize.cwm$trt_type, treatment = summarize.cwm$treatment,  dist = temp.beta$dist, plot_mani = summarize.cwm$plot_mani, treatment_year = summarize.cwm$treatment_year)



mod <- feols(dist~plot_mani | site + expgroup + treatment_year, data = tdistances_temp)
summary(mod)

x <- ggpredict(mod, "plot_mani")
ggplot(tdistances_temp, aes(plot_mani, dist))+
  geom_point(aes(color = trt_type), alpha = 0.1)+
  geom_smooth(data=x, aes(x=x, y=predicted), se = FALSE)+
  #geom_smooth(method = "loess")+
  theme_base()

mod <- feols(dist~plot_mani*treatment_year | site + expgroup , data = tdistances_temp)
summary(mod)

x <- ggpredict(mod, c("plot_mani", "treatment_year"))
ggplot(tdistances_temp, aes(plot_mani, dist))+
  facet_wrap(~treatment_year)+
  geom_point(aes(color = trt_type), alpha = 0.1)+
  #geom_smooth(data=x, aes(x=x, y=predicted), se = FALSE)+
  geom_smooth(method = "loess")+
  theme_base()
  