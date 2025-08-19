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
traits_cat <- read.csv('https://pasta.lternet.edu/package/data/eml/edi/1533/3/5ebbc389897a6a65dd0865094a8d0ffd')%>%
  dplyr::select(-family, -source, -error_risk_overall)%>%
  pivot_wider(names_from = trait, values_from = trait_value)%>%
  dplyr::rename(species_matched = species)
#categorical trait data

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



df <- left_join(test1, traits, by = "species_matched", keep = FALSE)

df <- unite(df, expgroup, c("site_code", "project_name", "community_type"), sep = "::")

#a few lines to remove NAs from the continuous trait data
df$ok <- complete.cases(df[,c("SLA", "LDMC", "leaf_N", "plant_height_vegetative", "seed_dry_mass", "SRL"
)])
df <- subset(df, ok == TRUE)


######
###The same stuff with traits but they include categorical traits
df$growth_form <- ifelse(df$growth_form == "graminoid", 1, 0)
df$photosynthetic_pathway <- ifelse(df$photosynthetic_pathway == "C4", 1, 0)
df$lifespan <- ifelse(df$lifespan == "perennial", 1, 0)
df$clonal <- ifelse(df$clonal == "yes", 1, 0)
df$mycorrhizal_type <- ifelse(df$mycorrhizal_type == "AM", 1, 0)
df$n_fixation_type <- ifelse(df$n_fixation_type == "rhizobial", 1, 0)


#the below chunk collates the raw trait data into CWMs and adds the categorical data
summarize.cwm <-   
  df %>%   # First step in the next string of statements
  dplyr::group_by(expgroup, treatment_year, trt_type, treatment, plot_mani) %>%   # Groups the summary file by Plot number
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


#loop
trt_vector <- unique(summarize.cwm$trt_type)

for(i in 1:length(trt_vector)) {
  
  
  site.list <- subset(summarize.cwm,  trt_type == trt_vector[i])%>%dplyr::select(expgroup)%>%unique()
  
temp.df <- subset(subset(summarize.cwm,  expgroup%in%site.list$expgroup), trt_type == "control"|trt_type == trt_vector[i])
  
temp.gow <- gowdis(temp.df[6:ncol(temp.df)])
temp.beta <- betadisper(temp.gow, group = temp.df$trt_type, type = "centroid")
tdistances_temp <- data.frame(site = separate(temp.df,expgroup, into = c("site", "project", "community"), sep = "::")$site,
  expgroup = temp.df$expgroup, trt_type = temp.df$trt_type, treatment = temp.df$treatment,  dist = temp.beta$dist, plot_mani = temp.df$plot_mani, treatment_year = temp.df$treatment_year)

assign(paste0("df", trt_vector[i]),tdistances_temp)
}


mod <- feols(dist~trt_type | site+expgroup +as.character(treatment_year), data = dfN)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(dist~as.factor(trt_type) , random = list(site = ~1,expgroup=~1, treatment_year=~1), data = dfN%>%
#             mutate(trt_type = fct_relevel(trt_type, "control")))
#summary(mod)

x <- ggpredict(mod, "trt_type")
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymin = predicted-std.error, ymax = predicted+std.error))+
  ylab("Distance among sites")+
  xlab("")+
  theme_base()
#dfN%>%
#  group_by(trt_type)%>%
#  dplyr::summarize(mean = mean(dist), se = sd(dist)/sqrt(n()), conf = se*1.96)%>%
#  ggplot(aes(trt_type, mean))+
#  geom_pointrange(aes(ymax = mean+conf,ymin = mean-conf))+
#  xlab("")+
#  ylab("Distance among sites")+
#  theme_base()


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/global_N.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 3.5,
  height = 3.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)


mod <- feols(dist~trt_type | site+expgroup +as.character(treatment_year), data = dfP)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(dist~as.factor(trt_type) , random = list(site = ~1,expgroup=~1, treatment_year=~1), data = dfP%>%
#             mutate(trt_type = fct_relevel(trt_type, "control")))
#summary(mod)

x <- ggpredict(mod, "trt_type")
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymin = predicted-std.error, ymax = predicted+std.error))+
  ylab("Distance among sites")+
  xlab("")+
  theme_base()

#dfP%>%
#  group_by(trt_type)%>%
#  dplyr::summarize(mean = mean(dist), se = sd(dist)/sqrt(n()), conf = se*1.96)%>%
#  ggplot(aes(trt_type, mean))+
#  geom_pointrange(aes(ymax = mean+conf,ymin = mean-conf))+
#  xlab("")+
#  ylab("Distance among sites")+
#  theme_base()


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/global_P.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 3.5,
  height = 3.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)

mod <- feols(dist~trt_type | site+expgroup +as.character(treatment_year), data = dfmult_nutrient)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(dist~as.factor(trt_type) , random = list(site = ~1,expgroup=~1, treatment_year=~1), data = dfmult_nutrient%>%
#             mutate(trt_type = fct_relevel(trt_type, "control")))
#summary(mod)

x <- ggpredict(mod, "trt_type")
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymin = predicted-std.error, ymax = predicted+std.error))+
  ylab("Distance among sites")+
  xlab("")+
  theme_base()

#dfmult_nutrient%>%
#  group_by(trt_type)%>%
#  dplyr::summarize(mean = mean(dist), se = sd(dist)/sqrt(n()), conf = se*1.96)%>%
#  ggplot(aes(trt_type, mean))+
#  geom_pointrange(aes(ymax = mean+conf,ymin = mean-conf))+
#  xlab("")+
#  ylab("Distance among sites")+
#  theme_base()


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/global_mult.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 3.5,
  height = 3.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)

mod <- feols(dist~trt_type | site+expgroup +as.character(treatment_year), data = dfirr)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(dist~as.factor(trt_type) , random = list(site = ~1,expgroup=~1, treatment_year=~1), data = dfirr%>%
#             mutate(trt_type = fct_relevel(trt_type, "control")))
#summary(mod)

x <- ggpredict(mod, "trt_type")
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymin = predicted-std.error, ymax = predicted+std.error))+
  ylab("Distance among sites")+
  xlab("")+
  theme_base()

#dfirr%>%
#  group_by(trt_type)%>%
#  dplyr::summarize(mean = mean(dist), se = sd(dist)/sqrt(n()), conf = se*1.96)%>%
#  ggplot(aes(trt_type, mean))+
#  geom_pointrange(aes(ymax = mean+conf,ymin = mean-conf))+
#  xlab("")+
#  ylab("Distance among sites")+
#  theme_base()


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/global_irr.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 3.5,
  height = 3.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)


mod <- feols(dist~trt_type | site+expgroup +as.character(treatment_year), data = dfCO2%>%
               mutate(trt_type = fct_relevel(trt_type, "control")))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(dist~as.factor(trt_type) , random = list(site = ~1,expgroup=~1, treatment_year=~1), data = dfCO2%>%
#             mutate(trt_type = fct_relevel(trt_type, "control")))
#summary(mod)

x <- ggpredict(mod, "trt_type")
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymin = predicted-std.error, ymax = predicted+std.error))+
  ylab("Distance among sites")+
  xlab("")+
  theme_base()

#dfCO2%>%
#  group_by(trt_type)%>%
#  dplyr::summarize(mean = mean(dist), se = sd(dist)/sqrt(n()), conf = se*1.96)%>%
#  ggplot(aes(trt_type, mean))+
#  geom_pointrange(aes(ymax = mean+conf,ymin = mean-conf))+
#  xlab("")+
#  ylab("Distance among sites")+
#  theme_base()


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/global_co2.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 3.5,
  height = 3.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)

mod <- feols(dist~trt_type | site+expgroup +as.character(treatment_year), data = `dfN*irr`)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(dist~as.factor(trt_type) , random = list(site = ~1,expgroup=~1, treatment_year=~1), data = `dfN*irr`%>%
#             mutate(trt_type = fct_relevel(trt_type, "control")))
#summary(mod)

x <- ggpredict(mod, "trt_type")
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymin = predicted-std.error, ymax = predicted+std.error))+
  ylab("Distance among sites")+
  xlab("")+
  theme_base()


#`dfN*irr`%>%
#  group_by(trt_type)%>%
#  dplyr::summarize(mean = mean(dist), se = sd(dist)/sqrt(n()), conf = se*1.96)%>%
#  ggplot(aes(trt_type, mean))+
#  geom_pointrange(aes(ymax = mean+conf,ymin = mean-conf))+
#  xlab("")+
#  ylab("Distance among sites")+
#  theme_base()


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/global_Nirr.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 3.5,
  height = 3.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)




###All together now!
temp.gow <- gowdis(summarize.cwm[6:ncol(summarize.cwm)])
temp.beta <- betadisper(temp.gow, group = summarize.cwm$trt_type, type = "centroid")
tdistances_temp <- data.frame(site = separate(summarize.cwm,expgroup, into = c("site", "project", "community"), sep = "::")$site,
                              expgroup = summarize.cwm$expgroup, trt_type = summarize.cwm$trt_type, treatment = summarize.cwm$treatment,  dist = temp.beta$dist, plot_mani = summarize.cwm$plot_mani, treatment_year = summarize.cwm$treatment_year)



mod <- feols(dist~plot_mani | site + expgroup + as.character(treatment_year), data = tdistances_temp)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(dist~plot_mani , random = list(site = ~1,expgroup=~1, treatment_year=~1), data = tdistances_temp)
#summary(mod)
#mod2 <- lme(dist~plot_mani+poly(plot_mani,2) , random = list(site = ~1,expgroup=~1, treatment_year=~1), data = tdistances_temp)
#summary(mod2)

x <- ggpredict(mod, "plot_mani")
tdistances_temp%>%
  group_by(plot_mani)%>%
  dplyr::summarize(mean = mean(dist), se = sd(dist)/sqrt(n()), sd = sd(dist), conf = se*1.96)%>%
ggplot( aes(plot_mani, mean))+
  #geom_pointrange( aes(ymax=mean+conf, ymin=mean-conf))+
  geom_smooth(data=x, aes(x=x, y=predicted), se = FALSE, color = "black")+
  geom_smooth(data=tdistances_temp, aes(y = dist,group = site),method = "lm",size = 0.5, se = FALSE, color = "lightgrey")+
#  geom_smooth(data=x, aes(x=x, y=predicted-std.error), se = FALSE, linetype = "dashed")+
  #geom_smooth(method = "loess")+
  xlab("Number of manipulations")+
  ylab("Distance among sites")+
  theme_base()


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/global_plotmani.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 4.5,
  height = 4,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)



mod <- feols(dist~plot_mani*treatment_year | site + expgroup , data = tdistances_temp)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#x <- ggpredict(mod, c("plot_mani", "treatment_year"))
#tdistances_temp%>%
#  group_by(plot_mani)%>%
#  dplyr::summarize(mean = mean(dist), se = sd(dist)/sqrt(n()))%>%
#ggplot( aes(plot_mani, mean))+
#  facet_wrap(~treatment_year)+
#  geom_pointrange(alpha = 0.1)+
  #geom_smooth(data=x, aes(x=x, y=predicted), se = FALSE)+
#  geom_smooth(method = "loess")+
#  theme_base()
  

tdistances_temp$any.treatment <-revalue(tdistances_temp$trt_type, c(N = "treatment",P = "treatment",irr = "treatment",mult_nutrient = "treatment",`irr*CO2` = "treatment",`N*irr*CO2` = "treatment",`mult_nutrient*irr` = "treatment",`N*CO2` = "treatment", `N*irr` = "treatment", CO2 = "treatment"
))
mod <- feols(dist~any.treatment | site + expgroup + as.character(treatment_year), data = tdistances_temp)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))

x <- ggpredict(mod, "any.treatment")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  xlab("")+
  ylab("Distance among sites")+
  theme_base()

#tdistances_temp%>%
#  group_by(any.treatment)%>%
#  dplyr::summarize(mean = mean(dist), se = sd(dist)/sqrt(n()), conf = se*1.96)%>%
#ggplot(aes(any.treatment, mean))+
#  geom_pointrange(aes(ymax = mean+conf,ymin = mean-conf))+
#  xlab("")+
#  ylab("Distance among sites")+
#  theme_base()

ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/global_overall.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 3.5,
  height = 3.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)

