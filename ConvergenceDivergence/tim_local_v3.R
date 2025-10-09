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

##must make dfs of sites with certain treatments so you can much up controls with only the appropriate sites for each treatment type
sites.n <- subset(mean.dist.df,  trt_type == "N")%>%dplyr::select(expgroup)%>%unique()
sites.p <- subset(mean.dist.df,  trt_type == "P")%>%dplyr::select(expgroup)%>%unique()
sites.multnutrient <- subset(mean.dist.df,  trt_type == "mult_nutrient")%>%dplyr::select(expgroup)%>%unique()
sites.irr <- subset(mean.dist.df,  trt_type == "irr")%>%dplyr::select(expgroup)%>%unique()
sites.co2 <- subset(mean.dist.df,  trt_type == "CO2")%>%dplyr::select(expgroup)%>%unique()
sites.nirr <- subset(mean.dist.df,  trt_type == "N*irr")%>%dplyr::select(expgroup)%>%unique()

#models to test results -
#stats for overall (any treatment)
mean.dist.df$any.treatment <-revalue(mean.dist.df$trt_type, c(N = "treatment",P = "treatment",irr = "treatment",mult_nutrient = "treatment",`irr*CO2` = "treatment",`N*irr*CO2` = "treatment",`mult_nutrient*irr` = "treatment",`N*CO2` = "treatment", `N*irr` = "treatment", CO2 = "treatment"
))
mod <- feols(mean_dist~any.treatment | site + expgroup +as.character(treatment_year)  ,data = subset(mean.dist.df, treatment_year != 0))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))


x <- ggpredict(mod, "any.treatment")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  #geom_violin(data = subset(mean.dist.df, treatment_year != 0), aes(any.treatment,mean_dist))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error, shape = x))+
  scale_shape_manual(values = c(2,17))+
  
    xlab("")+
    ylab("Distance among replicates within sites")+
    theme_base()+
    theme(legend.position = "none")



ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_overall_comp.pdf",
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



#stats for nitrogen treatment
mod <- feols(mean_dist~trt_type | site + expgroup +as.character(treatment_year)  ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.n$expgroup), treatment_year != 0), trt_type == "N"|trt_type=="control"))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))


x <- ggpredict(mod, "trt_type")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()

n_comp_stats <- x
n_comp_stats$set <- "N"


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_N_comp.pdf",
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


#stats for phosphorus treatment
mod <- feols(mean_dist~trt_type | site + expgroup +as.character(treatment_year)  , data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.p$expgroup), treatment_year != 0), trt_type == "P"|trt_type=="control"))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))

x <- ggpredict(mod, "trt_type")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()

p_comp_stats <- x
p_comp_stats$set <- "P"

ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_P_comp.pdf",
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

#stats for multiple nutrient addition
mod <- feols(mean_dist~trt_type | site + expgroup +as.character(treatment_year),data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.multnutrient$expgroup), treatment_year != 0), trt_type == "mult_nutrient"|trt_type=="control"))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))

x <- ggpredict(mod, "trt_type")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()

mult_comp_stats <- x
mult_comp_stats$set <- "mult"


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_mult_comp.pdf",
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

#stats for irrigation
mod <- feols(mean_dist~trt_type | site + expgroup +as.character(treatment_year),data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.irr$expgroup), treatment_year != 0), trt_type == "irr"|trt_type=="control"))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))

x <- ggpredict(mod, "trt_type")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()

irr_comp_stats <- x
irr_comp_stats$set <- "irr"


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_irr_comp.pdf",
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


#stats for co2
mod <- feols(mean_dist~trt_type | site + expgroup +as.character(treatment_year),data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.co2$expgroup), treatment_year != 0), trt_type == "CO2"|trt_type=="control")%>%
               mutate(trt_type = fct_relevel(trt_type, "control"))
)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))


x <- ggpredict(mod, "trt_type")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()

co2_comp_stats <- x
co2_comp_stats$set <- "co2"



ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_co2_comp.pdf",
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


#stats for N*irr
mod <- feols(mean_dist~trt_type | site + expgroup +as.character(treatment_year),data = 
               subset(subset(subset(mean.dist.df,  expgroup%in%sites.nirr$expgroup), treatment_year != 0), trt_type == "N*irr"|trt_type=="control")%>%
               mutate(trt_type = fct_relevel(trt_type, "control"))
)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))


x <- ggpredict(mod, "trt_type")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()

nxirr_comp_stats <- x
nxirr_comp_stats$set <- "Nxirr"


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_nirr_comp.pdf",
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


##combined figure for local composition

comp_stats <- rbind(n_comp_stats, p_comp_stats)%>%
  rbind(mult_comp_stats)%>%
  rbind(nxirr_comp_stats)%>%
  rbind(irr_comp_stats)%>%
  rbind(co2_comp_stats)%>%
  as.matrix()%>%
  as.data.frame()%>%
  mutate(predicted = as.numeric(predicted), std.error = as.numeric(std.error))

comp_stats%>%
   dplyr::mutate(set2 = factor(set, levels = c("N", "P", "mult", "irr", "Nxirr", "co2")))%>%
  dplyr::mutate(x = fct_relevel(x, "control"))%>%
ggplot( aes(set2, predicted, color = x))+
  geom_pointrange(aes(ymin = predicted - std.error, ymax = predicted+std.error, shape = x), position = position_dodge(0.3))+
  scale_color_manual(values = c("black", "#f2c300","#df0000", "darkorange1",  "#0099f6","purple", "#00b844"))+
  scale_shape_manual(values = c(2,17,17,17,17,17,17))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()+
  theme(legend.position = "None")


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_comp_alltreats.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 5,
  height = 3.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)


comp_stats%>%
  mutate(treatment = ifelse(x == "control", "control", "treatment"))%>%
  dplyr::select(predicted, std.error, set, treatment)%>%
  pivot_wider(names_from = "treatment", values_from = c("predicted", "std.error"))%>%
  mutate(treatment_minus_control = predicted_treatment-predicted_control)%>%
  dplyr::mutate(set2 = factor(set, levels = c("N", "P", "mult", "irr", "Nxirr", "co2")))%>%
ggplot(aes(x = set2, y = treatment_minus_control, color = set2))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_pointrange(aes(ymin = treatment_minus_control - std.error_treatment, ymax = treatment_minus_control+std.error_treatment), position = position_dodge(0.3))+
  scale_color_manual(values = c("#f2c300","#df0000", "darkorange1",  "#0099f6","purple", "#00b844"))+
  xlab("")+
  ylab("Distance among replicates within sites (treatment-control)")+
  theme_base()+
  theme(legend.position = "None")



ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_comp_alltreats_trtminuscon.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 5,
  height = 3.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)






##Stats about change over time: Nitrogen
mod <- feols(mean_dist~trt_type*treatment_year | site + expgroup ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.n$expgroup), treatment_year != 0), trt_type == "N"|trt_type=="control"))
summary(mod)


##Stats about change over time: Phosphorus
mod <- feols(mean_dist~trt_type*treatment_year | site + expgroup ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.p$expgroup), treatment_year != 0), trt_type == "P"|trt_type=="control"))
summary(mod)


subset(subset(subset(mean.dist.df,  expgroup%in%sites.p$expgroup), treatment_year != 0), trt_type == "P"|trt_type=="control")%>%
  ggplot(aes(treatment_year, mean_dist, color = trt_type))+
  geom_point()+
  geom_smooth()+
  theme_base()


##Stats about change over time: multiple nutrient addition
mod <- feols(mean_dist~trt_type*treatment_year | site + expgroup ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.multnutrient$expgroup), treatment_year != 0), trt_type == "mult_nutrient"|trt_type=="control"))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

subset(subset(subset(mean.dist.df,  expgroup%in%sites.multnutrient$expgroup), treatment_year != 0), trt_type == "mult_nutrient"|trt_type=="control")%>%
  ggplot(aes(treatment_year, mean_dist, color = trt_type))+
  geom_point()+
  geom_smooth()+
  theme_base()


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

sites.n <- subset(mean.dist.df,  trt_type == "N")%>%dplyr::select(expgroup)%>%unique()
sites.p <- subset(mean.dist.df,  trt_type == "P")%>%dplyr::select(expgroup)%>%unique()
sites.multnutrient <- subset(mean.dist.df,  trt_type == "mult_nutrient")%>%dplyr::select(expgroup)%>%unique()
sites.irr <- subset(mean.dist.df,  trt_type == "irr")%>%dplyr::select(expgroup)%>%unique()
sites.co2 <- subset(mean.dist.df,  trt_type == "CO2")%>%dplyr::select(expgroup)%>%unique()
sites.nirr <- subset(mean.dist.df,  trt_type == "N*irr")%>%dplyr::select(expgroup)%>%unique()


#stats for overall (any treatment)
mean.dist.df$any.treatment <-revalue(mean.dist.df$trt_type, c(N = "treatment",P = "treatment",irr = "treatment",mult_nutrient = "treatment",`irr*CO2` = "treatment",`N*irr*CO2` = "treatment",`mult_nutrient*irr` = "treatment",`N*CO2` = "treatment", `N*irr` = "treatment", CO2 = "treatment"
))
mod <- feols(mean_dist~any.treatment | site + expgroup +as.character(treatment_year) ,data = subset(mean.dist.df, treatment_year != 0))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))


x <- ggpredict(mod, "any.treatment")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error, shape = x))+
  scale_shape_manual(values = c(2,17))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()+
  theme(legend.position = "none")


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_overall_trait.pdf",
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

#stats for Nitrogen effect (traits)
mod <- feols(mean_dist~trt_type | site + expgroup+as.character(treatment_year) ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.n$expgroup), treatment_year != 0), trt_type == "N"|trt_type=="control"))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))

x <- ggpredict(mod, "trt_type")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()

n_trait_stats <- x
n_trait_stats$set <- "N"


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_N_trait.pdf",
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

#stats for phosphorus effect (traits)
mod <- feols(mean_dist~trt_type | site+ expgroup+as.character(treatment_year) ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.p$expgroup), treatment_year != 0), trt_type == "P"|trt_type=="control"))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))

x <- ggpredict(mod, "trt_type")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()

p_trait_stats <- x
p_trait_stats$set <- "P"


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_P_trait.pdf",
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

#stats for multiple nutrient addition effect (traits)
mod <- feols(mean_dist~trt_type | site + expgroup+as.character(treatment_year) ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.multnutrient$expgroup), treatment_year != 0), trt_type == "mult_nutrient"|trt_type=="control"))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))


x <- ggpredict(mod, "trt_type")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()

mult_trait_stats <- x
mult_trait_stats$set <- "mult"


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_mult_trait.pdf",
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

#stats for irrigation
mod <- feols(mean_dist~trt_type | site + expgroup +as.character(treatment_year),data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.irr$expgroup), treatment_year != 0), trt_type == "irr"|trt_type=="control"))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))


x <- ggpredict(mod, "trt_type")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()

irr_trait_stats <- x
irr_trait_stats$set <- "irr"

ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_irr_trait.pdf",
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

#stats for co2
mod <- feols(mean_dist~trt_type | site + expgroup +as.character(treatment_year),data = 
               subset(subset(subset(mean.dist.df,  expgroup%in%sites.co2$expgroup), treatment_year != 0), trt_type == "CO2"|trt_type=="control")%>%
               mutate(trt_type = fct_relevel(trt_type, "control"))
)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))

x <- ggpredict(mod, "trt_type")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()

co2_trait_stats <- x
co2_trait_stats$set <- "co2"

ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_co2_trait.pdf",
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


#stats for N*irr
mod <- feols(mean_dist~trt_type | site + expgroup +as.character(treatment_year),data = 
               subset(subset(subset(mean.dist.df,  expgroup%in%sites.nirr$expgroup), treatment_year != 0), trt_type == "N*irr"|trt_type=="control")%>%
               mutate(trt_type = fct_relevel(trt_type, "control"))
)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
control_se <- sd(fixef(mod)$expgroup)/sqrt(length(fixef(mod)$expgroup))

x <- ggpredict(mod, "trt_type")
x$std.error <- ifelse(x$x == "control", control_se[1], x$std.error)
ggplot(x, aes(x, predicted))+
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()

nxirr_trait_stats <- x
nxirr_trait_stats$set <- "nxirr"

ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_Nirr_trait.pdf",
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



#figure of all local traits
trait_stats <- rbind(n_trait_stats, p_trait_stats)%>%
  rbind(mult_trait_stats)%>%
  rbind(nxirr_trait_stats)%>%
  rbind(irr_trait_stats)%>%
  rbind(co2_trait_stats)%>%
  as.matrix()%>%
  as.data.frame()%>%
  mutate(predicted = as.numeric(predicted), std.error = as.numeric(std.error))

trait_stats%>%
  dplyr::mutate(set2 = factor(set, levels = c("N", "P", "mult", "irr", "nxirr", "co2")))%>%
  dplyr::mutate(x = fct_relevel(x, "control"))%>%
  ggplot( aes(set2, predicted, color = x))+
  geom_pointrange(aes(ymin = predicted - std.error, ymax = predicted+std.error, shape = x), position = position_dodge(0.3))+
  scale_color_manual(values = c("black","gold",  "blue", "pink", "purple", "green", "orange"))+
  scale_shape_manual(values = c(2,17,17,17,17,17,17))+
  xlab("")+
  ylab("Distance among replicates within sites")+
  theme_base()+
  theme(legend.position = "None")



ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_trait_alltreats.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 5,
  height = 3.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)


trait_stats%>%
  mutate(treatment = ifelse(x == "control", "control", "treatment"))%>%
  dplyr::select(predicted, std.error, set, treatment)%>%
  pivot_wider(names_from = "treatment", values_from = c("predicted", "std.error"))%>%
  mutate(treatment_minus_control = predicted_treatment-predicted_control)%>%
  dplyr::mutate(set2 = factor(set, levels = c("N", "P", "mult", "irr", "nxirr", "co2")))%>%
  ggplot(aes(x = set2, y = treatment_minus_control, color = set2))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_pointrange(aes(ymin = treatment_minus_control - std.error_treatment, ymax = treatment_minus_control+std.error_treatment), position = position_dodge(0.3))+
  scale_color_manual(values = c("#f2c300","#df0000", "darkorange1",  "#0099f6","purple", "#00b844"))+
  xlab("")+
  ylab("Distance among replicates within sites (treatment-control)")+
  theme_base()+
  theme(legend.position = "None")

ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_trait_alltreats_trtminuscon.pdf",
  plot = last_plot(),
  device = "pdf",
  path = NULL,
  scale = 1,
  width = 5,
  height = 3.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)



#stats for Nitrogen effect over time (traits)
mod <- feols(mean_dist~trt_type*treatment_year | site + expgroup ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.n$expgroup), treatment_year != 0), trt_type == "N"|trt_type=="control"))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))


#stats for phosphorus effect (traits)
mod <- feols(mean_dist~trt_type*treatment_year | site + expgroup ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.p$expgroup), treatment_year != 0), trt_type == "P"|trt_type=="control"))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

#stats for multiple nutrient effect over time (traits)
mod <- feols(mean_dist~trt_type*treatment_year | site + expgroup ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.multnutrient$expgroup), treatment_year != 0), trt_type == "mult_nutrient"|trt_type=="control"))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

#####################
###Compare species and trait responses
mean.dist.both <- left_join(mean.dist.comp, mean.dist.trait, by = c("site", "project", "community","expgroup","trt_type", "treatment", "plot_mani","treatment_year"))%>%
  mutate(mean_dist.comp = mean_dist.x,mean_dist.trait = mean_dist.y)


############
##BRING IN COVARIATES AND SEE IF THEY EXPLAIN BETA DIVERSITY
CoRRE_siteLocationClimate_Dec2021 <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/environmental data/CoRRE_siteLocationClimate_Dec2021.csv")

CoRRE_project_summary <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/CoRRE_project_summary.csv")
CoRRE_project_summary$project <- CoRRE_project_summary$project_name
CoRRE_project_summary$community <- CoRRE_project_summary$community_type
CoRRE_project_summary <- CoRRE_project_summary %>% dplyr::select(-c(project_name, community_type))



dist.both <- mean.dist.both%>%
  mutate(site_code = site)%>%
  left_join(CoRRE_siteLocationClimate_Dec2021, by = "site_code")

sites.n <- subset(dist.both,  trt_type == "N")%>%dplyr::select(expgroup)%>%unique()
sites.p <- subset(dist.both,  trt_type == "P")%>%dplyr::select(expgroup)%>%unique()
sites.multnutrient <- subset(dist.both,  trt_type == "mult_nutrient")%>%dplyr::select(expgroup)%>%unique()


n.df <- dist.both%>%
  subset(expgroup %in% sites.n$expgroup)%>%
  subset(trt_type=="N"|trt_type=="control")


p.df <- dist.both%>%
  subset(expgroup %in% sites.p$expgroup)%>%
  subset(trt_type=="P"|trt_type=="control")

mult.df <- dist.both%>%
  subset(expgroup %in% sites.multnutrient$expgroup)%>%
  subset(trt_type=="mult_nutrient"|trt_type=="control")


#NITROGEN

##NITROGEN
mod <- feols(mean_dist.comp~trt_type+trt_type*MAP|site + expgroup+as.character(treatment_year), data = n.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

mod <- feols(mean_dist.trait~trt_type+trt_type*MAP|site + expgroup+as.character(treatment_year), data = n.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

mod <- feols(mean_dist.comp~trt_type+trt_type*MAT| site+expgroup+as.character(treatment_year), data = n.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

mod <- feols(mean_dist.trait~trt_type+trt_type*MAT|site+expgroup+as.character(treatment_year), data = n.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))


#PHOSPHORUS
mod <- feols(mean_dist.comp~trt_type+trt_type*MAP|site+expgroup+as.character(treatment_year), data = p.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

mod <- feols(mean_dist.trait~trt_type+trt_type*MAP|site+expgroup+as.character(treatment_year), data = p.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

mod <- feols(mean_dist.comp~trt_type+trt_type*MAT|site+expgroup+as.character(treatment_year), data = p.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

mod <- feols(mean_dist.trait~trt_type+trt_type*MAT|site+expgroup+as.character(treatment_year), data = p.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

#mult nutrient
mod <- feols(mean_dist.comp~trt_type+trt_type*MAP|site+expgroup+as.character(treatment_year), data = mult.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

mod <- feols(mean_dist.trait~trt_type+trt_type*MAP|site+expgroup+as.character(treatment_year), data = mult.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

mod <- feols(mean_dist.comp~trt_type+trt_type*MAT|site+expgroup+as.character(treatment_year), data = mult.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

mod <- feols(mean_dist.trait~trt_type+trt_type*MAT|site+expgroup+as.character(treatment_year), data = mult.df)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
####could do the same as above for irr, CO2, N*irr ^^^


##with treatment information
treatment_info <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/basic dataset info/ExperimentInfo.csv")
treatment_info$trt_type <- revalue(treatment_info$trt_type, c("N*P" = "mult_nutrient"))
treatment_info$project <- treatment_info$project_name
treatment_info$community <- treatment_info$community_type
treatment_info <- treatment_info[,c("site_code", "project", "community", "plot_mani", "trt_type", "treatment", "n", "p",  "CO2", "precip", "temp")]
treatment_info <- unique(treatment_info)%>%
                  unite( expgroup, c("site_code", "project", "community"), sep = "::")%>%
                  group_by(expgroup, plot_mani, trt_type, treatment)%>%
                dplyr::summarize(n = mean(n), p = mean(p), CO2 = mean(CO2), precip = mean(precip), temp = mean(temp))


sites.n <- subset(mean.dist.both,  trt_type == "N")%>%dplyr::select(expgroup)%>%unique()
sites.p <- subset(mean.dist.both,  trt_type == "P")%>%dplyr::select(expgroup)%>%unique()
sites.multnutrient <- subset(mean.dist.both,  trt_type == "mult_nutrient")%>%dplyr::select(expgroup)%>%unique()
sites.irr <- subset(mean.dist.both,  trt_type == "irr")%>%dplyr::select(expgroup)%>%unique()
sites.co2 <- subset(mean.dist.both,  trt_type == "CO2")%>%dplyr::select(expgroup)%>%unique()
sites.nirr <- subset(mean.dist.both,  trt_type == "N*irr")%>%dplyr::select(expgroup)%>%unique()


n.df <- mean.dist.both%>%
  subset(expgroup %in% sites.n$expgroup)%>%
  subset(trt_type=="N"|trt_type=="control")%>%
  left_join( treatment_info, by = c("expgroup", "plot_mani", "trt_type", "treatment"))


p.df <- mean.dist.both%>%
  subset(expgroup %in% sites.p$expgroup)%>%
  subset(trt_type=="P"|trt_type=="control")%>%
  left_join( treatment_info, by = c("expgroup", "plot_mani", "trt_type", "treatment"))

mult.df <- mean.dist.both%>%
  subset(expgroup %in% sites.multnutrient$expgroup)%>%
  subset(trt_type=="mult_nutrient"|trt_type=="control")%>%
  left_join( treatment_info, by = c("expgroup", "plot_mani", "trt_type", "treatment"))

irr.df <- mean.dist.both%>%
  subset(expgroup %in% sites.irr$expgroup)%>%
  subset(trt_type=="irr"|trt_type=="control")%>%
  left_join( treatment_info, by = c("expgroup", "plot_mani", "trt_type", "treatment"))

co2.df <- mean.dist.both%>%
  subset(expgroup %in% sites.co2$expgroup)%>%
  subset(trt_type=="CO2"|trt_type=="control")%>%
  left_join( treatment_info, by = c("expgroup", "plot_mani", "trt_type", "treatment"))

nirr.df <- mean.dist.both%>%
  subset(expgroup %in% sites.nirr$expgroup)%>%
  subset(trt_type=="N*irr"|trt_type=="control")%>%
  left_join( treatment_info, by = c("expgroup", "plot_mani", "trt_type", "treatment"))


##Nitrogen gradient

mod <- feols(mean_dist.comp~ n | site+expgroup+as.character(treatment_year), data = subset(n.df, n != 0))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(mean_dist.comp~ n, random = list(site=~1,site=~1,expgroup=~1,treatment_year=~1), data = subset(n.df, n != 0))
#summary(mod)

mod <- feols(mean_dist.trait~ n | site+expgroup+as.character(treatment_year), data = subset(n.df, n != 0))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(mean_dist.trait~ n, random = list(site=~1,site=~1,expgroup=~1,treatment_year=~1), data = subset(n.df, n != 0))
#summary(mod)

ggplot(n.df, aes(x=n, y=mean_dist.comp, color=trt_type))+
  #facet_wrap(~treatment_year)+
  geom_point()+
  scale_color_manual(values = c("black", "#0099f6"))+
  geom_hline(yintercept = 0)+
  theme_base()

ggplot(n.df, aes(x=n, y=mean_dist.trait))+
  #facet_wrap(~treatment_year)+
  geom_point()+
  geom_smooth(method="lm", se = FALSE)+
  scale_color_manual(values = c("black", "#0099f6"))+
  geom_hline(yintercept = 0)+
  theme_base()


##P gradient

mod <- feols(mean_dist.comp~ p | site+expgroup+as.character(treatment_year), data = subset(p.df, p != 0))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(mean_dist.comp~ p, random = list(site=~1,site=~1,expgroup=~1,treatment_year=~1), data = subset(p.df, p != 0))
#summary(mod)

mod <- feols(mean_dist.trait~ p | site+expgroup+as.character(treatment_year), data = subset(p.df, p != 0))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(mean_dist.trait~ p, random = list(site=~1,site=~1,expgroup=~1,treatment_year=~1), data = subset(p.df, p != 0))
#summary(mod)

ggplot(p.df, aes(x=p, y=mean_dist.comp, color = trt_type))+
  #facet_wrap(~treatment_year)+
  geom_point()+
  geom_smooth(method="lm", se = FALSE)+
  scale_color_manual(values = c("black", "#00b844"))+
  geom_hline(yintercept = 0)+
  theme_base()

ggplot(p.df, aes(x=p, y=mean_dist.trait, color = trt_type))+
  #facet_wrap(~treatment_year)+
  geom_point()+
  scale_color_manual(values = c("black", "#00b844"))+
  geom_hline(yintercept = 0)+
  theme_base()


##irr 
mod <- feols(mean_dist.comp~ precip | site+expgroup+as.character(treatment_year), data = subset(irr.df, plot_mani != 0))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(mean_dist.comp~ precip, random = list(site=~1,site=~1,expgroup=~1,treatment_year=~1), data = subset(irr.df, plot_mani != 0))
#summary(mod)

mod <- feols(mean_dist.trait~ precip | site+expgroup+as.character(treatment_year), data = subset(irr.df, plot_mani != 0))
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))

ggplot(irr.df, aes(x=precip, y=mean_dist.comp, color=trt_type))+
  facet_wrap(~treatment_year)+
  geom_point()+
  scale_color_manual(values = c("black", "#6305dc"))+
  geom_hline(yintercept = 0)+
  theme_base()

ggplot(irr.df, aes(x=precip, y=mean_dist.trait, color=trt_type))+
  facet_wrap(~treatment_year)+
  geom_point()+
  scale_color_manual(values = c("black", "#6305dc"))+
  geom_hline(yintercept = 0)+
  theme_base()


###BIG FIGURE PUT IN ALL TOGETHER BABY PLOT MANI FOR THE WIN
#comp stats
mod <- feols(mean_dist.comp~ plot_mani | site+expgroup+as.character(treatment_year), data = mean.dist.both)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(mean_dist.comp~ plot_mani, random = list(site=~1,expgroup=~1,treatment_year=~1), data = mean.dist.both)
#summary(mod)
#mod2 <- lme(mean_dist.comp~ plot_mani+poly(plot_mani,2), random = list(site=~1,expgroup=~1,treatment_year=~1), data = mean.dist.both)
#summary(mod2)

#local comp figure
x <- ggpredict(mod, "plot_mani")
mean.dist.both%>%
  group_by(plot_mani, site)%>%
  dplyr::summarize(mean = mean(mean_dist.comp))%>%
  ggplot( aes(plot_mani, mean))+
  geom_smooth(aes(group = site),method = "lm",size = 0.5, se = FALSE, color = "lightgrey")+
  #geom_pointrange( aes(ymax=mean+conf, ymin=mean-conf))+
  geom_smooth(data=x, aes(x=x, y=predicted), se = FALSE, color = "black", size = 1.5)+
  xlab("Number of manipulations")+
  ylab("Distance among replicates within sites")+
  #ylim(0,0.17)+
  theme_base()

ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_plotmani_comp.pdf",
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


#trait stats
mod <- feols(mean_dist.trait~ plot_mani |site+expgroup+as.character(treatment_year), data = mean.dist.both)
summary(mod)
coeftest(mod, vcov = kernHAC(mod, kernel = "Bartlett"))
coeftest(mod, vcov = NeweyWest(mod, lag = 4))
#mod <- lme(mean_dist.trait~ plot_mani, random = list(site =~1,expgroup=~1, treatment_year=~1), data = mean.dist.both)#lme says same thing as feols
#summary(mod)
#mod2 <- lme(mean_dist.trait~ plot_mani+poly(plot_mani,2), random = list(site=~1,expgroup=~1, treatment_year=~1), data = mean.dist.both)#lme says same thing as feols
#summary(mod2)


#local trait figure
x <- ggpredict(mod, "plot_mani")
#mean.dist.both%>%
#  group_by(plot_mani)%>%
#  dplyr::summarize(mean = mean(mean_dist.trait), se = sd(mean_dist.trait)/sqrt(n()), sd = sd(mean_dist.trait), conf = se*1.96)%>%
#mean.dist.both%>%
#  group_by(plot_mani, site)%>%
#  dplyr::summarize(mean = mean(mean_dist.trait))%>%
#  group_by(plot_mani)%>%
#  dplyr::summarize(mean = mean(mean), se = sd(mean)/sqrt(n()), sd = sd(mean), conf = se*1.96)%>%
mean.dist.both%>%
  group_by(plot_mani, site)%>%
  dplyr::summarize(mean = mean(mean_dist.trait))%>%
  ggplot( aes(plot_mani, mean))+
  geom_smooth(aes(group = site),method = "lm",size = 0.5, se = FALSE, color = "lightgrey")+
  #geom_pointrange( aes(ymax=mean+conf, ymin=mean-conf))+
  geom_smooth(data=x, aes(x=x, y=predicted), se = FALSE, color = "black", size = 1.5)+
  xlab("Number of manipulations")+
  ylab("Distance among replicates within sites")+
  #ylim(0,0.17)+
  theme_base()

ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_plotmani_trait.pdf",
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


#subset(mean.dist.both,treatment_year == 2)
ggplot(mean.dist.both, aes(x=plot_mani, y=mean_dist.trait))+
  #facet_wrap(~treatment_year)+
  geom_point(, alpha = 0.1)+
  geom_smooth(data = x, aes(x, predicted), se = FALSE)+
  #geom_smooth(aes( color = site), se = FALSE)+
  #  geom_smooth(method = "loess")+
  # scale_color_manual(values = c("black", "#00b844"))+
  geom_hline(yintercept = 0)+
  theme_base()

ggplot(mean.dist.both, aes(x=plot_mani, y=mean_dist.trait))+
  #facet_wrap(~site)+
  geom_point(aes(color = site), alpha = 0.1)+
  geom_smooth(data = x, aes(x, predicted), se = FALSE)+
  #geom_smooth(aes( color = site), se = FALSE)+
  geom_smooth(aes(color = site),method = "lm", se = FALSE)+
  # scale_color_manual(values = c("black", "#00b844"))+
  geom_hline(yintercept = 0)+
  theme_base()



#below tries some stuff over time but no big deal
mean.dist.both%>%
  group_by(plot_mani, treatment_year)%>%
  dplyr::summarize(mean = mean(mean_dist.comp), se = sd(mean_dist.comp)/sqrt(n()), sd = sd(mean_dist.comp), conf = se*1.96)%>%
  ggplot( aes(treatment_year, mean, color = factor(plot_mani)))+
  facet_wrap(~plot_mani)+
  geom_pointrange( aes(ymax=mean+conf, ymin=mean-conf))+
  geom_smooth(method = "loess")+
  xlab("Treatment year")+
  ylab("Distance among replciates within sites")+
  #ylim(0,0.5)+
  theme_base()

mean.dist.both%>%
 # group_by(plot_mani, treatment_year)%>%
  #dplyr::summarize(mean = mean(mean_dist.comp), se = sd(mean_dist.comp)/sqrt(n()), sd = sd(mean_dist.comp), conf = se*1.96)%>%
  ggplot( aes(treatment_year, mean_dist.comp, color = factor(site)))+
  facet_wrap(~plot_mani)+
  
  geom_point()+
  #geom_pointrange( aes(ymax=mean+conf, ymin=mean-conf))+
  geom_smooth(method = "lm", se = FALSE)+
  xlab("Treatment year")+
  ylab("Distance among replciates within sites")+
  #ylim(0,0.5)+
  theme_base()

