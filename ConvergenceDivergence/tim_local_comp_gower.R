#only local composition in this file but ith gower


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
  temp.distances <- vegdist(temp.wide[10:ncol(temp.wide)], method = "gower")
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
  geom_pointrange(aes(ymax = predicted+std.error, ymin = predicted-std.error, shape = x), size = 1.5)+
  scale_shape_manual(values = c(2,17))+
  xlab("")+
  ylab("Taxonomic disperion within sites")+
  theme_base()+
  theme(legend.position = "none")



ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_overall_comp_gower.pdf",
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
  ylab("Taxonomic dispersion within sites")+
  theme_base()

n_comp_stats <- x
n_comp_stats$set <- "N"


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_N_comp_gower.pdf",
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
  ylab("Taxonomic dispersion within sites")+
  theme_base()

p_comp_stats <- x
p_comp_stats$set <- "P"

ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_P_comp_gower.pdf",
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
  ylab("Taxonomic dispersion within sites")+
  theme_base()

mult_comp_stats <- x
mult_comp_stats$set <- "mult"


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_mult_comp_gower.pdf",
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
  ylab("Taxonomic dispersion within sites")+
  theme_base()

irr_comp_stats <- x
irr_comp_stats$set <- "irr"


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_irr_comp_gower.pdf",
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
  ylab("Taxonomic dispersion within sites")+
  theme_base()

co2_comp_stats <- x
co2_comp_stats$set <- "co2"



ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_co2_comp_gower.pdf",
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
  ylab("Taxonomic dispersion within sites")+
  theme_base()

nxirr_comp_stats <- x
nxirr_comp_stats$set <- "Nxirr"


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_nirr_comp_gower.pdf",
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
  ylab("Taxonomic dispersion within sites (trt-ctrl)")+
  theme_base()+
  theme(legend.position = "None")


ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_comp_alltreats_gower.pdf",
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
  ylab("Taxonomic dispersion within sites (trt-ctrl)")+
  theme_base()+
  theme(legend.position = "None")



ggsave(
  "C:/Users/ohler/Dropbox/Tim Work/sCoRRE/Beta div/figures/local_comp_alltreats_trtminuscon_gower.pdf",
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
