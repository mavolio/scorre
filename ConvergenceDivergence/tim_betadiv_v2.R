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
#  subset(treatment_year <=5& treatment_year > 0) #only use treatment data and subset the number of years to be used

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

#Last year of experiment
#lastyear <- ddply(experimentinfo, .(site_code, project_name, community_type),
#                  function(x)data.frame(
#                    last_trt_yr = max(x$calendar_year)
#                  ))

#minimum treatment length
#nyear <- experimentinfo[c("site_code", "project_name", "community_type", "treatment_year")] %>%
#  unique()%>%
 # subset(treatment_year != 0)%>%
 # ddply(.(site_code, project_name, community_type),
 #       function(x)data.frame(
 #         n.trt.yrs = max(x$treatment_year)
 #       ))

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

test <- test%>%
  subset( trt_type == "control" | trt_type == "N" | trt_type == "P" | trt_type == "irr" |# trt_type == "drought"  | trt_type == "temp"| 
            trt_type == "mult_nutrient" #|trt_type == "mult_GCD"
          | trt_type == "CO2" | trt_type == "irr*CO2"  |trt_type == "N*irr*CO2" | trt_type == "mult_nutrient*irr" |trt_type == "N*CO2"|trt_type == "N*irr"
  )  #keep only the focal treatments

#Set minimum treatment years. Note that criteria is relaxed for drought experiments for: reasons
#N <-  test[test$trt_type %in% "N",]
#P <-  test[test$trt_type %in% "P",]
#irr <-  test[test$trt_type %in% "irr",]
#CO2 <-  subset(test[test$trt_type %in% "CO2",], n.trt.yrs >= 6)
#temp <-  test[test$trt_type %in% "temp",]
#mult_nutrient <-  test[test$trt_type %in% "mult_nutrient",]
#mult_GCD <-  subset(test[test$trt_type %in% "mult_GCD",], n.trt.yrs >= 6)
#drought <-  test[test$trt_type %in% "drought",]
#control <-  test[test$trt_type %in% "control",]

#test <- bind_rows(N, P, #irr, temp, 
#                  mult_nutrient, #drought, 
#                  control
                  #, mult_GCD, CO2
#)

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

sites <- unite(sites, temp, c("project_name", "community_type"), sep = "::", remove = FALSE)

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


trt.df <- subset(mean.dist.df, plot_mani >= 1)%>% #treatment data
  dplyr::rename(dist.trt = mean_dist)
con.df <- subset(mean.dist.df, plot_mani == 0)%>% #control data
  dplyr::rename(dist.con = mean_dist)%>%
  dplyr::select(expgroup, dist.con, treatment_year)

#lrr.df <- left_join(trt.df, con.df, by = c("expgroup","treatment_year"))%>%
#  mutate(lrr = log(dist.trt/dist.con))%>%
#  mutate(con_minus_trt = dist.con-dist.trt)


            
            
#lrr.df.conf <- lrr.df%>%
#  ddply(.(trt_type, treatment_year), function(x)data.frame(
#    lrr.mean = mean(x$lrr),
#    lrr.error = qt(0.975, df=length(x$trt_type)-1)*sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type)-1),
#    lrr.se = sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type)),
#    num_experiments = length(x$expgroup)
#  ))

#lrr.df.conf$trt_type <- factor(lrr.df.conf$trt_type, levels = c("drought", "irr", "temp", "N", "P", "mult_nutrient"#, "mult_GCD", "CO2"
#))
#lrr.df.conf$min <- lrr.df.conf$lrr.mean-lrr.df.conf$lrr.error
#lrr.df.conf$max <- lrr.df.conf$lrr.mean+lrr.df.conf$lrr.error 


#visualize
#lrr.df.conf <- lrr.df%>%
#  ddply(.(trt_type, treatment_year), function(x)data.frame(
#    lrr.mean = mean(x$lrr),
#    lrr.error = qt(0.975, df=length(x$trt_type)-1)*sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type)-1),
#    lrr.se = sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type)),
#    num_experiments = length(x$expgroup)
#  ))

#lrr.df.conf$trt_type <- factor(lrr.df.conf$trt_type, levels = c("drought", "irr", "temp", "N", "P", "mult_nutrient"#, "mult_GCD", "CO2"
#))
#lrr.df.conf$min <- lrr.df.conf$lrr.mean-lrr.df.conf$lrr.error
#lrr.df.conf$max <- lrr.df.conf$lrr.mean+lrr.df.conf$lrr.error



#ggplot(subset(subset(lrr.df.conf, treatment_year >=1), trt_type == "N"|trt_type == "P"|trt_type == "mult_nutrient"), aes(trt_type, lrr.mean, color = trt_type))+
#  facet_wrap(~treatment_year)+
#  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
#  geom_pointrange(aes(ymin = lrr.mean-lrr.error, ymax = lrr.mean+lrr.error), size = 1.5)+
#  xlab("")+
#  ylab("Species composition LRR distance between plots within treatment")+
#  scale_color_manual(values = c("#df0000","#0099f6", "orange", "#00b844","#f2c300","#6305dc", "black"))+
#  theme_base()+
#  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


#ggsave(
#  "C:/Users/ohler/Documents/converge-diverge/fig1_comp.pdf",
#  plot = last_plot(),
#  device = "pdf",
#  path = NULL,
#  scale = 1,
#  width = 7,
#  height = 6,
#  units = c("in"),
#  dpi = 600,
#  limitsize = TRUE
#)

##must make dfs of sites with certain treatments so you can much up controls with only the appropriate sites for each treatment type
sites.n <- subset(mean.dist.df,  trt_type == "N")%>%dplyr::select(expgroup)%>%unique()
sites.p <- subset(mean.dist.df,  trt_type == "P")%>%dplyr::select(expgroup)%>%unique()
sites.multnutrient <- subset(mean.dist.df,  trt_type == "mult_nutrient")%>%dplyr::select(expgroup)%>%unique()
sites.irr <- subset(mean.dist.df,  trt_type == "irr")%>%dplyr::select(expgroup)%>%unique()
sites.co2 <- subset(mean.dist.df,  trt_type == "CO2")%>%dplyr::select(expgroup)%>%unique()
sites.nirr <- subset(mean.dist.df,  trt_type == "N*irr")%>%dplyr::select(expgroup)%>%unique()

#models to test results -
#mod <- lmer(lrr~0+ trt_type + (1|expgroup)+ (1|treatment_year), data = subset(subset(lrr.df,  trt_type == "N"|trt_type =="mult_nutrient"|trt_type=="P"), treatment_year != 0))
#summary(mod)



#stats for nitrogen treatment
mod <- feols(mean_dist~trt_type | site + expgroup +treatment_year  ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.n$expgroup), treatment_year != 0), trt_type == "N"|trt_type=="control"))
summary(mod)

mean.dist.df%>%
  subset( expgroup%in%sites.n$expgroup)%>%
    subset(trt_type == "N"|trt_type=="control")%>%
  subset( treatment_year != 0)%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(trt_type, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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
mod <- feols(mean_dist~trt_type | site + expgroup +treatment_year  ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.p$expgroup), treatment_year != 0), trt_type == "P"|trt_type=="control"))
summary(mod)

mean.dist.df%>%
  subset( expgroup%in%sites.p$expgroup)%>%
  subset(trt_type == "P"|trt_type=="control")%>%
  subset( treatment_year != 0)%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(trt_type, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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
mod <- feols(mean_dist~trt_type | site + expgroup +treatment_year  ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.multnutrient$expgroup), treatment_year != 0), trt_type == "mult_nutrient"|trt_type=="control"))
summary(mod)

mean.dist.df%>%
  subset( expgroup%in%sites.multnutrient$expgroup)%>%
  subset(trt_type == "mult_nutrient"|trt_type=="control")%>%
  subset( treatment_year != 0)%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(trt_type, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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
mod <- feols(mean_dist~trt_type | site + expgroup +treatment_year  ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.irr$expgroup), treatment_year != 0), trt_type == "irr"|trt_type=="control"))
summary(mod)

mean.dist.df%>%
  subset( expgroup%in%sites.irr$expgroup)%>%
  subset(trt_type == "irr"|trt_type=="control")%>%
  subset( treatment_year != 0)%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(trt_type, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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
mod <- feols(mean_dist~trt_type | site + expgroup +treatment_year  ,data = 
               subset(subset(subset(mean.dist.df,  expgroup%in%sites.co2$expgroup), treatment_year != 0), trt_type == "CO2"|trt_type=="control")%>%
               mutate(trt_type = fct_relevel(trt_type, "control"))
             )
summary(mod)

mean.dist.df%>%
  subset( expgroup%in%sites.co2$expgroup)%>%
  subset(trt_type == "CO2"|trt_type=="control")%>%
  mutate(trt_type = fct_relevel(trt_type, "control"))%>%
  subset( treatment_year != 0)%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(trt_type, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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
mod <- feols(mean_dist~trt_type | site + expgroup +treatment_year  ,data = 
               subset(subset(subset(mean.dist.df,  expgroup%in%sites.nirr$expgroup), treatment_year != 0), trt_type == "N*irr"|trt_type=="control")%>%
               mutate(trt_type = fct_relevel(trt_type, "control"))
)
summary(mod)

mean.dist.df%>%
  subset( expgroup%in%sites.nirr$expgroup)%>%
  subset(trt_type == "N*irr"|trt_type=="control")%>%
  subset( treatment_year != 0)%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(trt_type, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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


##dataframes and plots for figures
#n.df <- subset(subset(subset(mean.dist.df,  expgroup%in%sites.n$expgroup), treatment_year != 0), trt_type == "N"|trt_type=="control")%>%
#        group_by(trt_type)%>%
#        dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()))

#p.df <- subset(subset(subset(mean.dist.df,  expgroup%in%sites.p$expgroup), treatment_year != 0), trt_type == "P"|trt_type=="control")%>%
#  group_by(trt_type)%>%
#  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()))

#mult.df <- subset(subset(subset(mean.dist.df,  expgroup%in%sites.multnutrient$expgroup), treatment_year != 0), trt_type == "mult_nutrient"|trt_type=="control")%>%
#  group_by(trt_type)%>%
#  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()))

#irr.df <- subset(subset(subset(mean.dist.df,  expgroup%in%sites.irr$expgroup), treatment_year != 0), trt_type == "irr"|trt_type=="control")%>%
#  group_by(trt_type)%>%
#  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()))

#co2.df <- subset(subset(subset(mean.dist.df,  expgroup%in%sites.co2$expgroup), treatment_year != 0), trt_type == "CO2"|trt_type=="control")%>%
#  group_by(trt_type)%>%
#  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()))

#nirr.df <- subset(subset(subset(mean.dist.df,  expgroup%in%sites.nirr$expgroup), treatment_year != 0), trt_type == "N*irr"|trt_type=="control")%>%
#  group_by(trt_type)%>%
#  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()))

  
#ggplot(n.df, aes(trt_type, mean, color = trt_type))+
#  geom_pointrange(aes(ymin = mean-se, ymax = mean+se, color = trt_type ))+
#  scale_color_manual(values = c("black", "#0099f6"))+
#  xlab("")+
#  ylab("Beta diversity")+
#  theme_base()+
#  theme(legend.position="none")

#"#df0000","#0099f6", "orange", "#00b844","#f2c300","#6305dc"

#ggplot(p.df, aes(trt_type, mean, color = trt_type))+
#  geom_pointrange(aes(ymin = mean-se, ymax = mean+se, color = trt_type ))+
#  scale_color_manual(values = c("black", "#00b844"))+
#  xlab("")+
#  ylab("Beta diversity")+
#  theme_base()+
#  theme(legend.position="none")

#ggplot(mult.df, aes(trt_type, mean, color = trt_type))+
#  geom_pointrange(aes(ymin = mean-se, ymax = mean+se, color = trt_type ))+
#  scale_color_manual(values = c("black", "#6305dc"))+
#  xlab("")+
#  ylab("Beta diversity")+
#  theme_base()+
#  theme(legend.position="none")

#ggplot(irr.df, aes(trt_type, mean, color = trt_type))+
#  geom_pointrange(aes(ymin = mean-se, ymax = mean+se, color = trt_type ))+
#  scale_color_manual(values = c("black", "#6305dc"))+
#  xlab("")+
#  ylab("Beta diversity")+
#  theme_base()+
#  theme(legend.position="none")

#ggplot(co2.df, aes(trt_type, mean, color = trt_type))+
#  geom_pointrange(aes(ymin = mean-se, ymax = mean+se, color = trt_type ))+
#  scale_color_manual(values = c("black", "#6305dc"))+
#  xlab("")+
#  ylab("Beta diversity")+
#  theme_base()+
#  theme(legend.position="none")

#ggplot(nirr.df, aes(trt_type, mean, color = trt_type))+
#  geom_pointrange(aes(ymin = mean-se, ymax = mean+se, color = trt_type ))+
#  scale_color_manual(values = c("black", "#6305dc"))+
#  xlab("")+
#  ylab("Beta diversity")+
#  theme_base()+
#  theme(legend.position="none")



#library(ggeffects)
#x <- ggpredict(mod, c("trt_type"))
#plot(x, ci=FALSE)#+
  #ylim(0, 200)+
#ggplot(x , aes(x, predicted))+
#  geom_pointrange(aes(ymax = conf.high, ymin = conf.low))+
#  geom_hline(yintercept = 0)+
#  theme_base()



##Stats about change over time: Nitrogen
mod <- feols(mean_dist~trt_type*treatment_year | site + expgroup ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.n$expgroup), treatment_year != 0), trt_type == "N"|trt_type=="control"))
summary(mod)


#stats for overall (any treatment)
mean.dist.df$any.treatment <-revalue(mean.dist.df$trt_type, c(N = "treatment",P = "treatment",irr = "treatment",mult_nutrient = "treatment",`irr*CO2` = "treatment",`N*irr*CO2` = "treatment",`mult_nutrient*irr` = "treatment",`N*CO2` = "treatment", `N*irr` = "treatment", CO2 = "treatment"
))
mod <- feols(mean_dist~any.treatment | site + expgroup +treatment_year  ,data = subset(mean.dist.df, treatment_year != 0))
summary(mod)

tidy_mod <- broom::tidy(mod)
coefplot(mod, keep = c("Intercept","control", "treatment"),intercept = TRUE)

mean.dist.df%>%
  subset( treatment_year != 0)%>%
  group_by(any.treatment)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(any.treatment, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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

##Stats about change over time: Phosphorus
mod <- feols(mean_dist~trt_type*treatment_year | site + expgroup ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.p$expgroup), treatment_year != 0), trt_type == "P"|trt_type=="control"))
summary(mod)


subset(subset(subset(mean.dist.df,  expgroup%in%sites.p$expgroup), treatment_year != 0), trt_type == "P"|trt_type=="control")%>%
  ggplot(aes(treatment_year, mean_dist, color = trt_type))+
  geom_point()+
  geom_smooth()+
  theme_base()

#x <- ggpredict(mod, c("trt_type", "treatment_year"))
#ggplot(x , aes(x = group, y= predicted, color=x))+
#  geom_pointrange(aes(ymax = conf.high, ymin = conf.low),position= position_dodge(width = 0.2))+
#  scale_color_manual(values = c("black", "#00b844"))+
#  xlab("Treatment year")+
#  ylab("Beta diversity")+
#  theme_base()

##Stats about change over time: multiple nutrient addition
mod <- feols(mean_dist~trt_type*treatment_year | site + expgroup ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.multnutrient$expgroup), treatment_year != 0), trt_type == "mult_nutrient"|trt_type=="control"))
summary(mod)


subset(subset(subset(mean.dist.df,  expgroup%in%sites.multnutrient$expgroup), treatment_year != 0), trt_type == "mult_nutrient"|trt_type=="control")%>%
  ggplot(aes(treatment_year, mean_dist, color = trt_type))+
  geom_point()+
  geom_smooth()+
  theme_base()

#x <- ggpredict(mod, c("trt_type", "treatment_year"))
#ggplot(x , aes(x = group, y= predicted, color=x))+
#  geom_pointrange(aes(ymax = conf.high, ymin = conf.low),position= position_dodge(width = 0.2))+
#  scale_color_manual(values = c("black", "#6305dc"))+
#  xlab("Treatment year")+
#  ylab("Beta diversity")+
#  theme_base()

#ggsave(
#  "C:/Users/ohler/Documents/converge-diverge/fig1_comp.pdf",
#  plot = last_plot(),
#  device = "pdf",
#  path = NULL,
#  scale = 1,
#  width = 7,
#  height = 6,
#  units = c("in"),
#  dpi = 600,
#  limitsize = TRUE
#)



#lrr.df_species <- lrr.df

###Summarize sites being used
sites <- test%>%
  dplyr::select(site_code, project_name, community_type, trt_type, treatment)%>%
  tidyr::unite("expgroup", c("site_code", "project_name", "community_type"))%>%
  unique()%>%
  subset(trt_type != "control")

n <- sites%>%
  
  ddply(.(trt_type), function(x)data.frame(n = length(x$expgroup)))

#trt.by.year <- sites%>%
  #  tidyr::unite("expgroup", c("site_code", "project_name", "community_type"))%>%
#  ddply(.( trt_type, treatment_year), function(x)data.frame(n = length(x$expgroup)))



######
###The same stuff with traits but they include categorical traits
#CoRRE_CWMtraits <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/paper 2_PD and FD responses/data/CoRRE_CWMtraits_12142022.csv") #for now I'll just use this for categorical traits
#CoRRE_CWMtraits_cat <- CoRRE_CWMtraits[, c(   "site_code", "project_name","community_type", "plot_id", "treatment_year", "CWM.growth_form", "CWM.photosynthetic_pathway", "CWM.lifespan", "CWM.clonal", "CWM.mycorrhizal_type", "CWM.n_fixation")]

#CoRRE_CWMtraits_cat <- tidyr::unite(CoRRE_CWMtraits_cat, "rep", c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = TRUE)

#CoRRE_CWMtraits_cat$is.graminoid <- ifelse(CoRRE_CWMtraits_cat$CWM.growth_form == "graminoid", 1, 0)
#CoRRE_CWMtraits_cat$is.C4 <- ifelse(CoRRE_CWMtraits_cat$CWM.photosynthetic_pathway == "C$", 1, 0)
#CoRRE_CWMtraits_cat$is.perennial <- ifelse(CoRRE_CWMtraits_cat$CWM.lifespan == "perennial", 1, 0)
#CoRRE_CWMtraits_cat$is.clonal <- ifelse(CoRRE_CWMtraits_cat$CWM.clonal == "yes", 1, 0)
#CoRRE_CWMtraits_cat$is.AM <- ifelse(CoRRE_CWMtraits_cat$CWM.mycorrhizal_type == "arbuscular", 1, 0)
#CoRRE_CWMtraits_cat$is.n_fixer <- ifelse(CoRRE_CWMtraits_cat$CWM.n_fixation == "yes", 1, 0)
#CoRRE_CWMtraits_cat <- CoRRE_CWMtraits_cat[,c("rep", "treatment_year", "is.graminoid", "is.C4", "is.perennial", "is.clonal", "is.AM", "is.n_fixer")]

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
    
    
  )#%>%
  #left_join(CoRRE_CWMtraits_cat, by = c("rep", "treatment_year"))


#summarize.traits.continuous <- traits[,c("species_matched", "SLA", "LDMC", "leaf_N","plant_height_vegetative", "seed_dry_mass", "SRL")]
#summarize.traits.continuous <- unique(summarize.traits.continuous)
#summarize.traits.categorical <- traits[,c("species_matched", "growth_form", "photosynthetic_pathway", "lifespan", "clonal", "mycorrhizal_type", "n_fixation")]
#summarize.traits.categorical <- subset(summarize.traits.categorical, photosynthetic_pathway == "C3" | photosynthetic_pathway == "C4" | photosynthetic_pathway == "CAM")

#summarize.traits <- left_join(summarize.traits.continuous, summarize.traits.categorical, by = "species_matched")                                       
# reassigning row names
#summarize.traits <- unique(summarize.traits)




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


trt.df <- subset(mean.dist.df, plot_mani >= 1)%>%
  dplyr::rename(dist.trt = mean_dist)
con.df <- subset(mean.dist.df, plot_mani == 0)%>%
  dplyr::rename(dist.con = mean_dist)%>%
  dplyr::select(expgroup, dist.con, treatment_year)

#lrr.df <- left_join(trt.df, con.df, by = c("expgroup", "treatment_year"))%>%
#  mutate(lrr = log(dist.trt/dist.con))%>%
#  mutate(con_minus_trt = dist.trt/dist.con)

#lrr.df.conf <- lrr.df%>%
#  ddply(.(trt_type, treatment_year), function(x)data.frame(
#    lrr.mean = mean(x$lrr),
#    lrr.error = qt(0.975, df=length(x$trt_type)-1)*sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type)-1),
 #   lrr.se = sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type)),
#    num_experiments = length(x$expgroup)
#  ))


#lrr.df.conf$min <- lrr.df.conf$lrr.mean-lrr.df.conf$lrr.error
#lrr.df.conf$max <- lrr.df.conf$lrr.mean+lrr.df.conf$lrr.error



#visualize
#lrr.df.conf$trt_type <- factor(lrr.df.conf$trt_type, levels = c(#"drought", "irr", "temp", 
#  "N", "P", "mult_nutrient" 
                                                                #,"mult_GCD", "CO2"
#))
#ggplot(subset(subset(lrr.df.conf, treatment_year >= 1), trt_type== "N"| trt_type=="P"|trt_type == "mult_nutrient"), aes(trt_type, lrr.mean, color = trt_type))+
#facet_wrap(~treatment_year)+
#  geom_hline(yintercept = 0, size = 1, linetype = "dashed")+
#  geom_pointrange(aes(ymin = lrr.mean-lrr.error, ymax = lrr.mean+lrr.error), size = 1.5)+
#  xlab("")+
#  ylab("Trait LRR distance between plots within treatment")+
#  scale_color_manual(values = c("#df0000","#0099f6", "orange", "#00b844","#f2c300","#6305dc"))+
#  theme_base()+
#  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



#lrr.df_traits <- lrr.df






sites.n <- subset(mean.dist.df,  trt_type == "N")%>%dplyr::select(expgroup)%>%unique()
sites.p <- subset(mean.dist.df,  trt_type == "P")%>%dplyr::select(expgroup)%>%unique()
sites.multnutrient <- subset(mean.dist.df,  trt_type == "mult_nutrient")%>%dplyr::select(expgroup)%>%unique()
sites.irr <- subset(mean.dist.df,  trt_type == "irr")%>%dplyr::select(expgroup)%>%unique()
sites.co2 <- subset(mean.dist.df,  trt_type == "CO2")%>%dplyr::select(expgroup)%>%unique()
sites.nirr <- subset(mean.dist.df,  trt_type == "N*irr")%>%dplyr::select(expgroup)%>%unique()

#models to test results
#mod <- lmer(lrr~0+trt_type+ (1|expgroup), data = subset(subset(lrr.df_traits, trt_type == "N" | trt_type == "P"| trt_type == "mult_nutrient"), treatment_year != 0))
#summary(mod)

#mod <- feols(lrr~0+trt_type | expgroup +treatment_year  ,data = subset(subset(lrr.df_traits,  trt_type == "N"|trt_type =="mult_nutrient"|trt_type=="P"), treatment_year != 0))
#summary(mod)

n.df <- subset(subset(subset(mean.dist.df,  expgroup%in%sites.n$expgroup), treatment_year != 0), trt_type == "N"|trt_type=="control")%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()))

p.df <-  subset(subset(subset(mean.dist.df,  expgroup%in%sites.p$expgroup), treatment_year != 0), trt_type == "P"|trt_type=="control")%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()))

mult.df <-  subset(subset(subset(mean.dist.df,  expgroup%in%sites.multnutrient$expgroup), treatment_year != 0), trt_type == "mult_nutrient"|trt_type=="control")%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()))


irr.df <- subset(subset(subset(mean.dist.df,  expgroup%in%sites.irr$expgroup), treatment_year != 0), trt_type == "irr"|trt_type=="control")%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()))

co2.df <- subset(subset(subset(mean.dist.df,  expgroup%in%sites.co2$expgroup), treatment_year != 0), trt_type == "CO2"|trt_type=="control")%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()))

nirr.df <- subset(subset(subset(mean.dist.df,  expgroup%in%sites.nirr$expgroup), treatment_year != 0), trt_type == "N*irr"|trt_type=="control")%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()))






#stats for overall (any treatment)
mean.dist.df$any.treatment <-revalue(mean.dist.df$trt_type, c(N = "treatment",P = "treatment",irr = "treatment",mult_nutrient = "treatment",`irr*CO2` = "treatment",`N*irr*CO2` = "treatment",`mult_nutrient*irr` = "treatment",`N*CO2` = "treatment", `N*irr` = "treatment", CO2 = "treatment"
))
mod <- feols(mean_dist~any.treatment | site + expgroup +treatment_year  ,data = subset(mean.dist.df, treatment_year != 0))
summary(mod)

mean.dist.df%>%
  subset( treatment_year != 0)%>%
  group_by(any.treatment)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(any.treatment, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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
mod <- feols(mean_dist~trt_type | site + expgroup+treatment_year ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.n$expgroup), treatment_year != 0), trt_type == "N"|trt_type=="control"))
summary(mod)

mean.dist.df%>%
  subset( expgroup%in%sites.n$expgroup)%>%
  subset(trt_type == "N"|trt_type=="control")%>%
  subset( treatment_year != 0)%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(trt_type, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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
mod <- feols(mean_dist~trt_type | site+ expgroup+treatment_year ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.p$expgroup), treatment_year != 0), trt_type == "P"|trt_type=="control"))
summary(mod)

mean.dist.df%>%
  subset( expgroup%in%sites.p$expgroup)%>%
  subset(trt_type == "P"|trt_type=="control")%>%
  subset( treatment_year != 0)%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(trt_type, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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
mod <- feols(mean_dist~trt_type | site + expgroup+treatment_year ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.multnutrient$expgroup), treatment_year != 0), trt_type == "mult_nutrient"|trt_type=="control"))
summary(mod)

mean.dist.df%>%
  subset( expgroup%in%sites.multnutrient$expgroup)%>%
  subset(trt_type == "mult_nutrient"|trt_type=="control")%>%
  subset( treatment_year != 0)%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(trt_type, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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
mod <- feols(mean_dist~trt_type | site + expgroup +treatment_year  ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.irr$expgroup), treatment_year != 0), trt_type == "irr"|trt_type=="control"))
summary(mod)

mean.dist.df%>%
  subset( expgroup%in%sites.irr$expgroup)%>%
  subset(trt_type == "irr"|trt_type=="control")%>%
  subset( treatment_year != 0)%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(trt_type, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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
mod <- feols(mean_dist~trt_type | site + expgroup +treatment_year  ,data = 
               subset(subset(subset(mean.dist.df,  expgroup%in%sites.co2$expgroup), treatment_year != 0), trt_type == "CO2"|trt_type=="control")%>%
               mutate(trt_type = fct_relevel(trt_type, "control"))
)
summary(mod)

mean.dist.df%>%
  subset( expgroup%in%sites.co2$expgroup)%>%
  subset(trt_type == "CO2"|trt_type=="control")%>%
  mutate(trt_type = fct_relevel(trt_type, "control"))%>%
  subset( treatment_year != 0)%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(trt_type, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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
mod <- feols(mean_dist~trt_type | site + expgroup +treatment_year  ,data = 
               subset(subset(subset(mean.dist.df,  expgroup%in%sites.nirr$expgroup), treatment_year != 0), trt_type == "N*irr"|trt_type=="control")%>%
               mutate(trt_type = fct_relevel(trt_type, "control"))
)
summary(mod)

mean.dist.df%>%
  subset( expgroup%in%sites.nirr$expgroup)%>%
  subset(trt_type == "N*irr"|trt_type=="control")%>%
  subset( treatment_year != 0)%>%
  group_by(trt_type)%>%
  dplyr::summarize(mean = mean(mean_dist), se = sd(mean_dist)/sqrt(n()), conf = se*1.96)%>%
  ggplot(aes(trt_type, mean))+
  geom_pointrange(aes(ymax = mean +conf, ymin = mean-conf))+
  xlab("")+
  ylab("Distance between replicates within sites")+
  theme_base()

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





#ggplot(n.df, aes(trt_type, mean, color = trt_type))+
#  geom_pointrange(aes(ymin = mean-se, ymax = mean+se, color = trt_type ))+
#  scale_color_manual(values = c("black", "#0099f6"))+
#  xlab("")+
#  ylab("Beta diversity")+
#  theme_base()+
#  theme(legend.position="none")

#ggplot(p.df, aes(trt_type, mean, color = trt_type))+
#  geom_pointrange(aes(ymin = mean-se, ymax = mean+se, color = trt_type ))+
#  scale_color_manual(values = c("black", "#00b844"))+
#  xlab("")+
#  ylab("Beta diversity")+
#  theme_base()+
#  theme(legend.position="none")

#ggplot(mult.df, aes(trt_type, mean, color = trt_type))+
#  geom_pointrange(aes(ymin = mean-se, ymax = mean+se, color = trt_type ))+
#  scale_color_manual(values = c("black", "#6305dc"))+
#  xlab("")+
#  ylab("Beta diversity")+
#  theme_base()+
#  theme(legend.position="none")

#ggplot(irr.df, aes(trt_type, mean, color = trt_type))+
#  geom_pointrange(aes(ymin = mean-se, ymax = mean+se, color = trt_type ))+
#  scale_color_manual(values = c("black", "#6305dc"))+
#  xlab("")+
#  ylab("Beta diversity")+
#  theme_base()+
#  theme(legend.position="none")

#ggplot(co2.df, aes(trt_type, mean, color = trt_type))+
#  geom_pointrange(aes(ymin = mean-se, ymax = mean+se, color = trt_type ))+
#  scale_color_manual(values = c("black", "#6305dc"))+
#  xlab("")+
#  ylab("Beta diversity")+
#  theme_base()+
#  theme(legend.position="none")

#ggplot(nirr.df, aes(trt_type, mean, color = trt_type))+
#  geom_pointrange(aes(ymin = mean-se, ymax = mean+se, color = trt_type ))+
#  scale_color_manual(values = c("black", "#6305dc"))+
#  xlab("")+
#  ylab("Beta diversity")+
#  theme_base()+
#  theme(legend.position="none")


#N
#mod <- lmer(lrr~treatment_year+ (1|expgroup), data = subset(lrr.df_traits, trt_type == "N" ))
#summary(mod)

#stats for Nitrogen effect over time (traits)
mod <- feols(mean_dist~trt_type*treatment_year | site + expgroup ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.n$expgroup), treatment_year != 0), trt_type == "N"|trt_type=="control"))
summary(mod)

#x <- ggpredict(mod, c("trt_type", "treatment_year"))
#ggplot(x , aes(x = group, y= predicted, color=x))+
#  geom_pointrange(aes(ymax = conf.high, ymin = conf.low),position= position_dodge(width = 0.2))+
#  scale_color_manual(values = c("black", "#0099f6"))+
#  xlab("Treatment year")+
#  ylab("Beta diversity")+
#  theme_base()



#stats for phosphorus effect (traits)
mod <- feols(mean_dist~trt_type*treatment_year | site + expgroup ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.p$expgroup), treatment_year != 0), trt_type == "P"|trt_type=="control"))
summary(mod)


#x <- ggpredict(mod, c("trt_type", "treatment_year"))
#ggplot(x , aes(x = group, y= predicted, color=x))+
#  geom_pointrange(aes(ymax = conf.high, ymin = conf.low),position= position_dodge(width = 0.2))+
#  scale_color_manual(values = c("black", "#00b844"))+
#  xlab("Treatment year")+
#  ylab("Beta diversity")+
#  theme_base()



#stats for multiple nutrient effect over time (traits)
mod <- feols(mean_dist~trt_type*treatment_year | site + expgroup ,data = subset(subset(subset(mean.dist.df,  expgroup%in%sites.multnutrient$expgroup), treatment_year != 0), trt_type == "mult_nutrient"|trt_type=="control"))
summary(mod)


#x <- ggpredict(mod, c("trt_type", "treatment_year"))
#ggplot(x , aes(x = group, y= predicted, color=x))+
#  geom_pointrange(aes(ymax = conf.high, ymin = conf.low),position= position_dodge(width = 0.2))+
#  scale_color_manual(values = c("black", "#6305dc"))+
#  xlab("Treatment year")+
#  ylab("Beta diversity")+
#  theme_base()

#ggplot(subset(lrr.df_traits, trt_type == "N" |trt_type == "P" |trt_type == "mult_nutrient" ), aes(treatment_year, lrr, color = trt_type))+
#  geom_point()+
#  geom_smooth(method = "lm")+
#  theme_base()



#ggsave(
#  "C:/Users/ohler/Documents/converge-diverge/fig1_trait.pdf",
#  plot = last_plot(),
#  device = "pdf",
#  path = NULL,
#  scale = 1,
#  width = 7,
#  height = 6,
#  units = c("in"),
#  dpi = 600,
#  limitsize = TRUE
#)



#####################
###Compare species and trait responses


mean.dist.both <- left_join(mean.dist.comp, mean.dist.trait, by = c("site", "project", "community","expgroup","trt_type", "treatment", "plot_mani","treatment_year"))%>%
                  mutate(mean_dist.comp = mean_dist.x,mean_dist.trait = mean_dist.y)

#some stats below but I'm not sure what they're good for
#mod <- feols(mean_dist.trait~mean_dist.comp+trt_type*treatment_year|expgroup , data = subset(subset(mean.dist.both, trt_type == "N"|trt_type == "control"), treatment_year != 0))
#summary(mod)

#mod <- feols(mean_dist.trait~mean_dist.comp+trt_type*treatment_year|expgroup , data = subset(subset(mean.dist.both, trt_type == "P"|trt_type == "control"), treatment_year != 0))
#summary(mod)

#mod <- feols(mean_dist.trait~mean_dist.comp+trt_type*treatment_year|expgroup , data = subset(subset(mean.dist.both, trt_type == "mult_nutrient"|trt_type == "control"), treatment_year != 0))
#summary(mod)

#ggplot(data = subset(subset(mean.dist.both, trt_type == "N"|trt_type == "control"), treatment_year != 0) , aes(x = mean_dist.comp, y= mean_dist.trait, color=trt_type))+
#  facet_wrap(~treatment_year)+
#  geom_point()+
#  geom_smooth(method = "lm")+
#  #  geom_hline(yintercept = 0)+
#  theme_base()


#visualize
#lrr_sp.tr$trt_type <- factor(lrr_sp.tr$trt_type, levels = c("drought", "irr", "temp", "N", "P", "mult_nutrient" 
                                                            #,"mult_GCD", "CO2"
#))

#library(ggpmisc)
#ggplot(subset(lrr_sp.tr, trt_type == "N"|trt_type == "P"|trt_type == "mult_nutrient"),  aes(lrr.species, lrr.traits))+
#  facet_wrap(~trt_type)+
#  geom_abline(slope = 1, linetype = "dotted")+
#  geom_point()+
#  geom_smooth(aes(color = trt_type),method = "lm", se = FALSE)+
#  geom_hline(yintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
#  geom_vline(xintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
#  ylim(-1.1, 2.4)+
#  xlim(-1.1, 2.4)+
  #stat_fit_glance(method = 'lm',
  #method.args = list(formula = formula),
  #                geom = 'text',
  #                aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")))+
#  xlab("LRR species composition beta diversity")+
#  ylab("LRR trait composition beta diversity")+
#  scale_color_manual(values = c("#df0000","#0099f6", "orange", "#00b844","#f2c300","#6305dc"))+
#  theme_base()+
#  theme(legend.position = "none")

#ggsave(
#  "C:/Users/ohler/Documents/converge-diverge/comp-trait_correlation.pdf",
#  plot = last_plot(),
#  device = "pdf",
#  path = NULL,
#  scale = 1,
#  width = 6,
#  height = 4.5,
#  units = c("in"),
#  dpi = 600,
#  limitsize = TRUE
#)




#models to test results

#mod <- lmer(lrr.traits~1+lrr.species*treatment_year + (1|site_code), data = subset(lrr_sp.tr, trt_type == "N" ))
#summary(mod)
#r.squaredGLMM(mod)


#ggplot(subset(lrr_sp.tr, trt_type == "N" ), aes(lrr.species, lrr.traits))+
#  facet_wrap(~treatment_year)+
#  geom_point()+
#  geom_smooth(method = "lm")+
#  geom_abline(slope = 1, linetype = "dotted")+
#  geom_hline(yintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
#  geom_vline(xintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
#  theme_base()



#mod <- lmer(lrr.traits~lrr.species*treatment_year + (1|site_code), data = subset(lrr_sp.tr, trt_type == "P"))
#summary(mod)
#r.squaredGLMM(mod)

#mod <- feols(lrr.traits~lrr.species | expgroup + treatment_year ,data = subset(lrr_sp.tr, trt_type == "P"))
#summary(mod)

#ggplot(subset(lrr_sp.tr, trt_type == "P"), aes(lrr.species, lrr.traits))+
#  facet_wrap(~treatment_year)+
#  geom_point()+
#  geom_smooth(method = "lm")+
#  geom_abline(slope = 1, linetype = "dotted")+
#  geom_hline(yintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
#  geom_vline(xintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
#  theme_base()


#mod <- lmer(lrr.traits~lrr.species*treatment_year + (1|site_code), data = subset(lrr_sp.tr, trt_type == "mult_nutrient"))
#summary(mod)
#r.squaredGLMM(mod)

#mod <- feols(lrr.traits~lrr.species | expgroup + treatment_year ,data = subset(lrr_sp.tr, trt_type == "mult_nutrient"))
#summary(mod)


#ggplot(subset(lrr_sp.tr, trt_type == "mult_nutrient" & treatment_year <= 10), aes(lrr.species, lrr.traits))+
#  facet_wrap(~treatment_year)+
#  geom_point()+
#  geom_smooth(method = "lm")+
#  geom_abline(slope = 1, linetype = "dotted")+
#  geom_hline(yintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
#  geom_vline(xintercept = 0, size = 1, linetype = "dashed", alpha = 0.5)+
#  theme_base()


############
##BRING IN COVARIATES AND SEE IF THEY EXPLAIN BETA DIVERSITY
CoRRE_siteLocationClimate_Dec2021 <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/environmental data/CoRRE_siteLocationClimate_Dec2021.csv")

CoRRE_project_summary <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/CoRRE_project_summary.csv")
CoRRE_project_summary$project <- CoRRE_project_summary$project_name
CoRRE_project_summary$community <- CoRRE_project_summary$community_type
CoRRE_project_summary <- CoRRE_project_summary %>% dplyr::select(-c(project_name, community_type))


#lrr.df_species <- tidyr::separate(lrr.df_species, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)
#lrr.df_traits <- tidyr::separate(lrr.df_traits, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)

#lrr_covariate <- left_join(lrr.df_species, CoRRE_project_summary, by = c("site_code", "project", "community"))
#lrr_covariate_traits <- left_join(lrr.df_traits, CoRRE_project_summary, by = c("site_code", "project", "community"))



#library(grf)

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

#tau.forest <- causal_forest(x = matrix(n.df$MAP, n.df$MAT, n.df$rrich), y = n.df$mean_dist.comp, w = n.df$trt_type)

##NITROGEN
mod <- feols(mean_dist.comp~trt_type+trt_type*MAP|site + expgroup+treatment_year, data = n.df)
summary(mod)

mod <- feols(mean_dist.trait~trt_type+trt_type*MAP|site + expgroup+treatment_year, data = n.df)
summary(mod)

mod <- feols(mean_dist.comp~trt_type+trt_type*MAT| site+expgroup+treatment_year, data = n.df)
summary(mod)

mod <- feols(mean_dist.trait~trt_type+trt_type*MAT|site+expgroup+treatment_year, data = n.df)
summary(mod)


#PHOSPHORUS
mod <- feols(mean_dist.comp~trt_type+trt_type*MAP|site+expgroup+treatment_year, data = p.df)
summary(mod)

mod <- feols(mean_dist.trait~trt_type+trt_type*MAP|site+expgroup+treatment_year, data = p.df)
summary(mod)

mod <- feols(mean_dist.comp~trt_type+trt_type*MAT|site+expgroup+treatment_year, data = p.df)
summary(mod)

mod <- feols(mean_dist.trait~trt_type+trt_type*MAT|site+expgroup+treatment_year, data = p.df)
summary(mod)

#mult nutrient
mod <- feols(mean_dist.comp~trt_type+trt_type*MAP|site+expgroup+treatment_year, data = mult.df)
summary(mod)

mod <- feols(mean_dist.trait~trt_type+trt_type*MAP|site+expgroup+treatment_year, data = mult.df)
summary(mod)

mod <- feols(mean_dist.comp~trt_type+trt_type*MAT|site+expgroup+treatment_year, data = mult.df)
summary(mod)

mod <- feols(mean_dist.trait~trt_type+trt_type*MAT|site+expgroup+treatment_year, data = mult.df)
summary(mod)



##with treatment information
treatment_info <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/basic dataset info/ExperimentInfo.csv")
treatment_info$trt_type <- revalue(treatment_info$trt_type, c("N*P" = "mult_nutrient"))
treatment_info$project <- treatment_info$project_name
treatment_info$community <- treatment_info$community_type
treatment_info <- treatment_info[,c("site_code", "project", "community", "plot_mani", "trt_type", "treatment", "n", "p",  "CO2", "precip", "temp")]
treatment_info <- unique(treatment_info)

#lrr_treat_species <- left_join(lrr_covariate, treatment_info, by = c("site_code", "project", "community", "plot_mani", "trt_type", "treatment"))
#lrr_treat_traits <- left_join(lrr_covariate_traits, treatment_info, by = c("site_code", "project", "community", "plot_mani", "trt_type", "treatment"))


sites.n <- subset(mean.dist.both,  trt_type == "N")%>%dplyr::select(expgroup)%>%unique()
sites.p <- subset(mean.dist.both,  trt_type == "P")%>%dplyr::select(expgroup)%>%unique()
sites.multnutrient <- subset(mean.dist.both,  trt_type == "mult_nutrient")%>%dplyr::select(expgroup)%>%unique()


n.df <- mean.dist.both%>%
        subset(expgroup %in% sites.n$expgroup)%>%
        subset(trt_type=="N"|trt_type=="control")%>%
  left_join( unite(data= treatment_info, expgroup, c("site_code", "project", "community"), sep = "::"), by = c("expgroup", "plot_mani", "trt_type", "treatment"))


p.df <- mean.dist.both%>%
  subset(expgroup %in% sites.p$expgroup)%>%
  subset(trt_type=="P"|trt_type=="control")%>%
  left_join( unite(data= treatment_info, expgroup, c("site_code", "project", "community"), sep = "::"), by = c("expgroup", "plot_mani", "trt_type", "treatment"))

mult.df <- mean.dist.both%>%
  subset(expgroup %in% sites.multnutrient$expgroup)%>%
  subset(trt_type=="mult_nutrient"|trt_type=="control")


##Nitrogen gradient

mod <- feols(mean_dist.comp~ n | site+expgroup+treatment_year, data = subset(n.df, n != 0))
summary(mod)

mod <- feols(mean_dist.trait~ n | site+expgroup+treatment_year, data = subset(n.df, n != 0))
summary(mod)

ggplot(n.df, aes(x=n, y=mean_dist.comp, color=trt_type))+
  facet_wrap(~treatment_year)+
  geom_point()+
  scale_color_manual(values = c("black", "#0099f6"))+
  geom_hline(yintercept = 0)+
  theme_base()

ggplot(n.df, aes(x=n, y=mean_dist.trait))+
  facet_wrap(~treatment_year)+
  geom_point()+
  scale_color_manual(values = c("black", "#0099f6"))+
  geom_hline(yintercept = 0)+
  theme_base()


#ggsave(
#  "C:/Users/ohler/Documents/converge-diverge/N_gradient.pdf",
#  plot = last_plot(),
#  device = "pdf",
#  path = NULL,
#  scale = 1,
#  width = 5,
#  height = 5,
#  units = c("in"),
#  dpi = 600,
#  limitsize = TRUE
#)



##P gradient

mod <- feols(mean_dist.comp~ p | site+expgroup+treatment_year, data = subset(p.df, p != 0))
summary(mod)

mod <- feols(mean_dist.trait~ p | site+expgroup+treatment_year, data = subset(p.df, p != 0))
summary(mod)

ggplot(p.df, aes(x=p, y=mean_dist.comp, color = trt_type))+
  facet_wrap(~treatment_year)+
  geom_point()+
  scale_color_manual(values = c("black", "#00b844"))+
  geom_hline(yintercept = 0)+
  theme_base()

ggplot(p.df, aes(x=p, y=mean_dist.trait, color = trt_type))+
  facet_wrap(~treatment_year)+
  geom_point()+
  scale_color_manual(values = c("black", "#00b844"))+
  geom_hline(yintercept = 0)+
  theme_base()


##mult nutrient (number of nutrients)
mod <- feols(mean_dist.comp~ plot_mani | site+expgroup+treatment_year, data = subset(mult.df, plot_mani != 0))
summary(mod)

mod <- feols(mean_dist.trait~ plot_mani | site+expgroup+treatment_year, data = subset(mult.df, plot_mani != 0))
summary(mod)

ggplot(mult.df, aes(x=plot_mani, y=mean_dist.comp, color=trt_type))+
  facet_wrap(~treatment_year)+
  geom_point()+
  scale_color_manual(values = c("black", "#6305dc"))+
  geom_hline(yintercept = 0)+
  theme_base()

ggplot(mult.df, aes(x=plot_mani, y=mean_dist.trait, color=trt_type))+
  facet_wrap(~treatment_year)+
  geom_point()+
  scale_color_manual(values = c("black", "#6305dc"))+
  geom_hline(yintercept = 0)+
  theme_base()


###BIG FIGURE PUT IN ALL TOGETHER BABY PLOT MANI FOR THE WIN

mod <- feols(mean_dist.comp~ plot_mani | site+expgroup+treatment_year, data = mean.dist.both
              # subset(mean.dist.both, plot_mani != 0)
             )
summary(mod)

x <- ggpredict(mod, c("plot_mani"))
ggplot(mean.dist.both, aes(x=plot_mani, y=mean_dist.comp))+
  #facet_wrap(~treatment_year)+
  geom_point(aes( color = trt_type), alpha = 0.1)+
  geom_smooth(data = x, aes(x, predicted), se = FALSE)+
  # scale_color_manual(values = c("black", "#00b844"))+
  geom_hline(yintercept = 0)+
  theme_base()


mod <- feols(mean_dist.trait~ plot_mani | site+expgroup+treatment_year, data = mean.dist.both
             # subset(mean.dist.both, plot_mani != 0)
)
summary(mod)


x <- ggpredict(mod, c("plot_mani"))
ggplot(mean.dist.both, aes(x=plot_mani, y=mean_dist.trait))+
  #facet_wrap(~treatment_year)+
  geom_point(aes( color = trt_type), alpha = 0.1)+
  geom_smooth(data = x, aes(x, predicted), se = FALSE)+
#  geom_smooth(method = "loess")+
  # scale_color_manual(values = c("black", "#00b844"))+
  geom_hline(yintercept = 0)+
  theme_base()


#scale_color_manual(values = c("black", "#0099f6"))+
#  scale_color_manual(values = c("black", "#00b844"))+
#  scale_color_manual(values = c("black", "#6305dc"))+















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


tempdf <- full_lrr.df#subset(full_lrr.df, trt_type == "drought" | trt_type == "N" | trt_type == "mult_nutrient")
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
#remotes::install_github("mastoffel/partR2") 
library(partR2)
library(tibble)

#tempdf <- subset(full_lrr.df, trt_type == "drought")
#mod <- lmer(lrr~sub.rich+sub.eve + sub.rank + sub.sp + (1|site_code), data = tempdf)
#summary(mod) #sig include marginal sub.rank
#r2 <- partR2(mod, data = tempdf, partvars = c("sub.rich","sub.eve" , "sub.rank" , "sub.sp"), R2_type = "marginal", nboot = 20)
#r2

#r2$R2%>%
#  subset(term == "sub.rich" | term == "sub.eve" | term == "sub.rank" | term == "sub.sp")%>%
#  mutate(term = factor(term, levels = c("sub.eve","sub.rich",  "sub.sp", "sub.rank"))) %>%
#  tibble::add_column(sig = c("0","0", "0", "0"))%>%
#  ggplot( aes(term, estimate, fill = sig))+
#  geom_bar(color = "black",stat = "identity")+
#  ylim(0,.2)+
#  ylab("Partial r-squared")+
#  scale_fill_manual(values = c( "white", "black"))+
#  coord_flip()+
#  ggtitle("DROUGHT")+
#  xlab("")+
#  theme_classic()+
#  theme(legend.position = "none")

#ggsave(
#  "C:/Users/ohler/Documents/converge-diverge/partialR2_drought.pdf",
#  plot = last_plot(),
#  device = "pdf",
#  path = NULL,
#  scale = 1,
#  width = 3,
#  height = 10,
#  units = c("in"),
#  dpi = 600,
#  limitsize = TRUE
#)


tempdf <- tempdf <- subset(full_lrr.df, trt_type == "N")%>%
  dplyr::select(lrr, sub.rich, sub.eve, sub.rank, sub.sp, expgroup)%>%
  filter(complete.cases(.))
mod <- lmer(lrr~sub.rich+sub.eve + sub.rank + sub.sp + (1|expgroup), data = tempdf)
summary(mod) #sig include sub.rank, sub.sp
r2 <- partR2(mod, data = tempdf, partvars = c("sub.rich","sub.eve" , "sub.rank" , "sub.sp"), R2_type = "marginal", nboot = 20)
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
r2 <- partR2(mod, data = tempdf, partvars = c("sub.rich","sub.eve" , "sub.rank" , "sub.sp"), R2_type = "marginal", nboot = 20)
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



tempdf <- subset(full_lrr.df, trt_type == "P")%>%
  dplyr::select(lrr, sub.rich, sub.eve, sub.rank, sub.sp, expgroup)%>%
  filter(complete.cases(.))
mod <- lmer(lrr~sub.rich+sub.eve + sub.rank + sub.sp + (1|expgroup), data = tempdf)
#mod <- lm(lrr~sub.rich+sub.eve + sub.rank + sub.sp, data = tempdf)
summary(mod) #sig include sub.rank, sub.eve
r2 <- partR2(mod, data = tempdf, partvars = c("sub.rich","sub.eve" , "sub.rank" , "sub.sp"), R2_type = "marginal", nboot = 20)
r2

r2$R2%>%
  subset(term == "sub.rich" | term == "sub.eve" | term == "sub.rank" | term == "sub.sp")%>%
  mutate(term = factor(term, levels = c("sub.sp","sub.eve",  "sub.rich", "sub.rank"))) %>%
  tibble::add_column(sig = c("0","0", "0", "0"))%>%
  ggplot( aes(term, estimate, fill = sig))+
  geom_bar(color = "black",stat = "identity")+
  ylim(0,.2)+
  ylab("Partial r-squared")+
  scale_fill_manual(values = c( "white", "black"))+
  coord_flip()+
  ggtitle("P")+
  xlab("")+
  theme_classic()+
  theme(legend.position = "none")

ggsave(
  "C:/Users/ohler/Documents/converge-diverge/partialR2_p.pdf",
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



#tempdf <- subset(full_lrr.df, trt_type == "irr")%>%
#  dplyr::select(lrr, sub.rich, sub.eve, sub.rank, sub.sp, expgroup)%>%
#  filter(complete.cases(.))
#mod <- lmer(lrr~sub.rich+sub.eve + sub.rank + sub.sp + (1|expgroup), data = tempdf)
#mod <- lm(lrr~sub.rich+sub.eve + sub.rank + sub.sp, data = tempdf)
#summary(mod) #sig include sub.rank, sub.eve
#r2 <- partR2(mod, data = tempdf, partvars = c("sub.rich","sub.eve" , "sub.rank" , "sub.sp"), R2_type = "marginal", nboot = 20)
#r2

#r2$R2%>%
#  subset(term == "sub.rich" | term == "sub.eve" | term == "sub.rank" | term == "sub.sp")%>%
#  mutate(term = factor(term, levels = c("sub.sp","sub.eve",  "sub.rich", "sub.rank"))) %>%
#  tibble::add_column(sig = c("0","0", "0", "0"))%>%
#  ggplot( aes(term, estimate, fill = sig))+
#  geom_bar(color = "black",stat = "identity")+
#  ylim(0,.2)+
#  ylab("Partial r-squared")+
#  scale_fill_manual(values = c( "white", "black"))+
#  coord_flip()+
#  ggtitle("irr")+
#  xlab("")+
#  theme_classic()+
#  theme(legend.position = "none")

#ggsave(
#  "C:/Users/ohler/Documents/converge-diverge/partialR2_irr.pdf",
#  plot = last_plot(),
#  device = "pdf",
#  path = NULL,
#  scale = 1,
#  width = 3,
#  height = 10,
#  units = c("in"),
#  dpi = 600,
#  limitsize = TRUE
#)











#library(sensemakr)
#tempdf <- subset(full_lrr.df, trt_type == "temp")%>%
#  dplyr::select(lrr, sub.rich, sub.eve, sub.rank, sub.sp, expgroup)%>%
#  filter(complete.cases(.))
#mod <- lmer(lrr~sub.rich+sub.eve + sub.rank + sub.sp + (1|expgroup), data = tempdf)
#mod <- lm(lrr~sub.rich+sub.eve + sub.rank + sub.sp, data = tempdf)
#summary(mod) #sig include sub.rank, sub.eve
#r2 <- as.data.frame(partial_r2(mod))
#r2$term <- rownames(r2)
#r2$estimate <- r2$`partial_r2(mod)`
#r2 <- partR2(mod, data = tempdf, partvars = c("sub.rich","sub.eve" , "sub.rank" , "sub.sp"), R2_type = "marginal", nboot = 10)
#r2

#r2%>%
  #subset(term == "sub.rich" | term == "sub.eve" | term == "sub.rank" | term == "sub.sp")%>%
  #mutate(term = factor(term, levels = c("sub.sp","sub.eve",  "sub.rich", "sub.rank"))) %>%
#  subset(term != "(Intercept)" )%>%
#  tibble::add_column(sig = c("0","0", "0", "0"))%>%
#  ggplot( aes(term, estimate, fill = sig))+
#  geom_bar(color = "black",stat = "identity")+
#  ylim(0,.2)+
#  ylab("Partial r-squared")+
#  scale_fill_manual(values = c( "white", "black"))+
#  coord_flip()+
#  ggtitle("TEMPERATURE")+
#  xlab("")+
#  theme_classic()+
#  theme(legend.position = "none")

ggsave(
  "C:/Users/ohler/Documents/converge-diverge/partialR2_temperature.pdf",
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


################
####SINGLE TRAIT VARIANCE AMONG REPLICATES (TRAIT)

trait_variance <- summarize.cwm%>%
  ddply(.(expgroup, trt_type, treatment, plot_mani), function(x)data.frame(
    seed_dry_mass.var = sd(x$seed_dry_mass.cwm)/mean(x$seed_dry_mass.cwm),
    LDMC.var =  sd(x$LDMC.cwm)/mean(x$LDMC.cwm),
    plant_height_vegetative.var =  sd(x$plant_height_vegetative.cwm)/mean(x$plant_height_vegetative.cwm),
    SLA.var =  sd(x$SLA.cwm)/mean(x$plant_height_vegetative.cwm),
    SRL.var = sd(x$SRL.cwm)/mean(x$SRL.cwm),
    leaf_N.var =  sd(x$leaf_N.cwm)/mean(x$leaf_N.cwm)
  ))


trt.df <- subset(trait_variance, plot_mani >= 1)%>%
  dplyr::rename(c(seed_dry_mass.var.trt = seed_dry_mass.var, 
                  LDMC.var.trt = LDMC.var, 
                  plant_height_vegetative.var.trt = plant_height_vegetative.var, 
                  SLA.var.trt = SLA.var, 
                  SRL.var.trt = SRL.var,
                  leaf_N.var.trt = leaf_N.var))
con.df <- subset(trait_variance, plot_mani == 0)%>%
  dplyr::rename(c(seed_dry_mass.var.con = seed_dry_mass.var, 
                  LDMC.var.con = LDMC.var, 
                  plant_height_vegetative.var.con = plant_height_vegetative.var, 
                  SLA.var.con = SLA.var, 
                  SRL.var.con = SRL.var,
                  leaf_N.var.con = leaf_N.var))%>%
  dplyr::select(!c(trt_type, treatment, plot_mani))

var.summary <- merge(trt.df, con.df, by = "expgroup", all.x = TRUE)

trait_variance <- tidyr::separate(trait_variance, expgroup, c("site_code", "project", "community"), sep = "::", remove = FALSE)


#Tim, start working on this chunk
#drought
#tempdf <- subset(trait_variance, trt_type == "control" |trt_type == "drought")%>%
#  subset(expgroup %in% expgroup_drought)

#mod <- lmer(seed_dry_mass.var~trt_type + (1|site_code/expgroup), data = tempdf)
#summary(mod)
#mod <- lmer(LDMC.var~trt_type + (1|site_code/expgroup), data = tempdf)
#summary(mod)
#mod <- lmer(plant_height_vegetative.var~trt_type + (1|site_code/expgroup), data = tempdf)
#summary(mod)
#mod <- lmer(SLA.var~trt_type + (1|site_code/expgroup), data = tempdf)
#summary(mod)
#mod <- lmer(SRL.var~trt_type + (1|site_code/expgroup), data = tempdf)
#summary(mod)
#mod <- lmer(leaf_N.var~trt_type + (1|site_code/expgroup), data = tempdf)
#summary(mod)


#irrigation
#tempdf <- subset(trait_variance, trt_type == "control" |trt_type == "irr")%>%
#  subset(expgroup %in% expgroup_irr)

#mod <- lmer(seed_dry_mass.var~trt_type + (1|site_code), data = tempdf)
#summary(mod)
#mod <- lmer(LDMC.var~trt_type + (1|site_code), data = tempdf)
#summary(mod)
#mod <- lmer(plant_height_vegetative.var~trt_type + (1|site_code), data = tempdf)
#summary(mod)
#mod <- lmer(SLA.var~trt_type + (1|site_code), data = tempdf)
#summary(mod)
#mod <- lmer(SRL.var~trt_type + (1|site_code), data = tempdf)
#summary(mod)
#mod <- lmer(leaf_N.var~trt_type + (1|site_code), data = tempdf)
#summary(mod)


#temperature
#tempdf <- subset(trait_variance, trt_type == "control" |trt_type == "temp")%>%
#  subset(expgroup %in% expgroup_temp)

#mod <- lmer(seed_dry_mass.var~trt_type + (1|site_code), data = tempdf)
#summary(mod)
#mod <- lmer(LDMC.var~trt_type + (1|site_code), data = tempdf)
#summary(mod)
#mod <- lmer(plant_height_vegetative.var~trt_type + (1|site_code), data = tempdf)
#summary(mod)
#mod <- lmer(SLA.var~trt_type + (1|site_code), data = tempdf)
#summary(mod)
#mod <- lmer(SRL.var~trt_type + (1|site_code), data = tempdf)
#summary(mod)
#mod <- lmer(leaf_N.var~trt_type + (1|site_code), data = tempdf)
#summary(mod)


#P
tempdf <- subset(trait_variance, trt_type == "control" |trt_type == "P")%>%
  subset(expgroup %in% expgroup_P)

mod <- lmer(seed_dry_mass.var~trt_type + (1|site_code), data = tempdf)
summary(mod)
mod <- lmer(LDMC.var~trt_type + (1|site_code), data = tempdf)
summary(mod)
mod <- lmer(plant_height_vegetative.var~trt_type + (1|site_code), data = tempdf)
summary(mod)
mod <- lmer(SLA.var~trt_type + (1|site_code), data = tempdf)
summary(mod)
mod <- lmer(SRL.var~trt_type + (1|site_code), data = tempdf)
summary(mod)
mod <- lmer(leaf_N.var~trt_type + (1|site_code), data = tempdf)
summary(mod)

#N
tempdf <- subset(trait_variance, trt_type == "control" |trt_type == "N")%>%
  subset(expgroup %in% expgroup_N)

mod <- lmer(seed_dry_mass.var~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)
mod <- lmer(LDMC.var~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)
mod <- lmer(plant_height_vegetative.var~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)
mod <- lmer(SLA.var~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)
mod <- lmer(SRL.var~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)
mod <- lmer(leaf_N.var~trt_type + (1|site_code/expgroup), data = tempdf)
summary(mod)

#mult_nutrient
tempdf <- subset(trait_variance, trt_type == "control" |trt_type == "mult_nutrient")%>%
  subset(expgroup %in% expgroup_mult_nutrient)

mod <- lmer(seed_dry_mass.var~trt_type + (1|site_code), data = tempdf)
summary(mod)
mod <- lmer(LDMC.var~trt_type + (1|site_code), data = tempdf)
summary(mod)
mod <- lmer(plant_height_vegetative.var~trt_type + (1|site_code), data = tempdf)
summary(mod)
mod <- lmer(SLA.var~trt_type + (1|site_code), data = tempdf)
summary(mod)
mod <- lmer(SRL.var~trt_type + (1|site_code), data = tempdf)
summary(mod)
mod <- lmer(leaf_N.var~trt_type + (1|site_code), data = tempdf)
summary(mod)


##

#mod <- lmer(seed_dry_mass.var~trt_type + (1|site_code/expgroup), data = trait_variance)
#summary(mod)


#mod <- lmer(LDMC.var~trt_type + (1|site_code/expgroup), data = trait_variance)
#summary(mod)

#mod <- lmer(plant_height_vegetative.var~trt_type + (1|site_code/expgroup), data = trait_variance)
#summary(mod)

#mod <- lmer(SLA.var~trt_type + (1|expgroup), data = trait_variance)
#summary(mod)

#mod <- lmer(rooting_depth.var~trt_type + (1|expgroup), data = trait_variance)
#summary(mod)

#mod <- lmer(leaf_C.N.var~trt_type + (1|expgroup), data = trait_variance)
#summary(mod)

