################################################################################
##  7_sCoRRE_phylogeneticDiversity_CWMtraits.R: Examining community weighted mean trait responses to treatments in the CoRRE database.
##
##  Author: Kimberly Komatsu
##  Date created: December 13, 2022
################################################################################

library(data.table)
library(codyn)
library(fixest)
library(lme4)
library(tidyverse)

#set working directory
setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\paper 2_PD and FD responses\\data')


#functional diversity data
fDiv <- read.csv('CoRRE_functionalDiversity_2022-12-13.csv')

#treatment data
trt <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\RawAbundance.csv') %>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id) %>%
  unique() %>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\ExperimentInfo.csv')) %>%
  group_by(site_code, project_name, community_type) %>%
  mutate(experiment_length=max(treatment_year)) %>%
  ungroup() %>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, trt_type, experiment_length, plot_mani, n, p, CO2, precip, temp)

#merge CWM and trt data
allDiv <- fDiv %>% #functional metrics
  select(plot_id, treatment_year, site_code, project_name, community_type, site_proj_comm, treatment, CWM.growth_form, CWM.photosynthetic_pathway, CWM.lifespan, CWM.clonal, CWM.mycorrhizal_type, CWM.n_fixation, CWM.leaf_C.N, CWM.LDMC, CWM.SLA, CWM.plant_height_vegetative, CWM.rooting_depth, CWM.seed_dry_mass) %>% 
  left_join(trt) %>% #treatments
  full_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteBiotic.csv')) %>% #site anpp and regional richness
  full_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteLocationClimate.csv')) %>% #site MAP and MAT
  mutate(site_proj_comm_trt=paste(site_proj_comm, treatment, sep='::')) %>% 
  select(-PubLat, -PubLong, -Offset, -Location) %>%
  filter(!(site_code %in% c('SERC', 'CAR', 'PIE', 'NANT'))) #remove wetlands, which have very few species and therefore don't nicely fit into these response types (about 2000 data points)

#selecting relevant treatments for analysis (high resource, high stress)
trt_analysis <- trt %>%
  mutate(alltrts=ifelse(trt_type %in% c("control", "CO2","CO2*temp", "mow_clip","burn","burn*graze","disturbance","burn*mow_clip","drought","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N","mult_nutrient","N*P","P","N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp","temp*mow_clip","drought*temp","irr*temp","irr"),1,0)) %>%
  filter(alltrts==1) %>%
  mutate(dist=ifelse(trt_type %in% c("mow_clip","burn","burn*graze","disturbance","burn*mow_clip"), 1, 0), #unify codes across datasets
         # tCO2=ifelse(trt_type %in% c("CO2"), 1, 0),
         drought=ifelse(trt_type %in% c("drought"), 1, 0),
         # therb_removal=ifelse(trt_type %in% c("herb_removal"), 1, 0),
         irg=ifelse(trt_type %in% c("irr"), 1, 0),
         # ttemp=ifelse(trt_type %in% c("temp"), 1, 0),
         # tn=ifelse(trt_type %in% c("N"), 1, 0),
         # tp=ifelse(trt_type %in% c("P"), 1, 0),
         multtrts=ifelse(trt_type %in% c("CO2*temp", "burn*graze","burn*mow_clip","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp*mow_clip","drought*temp","irr*temp","mult_nutrient","N*P"),1,0)) %>%
  mutate(trt_type2=ifelse(dist==1, 'disturbance', ifelse(multtrts==1, 'multiple trts', trt_type))) %>%
  select(site_code, project_name, community_type, treatment, alltrts, dist, drought, irg, multtrts, trt_type2) %>%
  unique()

#pick treatments here
allDivTrt <- allDiv %>%
  right_join(trt_analysis) %>%
  mutate(trt_binary=ifelse(plot_mani>0, 1, 0)) %>%
  filter(!(is.na(site_proj_comm)))

# write.csv(allDivTrt, 'CoRRE_CWMtraits_12142022.csv')


##### calculate response ratios #####
#filter control plots
controlCont <- allDivTrt %>% 
  filter(trt_type=='control') %>%
  group_by(site_code, project_name, community_type, treatment_year) %>%
  summarize_at(vars(CWM.LDMC, CWM.SLA, CWM.plant_height_vegetative, CWM.rooting_depth, CWM.seed_dry_mass), list(ctl=mean), na.rm=T) %>% #average aross plots and years
  ungroup()

controlCat <- allDivTrt %>% 
  filter(trt_type=='control') %>%
  group_by(site_code, project_name, community_type, treatment_year) %>%
  summarize_at(vars(CWM.growth_form, CWM.photosynthetic_pathway, CWM.lifespan, CWM.clonal, CWM.mycorrhizal_type, CWM.n_fixation), list(ctl=max), na.rm=T) %>% #average aross plots and years
  ungroup()

control <- controlCont %>% 
  left_join(controlCat)

allDivRRCont <- allDivTrt %>%
  filter(trt_type!='control') %>%
  left_join(controlCont) %>%
  mutate(CWM.LDMC_RR=(CWM.LDMC-CWM.LDMC_ctl)/CWM.LDMC_ctl, CWM.SLA_RR=(CWM.SLA-CWM.SLA_ctl)/CWM.SLA_ctl, CWM.plant_height_vegetative_RR=(CWM.plant_height_vegetative-CWM.plant_height_vegetative_ctl)/CWM.plant_height_vegetative_ctl, CWM.rooting_depth_RR=(CWM.rooting_depth-CWM.rooting_depth_ctl)/CWM.rooting_depth_ctl, CWM.seed_dry_mass_RR=(CWM.seed_dry_mass-CWM.seed_dry_mass_ctl)/CWM.seed_dry_mass_ctl) %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::')) %>% 
  group_by(site_proj_comm, site_code, project_name, community_type, treatment, trt_type2, plot_id) %>%
  summarise_at(vars(CWM.LDMC_RR, CWM.SLA_RR, CWM.plant_height_vegetative_RR, CWM.rooting_depth_RR, CWM.seed_dry_mass_RR), list(mean=mean), na.rm=T) %>%
  ungroup()

allDivRRCat <- allDivTrt 


##### checking data #####
hist(allDivRRCont$CWM.LDMC_RR_mean)
qqPlot(allDivRRCont$CWM.LDMC_RR_mean)

hist(allDivRRCont$CWM.SLA_RR_mean)
qqPlot(allDivRRCont$CWM.SLA_RR_mean)

hist(allDivRRCont$CWM.plant_height_vegetative_RR_mean)
qqPlot(allDivRRCont$CWM.plant_height_vegetative_RR_mean)

hist(allDivRRCont$CWM.rooting_depth_RR_mean)
qqPlot(allDivRRCont$CWM.rooting_depth_RR_mean)

hist(allDivRRCont$CWM.seed_dry_mass_RR_mean)
qqPlot(allDivRRCont$CWM.seed_dry_mass_RR_mean)


##### mixed models #####
library(nlme)
library(emmeans)
library(sjPlot)
library(performance)
library(grid)

options(contrasts=c('contr.sum','contr.poly')) 

#LDMC model
summary(ldmcModel <- lme(CWM.LDMC_RR_mean ~ as.factor(trt_type2),
                         data=na.omit(subset(allDivRRCont, trt_type2!='herb_removal')),
                         random=~1|site_proj_comm))
anova.lme(ldmcModel, type='sequential')
meansLDMCModel <- emmeans(ldmcModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansLDMCModelOutput <- as.data.frame(meansLDMCModel$emmeans)

ldmcFig <- ggplot(data=meansLDMCModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('CWM LDMC\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')

#SLA model
summary(slaModel <- lme(CWM.SLA_RR_mean ~ as.factor(trt_type2),
                         data=na.omit(subset(allDivRRCont, trt_type2!='herb_removal')),
                         random=~1|site_proj_comm))
anova.lme(slaModel, type='sequential')
meansSLAModel <- emmeans(slaModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansSLAModelOutput <- as.data.frame(meansSLAModel$emmeans)

slaFig <- ggplot(data=meansSLAModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('CWM SLA\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')

#height model
summary(heightModel <- lme(CWM.plant_height_vegetative_RR_mean ~ as.factor(trt_type2),
                        data=na.omit(subset(allDivRRCont, trt_type2!='herb_removal')),
                        random=~1|site_proj_comm))
anova.lme(heightModel, type='sequential')
meansHeightModel <- emmeans(heightModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansHeightModelOutput <- as.data.frame(meansHeightModel$emmeans)

heightFig <- ggplot(data=meansHeightModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('CWM Height\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')

#rooting depth model
summary(rootModel <- lme(CWM.rooting_depth_RR_mean ~ as.factor(trt_type2),
                           data=na.omit(subset(allDivRRCont, trt_type2!='herb_removal')),
                           random=~1|site_proj_comm))
anova.lme(rootModel, type='sequential')
meansRootModel <- emmeans(rootModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansRootModelOutput <- as.data.frame(meansRootModel$emmeans)

rootFig <- ggplot(data=meansRootModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('CWM Rooting Depth\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')

#seed mass model
summary(seedModel <- lme(CWM.seed_dry_mass_RR_mean ~ as.factor(trt_type2),
                         data=na.omit(subset(allDivRRCont, trt_type2!='herb_removal')),
                         random=~1|site_proj_comm))
anova.lme(seedModel, type='sequential')
meansSeedModel <- emmeans(seedModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansSeedModelOutput <- as.data.frame(meansSeedModel$emmeans)

seedFig <- ggplot(data=meansSeedModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('CWM Seed Dry Mass\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


pushViewport(viewport(layout=grid.layout(3,2)))
print(slaFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(ldmcFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(heightFig, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(rootFig, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(seedFig, vp=viewport(layout.pos.row=3, layout.pos.col=1))
#export at 1000x1500