################################################################################
##  5_sCoRRE_CWMtraits.R: Examining community weighted mean trait responses to treatments in the CoRRE database.
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
fDiv <- read.csv('CoRRE_allDiversityMetrics_phyFunAnalysis.csv')


##### calculate continuous response ratios #####
#filter control plots
controlCont <- fDiv %>% 
  filter(trt_type2=='control') %>%
  group_by(site_code, project_name, community_type, treatment_year) %>%
  summarize_at(vars(CWM.leaf_C.N, CWM.LDMC, CWM.SLA, CWM.plant_height_vegetative, CWM.rooting_depth, CWM.seed_dry_mass), list(ctl=mean), na.rm=T) %>% #average across plots and years
  ungroup()

allDivRRCont <- fDiv %>%
  filter(trt_type2!='control') %>%
  left_join(controlCont) %>%
  mutate(CWM.leaf_C.N_RR=(CWM.leaf_C.N-CWM.leaf_C.N_ctl)/CWM.leaf_C.N_ctl, CWM.LDMC_RR=(CWM.LDMC-CWM.LDMC_ctl)/CWM.LDMC_ctl, CWM.SLA_RR=(CWM.SLA-CWM.SLA_ctl)/CWM.SLA_ctl, CWM.plant_height_vegetative_RR=(CWM.plant_height_vegetative-CWM.plant_height_vegetative_ctl)/CWM.plant_height_vegetative_ctl, CWM.rooting_depth_RR=(CWM.rooting_depth-CWM.rooting_depth_ctl)/CWM.rooting_depth_ctl, CWM.seed_dry_mass_RR=(CWM.seed_dry_mass-CWM.seed_dry_mass_ctl)/CWM.seed_dry_mass_ctl) %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::')) %>% 
  group_by(site_proj_comm, site_code, project_name, community_type, treatment, trt_type2, plot_id) %>%
  summarise_at(vars(CWM.leaf_C.N_RR, CWM.LDMC_RR, CWM.SLA_RR, CWM.plant_height_vegetative_RR, CWM.rooting_depth_RR, CWM.seed_dry_mass_RR), list(mean=mean), na.rm=T) %>%
  ungroup() %>% 
  filter(!is.nan(CWM.leaf_C.N_RR_mean))

##### checking continuous trait data #####
hist(allDivRRCont$CWM.leaf_C.N_RR_mean)
qqPlot(allDivRRCont$CWM.leaf_C.N_RR_mean)

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


##### continuous trait mixed models #####
library(nlme)
library(emmeans)
library(sjPlot)
library(performance)
library(grid)

options(contrasts=c('contr.sum','contr.poly')) 

#leaf C:N model
summary(cnModel <- lme(CWM.leaf_C.N_RR_mean ~ as.factor(trt_type2),
                        data=na.omit(subset(allDivRRCont, trt_type2!='herb_removal')),
                        random=~1|site_proj_comm))
anova.lme(cnModel, type='sequential')
meansCNModel <- emmeans(cnModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansCNModelOutput <- as.data.frame(meansCNModel$emmeans)

cnFig <- ggplot(data=meansCNModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('CWM SLA\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')

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
print(cnFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(seedFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(slaFig, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(ldmcFig, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(heightFig, vp=viewport(layout.pos.row=3, layout.pos.col=1))
print(rootFig, vp=viewport(layout.pos.row=3, layout.pos.col=2))
#export at 1000x1500


##### calculate categorical trait responses #####
#filter control plots
controlCat <- fDiv %>% 
  filter(trt_type2=='control') %>%
  group_by(site_code, project_name, community_type, treatment_year) %>%
  summarize_at(vars(CWM.growth_form, CWM.photosynthetic_pathway, CWM.lifespan, CWM.clonal, CWM.mycorrhizal_type, CWM.n_fixation), list(ctl=max), na.rm=T) %>% #average aross plots and years
  ungroup()

allDivRRCat <- fDiv %>% 
  filter(!(trt_type2=='control')) %>% 
  left_join(controlCat) %>% 
  select(site_code, project_name, community_type, site_proj_comm, treatment_year, trt_type2, treatment, plot_id, CWM.growth_form:CWM.n_fixation, CWM.growth_form_ctl:CWM.n_fixation_ctl) %>% 
  rename(CWM.growthform_trt=CWM.growth_form, 
         CWM.growthform_ctl=CWM.growth_form_ctl,
         CWM.photosyntheticpathway_trt=CWM.photosynthetic_pathway, 
         CWM.photosyntheticpathway_ctl=CWM.photosynthetic_pathway_ctl,
         CWM.lifespan_trt=CWM.lifespan,
         CWM.clonal_trt=CWM.clonal,
         CWM.mycorrhizaltype_trt=CWM.mycorrhizal_type, 
         CWM.mycorrhizaltype_ctl=CWM.mycorrhizal_type_ctl,
         CWM.nfixation_trt=CWM.n_fixation,
         CWM.nfixation_ctl=CWM.n_fixation_ctl) %>% 
  group_by(site_code, project_name, community_type, site_proj_comm, trt_type2, treatment, plot_id) %>% 
  summarize_at(vars(CWM.growthform_trt, CWM.photosyntheticpathway_trt, CWM.lifespan_trt, CWM.clonal_trt, CWM.mycorrhizaltype_trt, CWM.nfixation_trt, CWM.growthform_ctl, CWM.photosyntheticpathway_ctl, CWM.lifespan_ctl, CWM.clonal_ctl, CWM.mycorrhizaltype_ctl, CWM.nfixation_ctl), list(max), na.rm=T) %>% #average aross plots and years
  ungroup()

##### individual categorical trait figures #####
### growth form ###
growthform <- allDivRRCat %>% 
  select(site_code, project_name, community_type, site_proj_comm, trt_type2, treatment, plot_id, CWM.growthform_trt, CWM.growthform_ctl)

#N
growthformN <- growthform %>% 
  filter(trt_type2=='N') %>% 
  rename(trt=CWM.growthform_trt, ctl=CWM.growthform_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='graminoid', 1, 0))

#P
growthformP <- growthform %>% 
  filter(trt_type2=='P') %>% 
  rename(trt=CWM.growthform_trt, ctl=CWM.growthform_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='graminoid', 1, 0))

#irr
growthformIrr <- growthform %>% 
  filter(trt_type2=='irr') %>% 
  rename(trt=CWM.growthform_trt, ctl=CWM.growthform_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='graminoid', 1, 0))

#CO2
growthformCO2 <- growthform %>% 
  filter(trt_type2=='CO2') %>% 
  rename(trt=CWM.growthform_trt, ctl=CWM.growthform_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='graminoid', 1, 0))

#drought
growthformDrought <- growthform %>% 
  filter(trt_type2=='drought') %>% 
  rename(trt=CWM.growthform_trt, ctl=CWM.growthform_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='graminoid', 1, 0))

#temp
growthformTemp <- growthform %>% 
  filter(trt_type2=='temp') %>% 
  rename(trt=CWM.growthform_trt, ctl=CWM.growthform_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='graminoid', 1, 0))

#disturbance
growthformDisturbance <- growthform %>% 
  filter(trt_type2=='disturbance') %>% 
  rename(trt=CWM.growthform_trt, ctl=CWM.growthform_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='graminoid', 1, 0))

#multiple trts
growthformMultiple <- growthform %>% 
  filter(trt_type2=='multiple trts') %>% 
  rename(trt=CWM.growthform_trt, ctl=CWM.growthform_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='graminoid', 1, 0))


#rbind for all growth form
growthform2 <- rbind(growthformN, growthformP, growthformIrr, growthformCO2, growthformDrought, growthformTemp, growthformDisturbance, growthformMultiple)

#figure
growthformFig <- ggplot(data=growthform2, aes(x=trt_type2, y=value2, color=trtctl)) +
  geom_jitter(height=0.1) +
  stat_summary(
    geom = "linerange",
    fun.data = mean_sdl, 
    fun.args = list(mult = 1),
    aes(colour = trtctl),
    position=position_dodge(0.4)
  ) +
  stat_summary(
    geom = "point",
    fun = mean,
    aes(colour = trtctl), size = 3,
    position=position_dodge(0.4)
  ) +
  ylab('Growth Form') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_color_manual(values=c('grey', 'dark orange'))

  
### lifespan ###
lifespan <- allDivRRCat %>% 
  select(site_code, project_name, community_type, site_proj_comm, trt_type2, treatment, plot_id, CWM.lifespan_trt, CWM.lifespan_ctl)

#N
lifespanN <- lifespan %>% 
  filter(trt_type2=='N') %>% 
  rename(trt=CWM.lifespan_trt, ctl=CWM.lifespan_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='perennial', 1, 0))

#P
lifespanP <- lifespan %>% 
  filter(trt_type2=='P') %>% 
  rename(trt=CWM.lifespan_trt, ctl=CWM.lifespan_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='perennial', 1, 0))

#irr
lifespanIrr <- lifespan %>% 
  filter(trt_type2=='irr') %>% 
  rename(trt=CWM.lifespan_trt, ctl=CWM.lifespan_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='perennial', 1, 0))

#CO2
lifespanCO2 <- lifespan %>% 
  filter(trt_type2=='CO2') %>% 
  rename(trt=CWM.lifespan_trt, ctl=CWM.lifespan_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='perennial', 1, 0))

#drought
lifespanDrought <- lifespan %>% 
  filter(trt_type2=='drought') %>% 
  rename(trt=CWM.lifespan_trt, ctl=CWM.lifespan_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='perennial', 1, 0))

#temp
lifespanTemp <- lifespan %>% 
  filter(trt_type2=='temp') %>% 
  rename(trt=CWM.lifespan_trt, ctl=CWM.lifespan_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='perennial', 1, 0))

#disturbance
lifespanDisturbance <- lifespan %>% 
  filter(trt_type2=='disturbance') %>% 
  rename(trt=CWM.lifespan_trt, ctl=CWM.lifespan_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='perennial', 1, 0))

#multiple trts
lifespanMultiple <- lifespan %>% 
  filter(trt_type2=='multiple trts') %>% 
  rename(trt=CWM.lifespan_trt, ctl=CWM.lifespan_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='perennial', 1, 0))


#rbind for all lifespan
lifespan2 <- rbind(lifespanN, lifespanP, lifespanIrr, lifespanCO2, lifespanDrought, lifespanTemp, lifespanDisturbance, lifespanMultiple)

#figure
lifespanFig <- ggplot(data=lifespan2, aes(x=trt_type2, y=value2, color=trtctl)) +
  geom_jitter(height=0.1) +
  stat_summary(
    geom = "linerange",
    fun.data = mean_sdl, 
    fun.args = list(mult = 1),
    aes(colour = trtctl),
    position=position_dodge(0.4)
  ) +
  stat_summary(
    geom = "point",
    fun = mean,
    aes(colour = trtctl), size = 3,
    position=position_dodge(0.4)
  ) +
  ylab('Lifespan') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_color_manual(values=c('grey', 'dark orange'))

  
### mycorrhizaltype ###
mycorrhizaltype <- allDivRRCat %>% 
  select(site_code, project_name, community_type, site_proj_comm, trt_type2, treatment, plot_id, CWM.mycorrhizaltype_trt, CWM.mycorrhizaltype_ctl)

#N
mycorrhizaltypeN <- mycorrhizaltype %>% 
  filter(trt_type2=='N') %>% 
  rename(trt=CWM.mycorrhizaltype_trt, ctl=CWM.mycorrhizaltype_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='yes', 1, 0))

#P
mycorrhizaltypeP <- mycorrhizaltype %>% 
  filter(trt_type2=='P') %>% 
  rename(trt=CWM.mycorrhizaltype_trt, ctl=CWM.mycorrhizaltype_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='yes', 1, 0))

#irr
mycorrhizaltypeIrr <- mycorrhizaltype %>% 
  filter(trt_type2=='irr') %>% 
  rename(trt=CWM.mycorrhizaltype_trt, ctl=CWM.mycorrhizaltype_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='yes', 1, 0))

#CO2
mycorrhizaltypeCO2 <- mycorrhizaltype %>% 
  filter(trt_type2=='CO2') %>% 
  rename(trt=CWM.mycorrhizaltype_trt, ctl=CWM.mycorrhizaltype_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='yes', 1, 0))

#drought
mycorrhizaltypeDrought <- mycorrhizaltype %>% 
  filter(trt_type2=='drought') %>% 
  rename(trt=CWM.mycorrhizaltype_trt, ctl=CWM.mycorrhizaltype_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='yes', 1, 0))

#temp
mycorrhizaltypeTemp <- mycorrhizaltype %>% 
  filter(trt_type2=='temp') %>% 
  rename(trt=CWM.mycorrhizaltype_trt, ctl=CWM.mycorrhizaltype_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='yes', 1, 0))

#disturbance
mycorrhizaltypeDisturbance <- mycorrhizaltype %>% 
  filter(trt_type2=='disturbance') %>% 
  rename(trt=CWM.mycorrhizaltype_trt, ctl=CWM.mycorrhizaltype_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='yes', 1, 0))

#multiple trts
mycorrhizaltypeMultiple <- mycorrhizaltype %>% 
  filter(trt_type2=='multiple trts') %>% 
  rename(trt=CWM.mycorrhizaltype_trt, ctl=CWM.mycorrhizaltype_ctl) %>% 
  pivot_longer(cols=trt:ctl, names_to='trtctl') %>% 
  mutate(value2=ifelse(value=='yes', 1, 0))


#rbind for all mycorrhizaltype
mycorrhizaltype2 <- rbind(mycorrhizaltypeN, mycorrhizaltypeP, mycorrhizaltypeIrr, mycorrhizaltypeCO2, mycorrhizaltypeDrought, mycorrhizaltypeTemp, mycorrhizaltypeDisturbance, mycorrhizaltypeMultiple)

#figure
mycorrhizaltypeFig <- ggplot(data=mycorrhizaltype2, aes(x=trt_type2, y=value2, color=trtctl)) +
  geom_jitter(height=0.1) +
  stat_summary(
    geom = "linerange",
    fun.data = mean_sdl, 
    fun.args = list(mult = 1),
    aes(colour = trtctl),
    position=position_dodge(0.4)
  ) +
  stat_summary(
    geom = "point",
    fun = mean,
    aes(colour = trtctl), size = 3,
    position=position_dodge(0.4)
  ) +
  ylab('mycorrhizaltype') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_color_manual(values=c('grey', 'dark orange'))

  
#didn't finish making these figs, there's three more traits to go
pushViewport(viewport(layout=grid.layout(1,3)))
print(growthformFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(lifespanFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(mycorrhizaltypeFig, vp=viewport(layout.pos.row=1, layout.pos.col=3))
#export 1800x800
