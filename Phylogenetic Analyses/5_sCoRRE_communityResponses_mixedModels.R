################################################################################
##  sCoRRE_communityResponses_mixedModels.R: Examining differences in
##  phylogenetic and functional diversity within the CoRRE database.
##
##  Author: Kimberly Komatsu
##  Date created: December 13, 2021
################################################################################

library(PerformanceAnalytics)
library(nlme)
library(emmeans)
library(sjPlot)
library(performance)
library(grid)
library(data.table)
library(fixest)
library(tidyverse)


setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared')  #kim's computer

##### functions and themes #####
###standard error function
se <- function(x, na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  

#model summary table (from Peter Wilfahrt)
summary.tablefunc <- function(mod) {  
  dat <- data.frame(summary(mod)$tTable) %>%
    tibble::rownames_to_column(var = 'Effect') %>%
    rename_with(stringr::str_replace, 
                pattern = "-", replacement = ".", 
                matches("Length")) %>% 
    dplyr::mutate(Estimate = signif(Value, digits = 3),
                  Std.Error = signif(Std.Error, digits = 2),
                  t.value = signif(t.value, digits = 2),
                  p.value = signif(p.value, digits = 2)) %>%
    dplyr::mutate(p.value = ifelse(p.value <= 0.001, '< 0.001', p.value)) %>% 
    dplyr::select(-Value) %>% 
    relocate(Estimate,.before = Std.Error)
  return(dat)
}


#theme set
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))



##### data #####
#treatment data
trt <- read.csv('CoRRE data\\CoRRE data\\community composition\\CoRRE_RawAbundance_Jan2023.csv') %>%
select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id) %>%
  unique() %>%
  left_join(read.csv('CoRRE data\\CoRRE data\\basic dataset info\\ExperimentInfo.csv')) %>%
  group_by(site_code, project_name, community_type) %>%
  mutate(experiment_length=max(treatment_year)) %>%
  ungroup() %>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, trt_type, experiment_length, plot_mani, n, p, CO2, precip, temp)

#phylogenetic diversity data
pDiv <- read.csv('paper 2_PD and FD responses\\data\\CoRRE_pd_metrics_non_weighted.csv') %>%
  separate(identifier, into=c("site_code", "project_name", "community_type", "treatment_year", "plot_id"), sep="::") %>%
  mutate(treatment_year=as.integer(treatment_year))

#functional diversity data
fDiv <- read.csv('paper 2_PD and FD responses\\data\\CoRRE_functionalDiversity_2023-03-21.csv')

#taxonomic diversity data
rDiv <- read.csv('paper 2_PD and FD responses\\data\\CoRRE_taxonomicDiversity_2023-03-21.csv')


#merge all data on diversity metrics (phylogenetic, functional, species), experimental treatments, and site characteristics
allDiv <- pDiv %>% #phylogenetic metrics
  full_join(fDiv) %>% #functional metrics
  full_join(rDiv) %>% #species metrics
  filter(treatment_year>0) %>% #removing pre-treatment data
  left_join(trt) %>% #treatments
  full_join(read.csv('CoRRE data\\CoRRE data\\environmental data\\CoRRE_siteBiotic_Dec2021.csv')) %>% #site anpp and regional richness
  full_join(read.csv('CoRRE data\\CoRRE data\\environmental data\\CoRRE_siteLocationClimate_July2022.csv')) %>% #site MAP and MAT
  mutate(site_proj_comm=paste(site_code,  project_name, community_type, sep='::')) %>%
  mutate(site_proj_comm_trt=paste(site_proj_comm, treatment, sep='::')) %>% 
  select(-pd.ses, -pd.raw, -pd.pval, -nbsp, -sing.sp, -FEve, -FDiv, -block, -data_type, -version, -Evar, -trt_type, -experiment_length, -Location, -Continent, -PubLat, -PubLong, -Offset) %>%
  filter(richness>3) #filter plots with richness less than 4 because they have volatile values for functional and phylogenetic diversity responses

##### determine which sites don't have adequate cover to species for these analyses #####
#species relative cover data
relCover <- read.csv('CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeCover_Jan2023.csv') %>%
  mutate(replicate=paste(site_code, project_name, community_type, plot_id, sep='::')) #creating identifying column of each plot

corre_to_try <- read.csv("CoRRE data\\trait data\\corre2trykey_2021.csv") %>%
  dplyr::select(genus_species, species_matched) %>%
  unique(.)

exptN <- relCover %>% 
  select(site_code, project_name, community_type, treatment, treatment_year) %>% 
  unique() %>%
  group_by(site_code, project_name, community_type, treatment) %>% 
  summarise(length=length(treatment_year)) %>% 
  ungroup() 

lowTraitCover <- relCover %>% 
  left_join(corre_to_try) %>% 
  left_join(trt) %>% 
  mutate(unknown_spp=ifelse(is.na(species_matched), 1, 0)) %>% 
  group_by(site_code, project_name, community_type, trt_type, treatment, treatment_year, plot_id, unknown_spp) %>% 
  summarise(relcover=sum(relcov)) %>% 
  ungroup() %>% 
  filter(unknown_spp==1, relcover>0.2) %>% 
  select(-relcover, -trt_type)
#2649 plot*yr combinations will be dropped from analysis for too little trait data going into calculations


##### selecting relevant treatments for analysis (high resource, high stress) #####
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
  mutate(drop=ifelse(site_code=="CDR"&treatment %in% c(2, 3, 4, 5, 7), 1, 0)) %>% #drop some of the CDR e001 and e002 treatments to prevent over-representation
  filter(drop==0) %>%
  filter(!(site_proj_comm %in% c('DCGS::gap::0'))) %>% #remove pulse light expt
  mutate_at(c('MAP', 'MAT', 'rrich', 'anpp'), funs(c(scale(.)))) %>% #scale site characteristics
  filter(!(is.na(site_proj_comm))) %>% 
  left_join(lowTraitCover) %>% 
  filter(is.na(unknown_spp), RaoQ>0) %>% #filter where there are too many unknown species to have trait data for plot, removes 2137 data points
  group_by(site_code, project_name, community_type, treatment, treatment_year, trt_type2) %>% 
  mutate(length=length(unique(plot_id))) %>% 
  ungroup() %>% 
  filter(length>2) #filter out any experimental treatment where N<3 for a given year, removes 673 data points


# write.csv(allDivTrt, 'paper 2_PD and FD responses\\data\\CoRRE_allDiversityMetrics_phyFunAnalysis.csv', row.names=F)


##### calculating response ratios #####
#filter control plots
control <- allDivTrt %>% 
  filter(trt_type2=='control') %>%
  rename(mpd.raw_ctl=mpd.raw, 
         mntd.raw_ctl=mntd.raw, 
         FDis_ctl=FDis, 
         richness_ctl=richness, 
         RaoQ_ctl=RaoQ, 
         MNTD_traits_ctl=MNTD_traits,
         hill_ctl=hill) %>%
  group_by(site_code, project_name, community_type, treatment_year) %>%
  summarize_at(vars(mpd.raw_ctl, 
                    mntd.raw_ctl, 
                    FDis_ctl, 
                    RaoQ_ctl, 
                    richness_ctl, 
                    MNTD_traits_ctl,
                    hill_ctl), list(mean=mean), na.rm=T) %>% #average across plots and years
  ungroup()

#merge on site characteristics
controlEnv <- control %>%
  group_by(site_code, project_name, community_type) %>%
  summarize_at(vars(mpd.raw_ctl_mean, mntd.raw_ctl_mean, FDis_ctl_mean, RaoQ_ctl_mean, richness_ctl_mean, MNTD_traits_ctl_mean, hill_ctl_mean), list(mean=mean), na.rm=T) %>% #average across plots and years
  ungroup() %>%
  left_join(read.csv('CoRRE data\\CoRRE data\\environmental data\\CoRRE_siteBiotic_Dec2021.csv')) %>%
  left_join(read.csv('CoRRE data\\CoRRE data\\environmental data\\CoRRE_siteLocationClimate_July2022.csv')) %>%
  gather(key='env_variable', value='env_value', rrich, anpp, MAP, MAT, aridityValues)

allDivRR <- allDivTrt %>%
  filter(trt_type2!='control') %>%
  left_join(control) %>%
  filter(!is.na(RaoQ_ctl_mean)) %>%  #remove lines where there was no control to compare to due to lack of spp cover for traits; lose 550 lines
  mutate(mpd_RR=log(mpd.raw/mpd.raw_ctl_mean), 
         mntd_RR=log(mntd.raw/mntd.raw_ctl_mean), 
         FDis_RR=log(FDis/(FDis_ctl_mean)), 
         RaoQ_RR=log(RaoQ/(RaoQ_ctl_mean)), 
         richness_RR=log(richness/(richness_ctl_mean)),
         hill_RR=log(hill/(hill_ctl_mean)),
         MNTD_traits_RR=log(MNTD_traits/(MNTD_traits_ctl_mean))) %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::')) %>% 
  select(site_proj_comm, site_code, project_name, community_type, treatment_year, calendar_year, treatment, trt_type2, plot_mani, plot_id, rrich, anpp, MAP, MAT, aridityValues, mpd_RR:MNTD_traits_RR)

allDivRRmean <- allDivRR %>% 
  group_by(site_proj_comm, site_code, project_name, community_type, treatment, trt_type2, plot_id) %>%
  summarise_at(vars(mpd_RR, mntd_RR, FDis_RR, RaoQ_RR, richness_RR, hill_RR, MNTD_traits_RR), list(mean=mean), na.rm=T) %>%
  ungroup()


##### check if richness difference is correlated with other diversity metrics #####

chart.Correlation(allDivTrt[c(6,7,9,10,12,13,29,30,31)]) #raw values

chart.Correlation(allDivRR[16:22]) #response ratios

chart.Correlation(allDivRRmean[8:14]) #mean response ratios over all years


##### mixed effects model #####
#NOTE: these models do not account for biotic or abiotic env drivers at a site or for trt magnitude (but do include a random effect of site)

options(contrasts=c('contr.sum','contr.poly')) 

summary(richModel <- lme(richness_RR_mean ~ as.factor(trt_type2),
                         data=na.omit(subset(allDivRRmean, trt_type2!='herb_removal')),
                         random=~1|site_proj_comm))
anova.lme(richModel, type='sequential') #significant effect of trt
meansRichModel <- emmeans(richModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansRichModelOutput <- as.data.frame(meansRichModel$emmeans)
plot_model(richModel, type = "pred", terms = c("trt_type2", "richness_RR_mean"))

richFig <- ggplot(data=meansRichModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0, size=3) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('Species Richness\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


summary(HillModel <- lme(hill_RR_mean ~ as.factor(trt_type2),
                         data=na.omit(subset(allDivRRmean, trt_type2!='herb_removal')),
                         random=~1|site_proj_comm))
anova.lme(HillModel, type='sequential') #significant trt effect
meansHillModel <- emmeans(HillModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansHillModelOutput <- as.data.frame(meansHillModel$emmeans)
plot_model(HillModel, type = "pred", terms = c("trt_type2", "hill_RR_mean"))

hillFig <- ggplot(data=meansHillModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0, size=3) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('Hill Number\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


summary(mpdModel <- lme(mpd_RR_mean ~ as.factor(trt_type2) + richness_RR_mean,
                         data=na.omit(subset(allDivRRmean, trt_type2!='herb_removal')),
                         random=~1|site_proj_comm))
anova.lme(mpdModel, type='sequential') #significant effect of trt
meansMPDModel <- emmeans(mpdModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansMPDModelOutput <- as.data.frame(meansMPDModel$emmeans)
plot_model(mpdModel, type = "pred", terms = c("trt_type2", "richness_RR_mean"))

mpdFig <- ggplot(data=meansMPDModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0, size=3) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('MPD\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


summary(RaoQModel <- lme(RaoQ_RR_mean ~ as.factor(trt_type2) + richness_RR_mean,
                         data=na.omit(subset(allDivRRmean, trt_type2!='herb_removal')),
                         random=~1|site_proj_comm))
anova.lme(RaoQModel, type='sequential') #significant effect of trt
meansRaoQModel <- emmeans(RaoQModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansRaoQModelOutput <- as.data.frame(meansRaoQModel$emmeans)
plot_model(RaoQModel, type = "pred", terms = c("trt_type2", "richness_RR_mean"))

RaoQFig <- ggplot(data=meansRaoQModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0, size=3) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('Rao Q\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


summary(MNTDModel <- lme(mntd_RR_mean ~ as.factor(trt_type2) + richness_RR_mean,
                         data=na.omit(subset(allDivRRmean, trt_type2!='herb_removal')),
                         random=~1|site_proj_comm))
anova.lme(MNTDModel, type='sequential') #significant effect of trt
meansMNTDModel <- emmeans(MNTDModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansMNTDModelOutput <- as.data.frame(meansMNTDModel$emmeans)
plot_model(MNTDModel, type = "pred", terms = c("trt_type2", "richness_RR_mean"))

MNTDFig <- ggplot(data=meansMNTDModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0, size=3) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('MNTD\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


summary(MNTDtraitsModel <- lme(MNTD_traits_RR_mean ~ as.factor(trt_type2) + richness_RR_mean,
                         data=na.omit(subset(allDivRRmean, trt_type2!='herb_removal')),
                         random=~1|site_proj_comm))
anova.lme(MNTDtraitsModel, type='sequential') #significant effect of trt
meansMNTDtraitsModel <- emmeans(MNTDtraitsModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansMNTDtraitsModelOutput <- as.data.frame(meansMNTDtraitsModel$emmeans)
plot_model(MNTDtraitsModel, type = "pred", terms = c("trt_type2", "richness_RR_mean"))

MNTDtraitsFig <- ggplot(data=meansMNTDtraitsModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0, size=3) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('Trait MNTD\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


pushViewport(viewport(layout=grid.layout(3,2)))
print(richFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(hillFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(mpdFig, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(RaoQFig, vp=viewport(layout.pos.row=3, layout.pos.col=1))
print(MNTDFig, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(MNTDtraitsFig, vp=viewport(layout.pos.row=3, layout.pos.col=2))
#export at 1000x1800

### see sCoRRE_dCCA_traits script to find PCA of case study examples of extreme responraw for each trt type
# N example of FDis and MNTD decreasing effect size: 	KUFS::E2::0::N1S0H0 or YMN::NitAdd::0::N80 (not CUL::Culardoch::0 N50 or maerc::fireplots::0 unuu)
# irrigation example of increasing FDis and MNTD effect size: KNZ::IRG::l i or SEV::WENNDEx::0 P or MNR::watfer::0::W
# drought example of increasing FDis and MNTD effect size: SFREC::GrazePrecip::G4 D or HAYS::Precip::0 reduction


##### treatment magnitude #####
allDivRRtrt <- allDivRRmean %>% 
  left_join(trt) %>% 
  select(-treatment_year, -calendar_year) %>% 
  unique()

#N additions
nDivRR <- allDivRRtrt %>% 
  filter(trt_type2=='N') %>% 
  mutate(n=ifelse(site_proj_comm=='DL::NSFC::0'&trt_type2=='N', 10, n))

summary(nRaoQModel <- lme(RaoQ_RR_mean ~ n, #tried polynomial, but linear was better
                          data=na.omit(subset(nDivRR)),
                          random=~1|site_proj_comm))
anova.lme(nRaoQModel, type='sequential')
coef(nRaoQModel)
summary.tablefunc(nRaoQModel)

nRaoQFig <- ggplot(data=nDivRR, aes(x=n, y=RaoQ_RR_mean)) +
  geom_point(size=2) +
  # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
  geom_smooth(method='lm', formula=y~x, color='black') +
  geom_hline(yintercept=0) +
  # coord_cartesian(ylim=c(-1,2.2)) +
  ylab('Rao Q\nEffect Size') + xlab(bquote('N added '(gm^-2)))


summary(nFMNTDModel <- lme(MNTD_traits_RR_mean ~ n, #tried polynomial, but linear was better
                          data=na.omit(subset(nDivRR)),
                          random=~1|site_proj_comm))
anova.lme(nFMNTDModel, type='sequential')
coef(nFMNTDModel)
summary.tablefunc(nFMNTDModel) #not significant

nFMNTDFig <- ggplot(data=nDivRR, aes(x=n, y=MNTD_traits_RR_mean)) +
  geom_point(size=2) +
  # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
  # geom_smooth(method='lm', formula=y~x, color='black') +
  geom_hline(yintercept=0) +
  # coord_cartesian(ylim=c(-1,2.2)) +
  ylab('Functional MNTD\nEffect Size') + xlab(bquote('N added '(gm^-2)))


summary(nMPDModel <- lme(mpd_RR_mean ~ n, #tried polynomial, but linear was better
                          data=subset(nDivRR, !is.na(mpd_RR_mean)),
                          random=~1|site_proj_comm))
anova.lme(nMPDModel, type='sequential')
coef(nMPDModel)
summary.tablefunc(nMPDModel)

nMPDFig <- ggplot(data=nDivRR, aes(x=n, y=mpd_RR_mean)) +
  geom_point(size=2)+
  # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
  geom_smooth(method='lm', formula=y~x, color='black') +
  geom_hline(yintercept=0) +
  # coord_cartesian(ylim=c(-1,2.2)) +
  ylab('MPD\nEffect Size') + xlab(bquote('N added '(gm^-2)))

summary(nMNTDModel <- lme(mntd_RR_mean ~ n, #tried polynomial, but linear was better
                         data=subset(nDivRR, !is.na(mntd_RR_mean)),
                         random=~1|site_proj_comm))
anova.lme(nMNTDModel, type='sequential')
coef(nMNTDModel)
summary.tablefunc(nMNTDModel) #not significant

nMNTDFig <- ggplot(data=nDivRR, aes(x=n, y=mntd_RR_mean)) +
  geom_point(size=2)+
  # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
  # geom_smooth(method='lm', formula=y~x, color='black') +
  geom_hline(yintercept=0) +
  # coord_cartesian(ylim=c(-1,2.2)) +
  ylab('MNTD\nEffect Size') + xlab(bquote('N added '(gm^-2)))

#N magnitude figure
pushViewport(viewport(layout=grid.layout(2,2)))
print(nMPDFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(nMNTDFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(nRaoQFig, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(nFMNTDFig, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export at 1000x1000



#precip
precipDivRR <- allDivRRtrt %>% 
  filter(precip!=0)

summary(precipFMNTDModel <- lme(MNTD_traits_RR_mean ~ precip, #tried a polynomial model, but linear was better fit
                          data=na.omit(subset(precipDivRR)),
                          random=~1|site_proj_comm))
anova.lme(precipFMNTDModel, type='sequential')
coef(precipFMNTDModel)
summary.tablefunc(precipFMNTDModel)

precipFMNTDFig <- ggplot(data=precipDivRR, aes(x=precip, y=MNTD_traits_RR_mean)) +
  geom_point(size=2)+
  # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
  geom_smooth(method='lm', formula=y~x, color='black') +
  geom_hline(yintercept=0) +
  # coord_cartesian(ylim=c(-1,2.2)) +
  ylab('Functional MNTD\nEffect Size') + xlab('Precipitation Manipulation (%)')


summary(precipRaoQModel <- lme(RaoQ_RR_mean ~ precip, #tried a polynomial model, but linear was better fit
                               data=na.omit(subset(precipDivRR)),
                               random=~1|site_proj_comm))
anova.lme(precipRaoQModel, type='sequential')
coef(precipRaoQModel)
summary.tablefunc(precipRaoQModel)

precipRaoQFig <- ggplot(data=precipDivRR, aes(x=precip, y=RaoQ_RR_mean)) +
  geom_point(size=2)+
  # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
  geom_smooth(method='lm', formula=y~x, color='black') +
  geom_hline(yintercept=0) +
  # coord_cartesian(ylim=c(-1,2.2)) +
  ylab('Rao Q\nEffect Size') + xlab('Precipitation Manipulation (%)')


summary(precipMNTDModel <- lme(mntd_RR_mean ~ precip,
                          data=precipDivRR,
                          random=~1|site_proj_comm))
anova.lme(precipMNTDModel, type='sequential')
coef(precipMNTDModel)
summary.tablefunc(precipMNTDModel)

precipMNTDFig <- ggplot(data=precipDivRR, aes(x=precip, y=mntd_RR_mean)) +
  geom_point(size=2)+
  # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
  geom_smooth(method='lm', formula=y~poly(x,2), color='black') + 
  geom_hline(yintercept=0) +
  # coord_cartesian(ylim=c(-1,2.2)) +
  ylab('MNTD\nEffect Size') + xlab('Precipitation Manipulation (%)')


summary(precipMPDModel <- lme(mpd_RR_mean ~ precip,
                               data=precipDivRR,
                               random=~1|site_proj_comm))
anova.lme(precipMPDModel, type='sequential')
coef(precipMPDModel)
summary.tablefunc(precipMPDModel)

precipMPDFig <- ggplot(data=precipDivRR, aes(x=precip, y=mpd_RR_mean)) +
  geom_point(size=2)+
  # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
  geom_smooth(method='lm', formula=y~poly(x,2), color='black') + 
  geom_hline(yintercept=0) +
  # coord_cartesian(ylim=c(-1,2.2)) +
  ylab('MPD\nEffect Size') + xlab('Precipitation Manipulation (%)')


#combined N and precip magnitude figure
pushViewport(viewport(layout=grid.layout(2,2)))
print(precipMPDFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(precipMNTDFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(precipRaoQFig, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(precipFMNTDFig, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export at 1000x1000

# #temp
# tempDivRR <- allDivRRtrt %>% 
#   filter(temp!=0)
# 
# summary(tempFDisModel <- lme(FDis_RR_mean ~ temp, #tried a polynomial model, but linear was better fit
#                                data=na.omit(subset(tempDivRR, FDis_RR_mean<5)),
#                                random=~1|site_proj_comm))
# anova.lme(tempFDisModel, type='sequential')
# coef(tempFDisModel)
# summary.tablefunc(tempFDisModel)
# 
# ggplot(data=tempDivRR, aes(x=temp, y=FDis_RR_mean)) +
#   geom_point(size=2)+
#   # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
#   geom_smooth(method='lm', formula=y~x, color='black') +
#   geom_hline(yintercept=0) +
#   coord_cartesian(ylim=c(-0.5,1.25)) +
#   ylab('Functional Dispersion\nEffect Size') + xlab('Temperature Manipulation')
# 
# 
# summary(tempMNTDModel <- lme(mntd_diff_mean ~ temp, #tried a polynomial model, but linear was better fit
#                                data=tempDivRR,
#                                random=~1|site_proj_comm))
# anova.lme(tempMNTDModel, type='sequential')
# coef(tempMNTDModel)
# summary.tablefunc(tempMNTDModel)
# 
# ggplot(data=tempDivRR, aes(x=temp, y=mntd_diff_mean)) +
#   geom_point(size=2)+
#   # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
#   geom_smooth(method='lm', formula=y~x, color='black') +
#   geom_hline(yintercept=0) +
#   coord_cartesian(ylim=c(-0.5,1.25)) +
#   ylab('Phylogenetic Diversity (MNTD)\nEffect Size') + xlab('Temperature Manipulation')