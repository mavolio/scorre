################################################################################
##  sCoRRE_communityResponses_mixedModels.R: Examining differences in
##  phylogenetic and functional diversity within the CoRRE database.
##
##  Author: Kimberly Komatsu
##  Date created: December 13, 2021
################################################################################

rm(list = ls())

library(PerformanceAnalytics)
library(nlme)
library(emmeans)
library(sjPlot)
library(performance)
library(grid)
library(data.table)
library(fixest)
library(car)
library(interactions)
library(tidyverse)



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
trt <- read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\RawAbundanceMarch2024.csv') %>%
select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, block, plot_id) %>%
  unique() %>%
  mutate(plot_id=ifelse(project_name=='IRG', paste(block, plot_id, sep='_'), plot_id)) %>% 
  select(-block) %>% 
  left_join(readRDS('PhyDiv_FuncDiv/trt_info.rds')) %>%
  group_by(site_code, project_name, community_type) %>%
  mutate(experiment_length=max(treatment_year)) %>%
  ungroup() %>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, trt_type2, experiment_length, plot_mani)

#phylogenetic diversity data
pDiv <- readRDS('PhyDiv_FuncDiv/phylogeneticDiversityMetrics.rds') %>%
  separate(plot_id2, into=c("site_code", "project_name", "community_type", "drop", "plot_id"), sep="::") %>%
  select(-drop, -site_proj_comm)

#functional diversity data
fDiv <- readRDS('PhyDiv_FuncDiv/functionalDiversityMetrics.rds') %>% 
  select(site_code, project_name, community_type, 
         calendar_year, plot_id, Q, RaoQ_ses) %>% 
  rename(RaoQ=Q)

#taxonomic diversity data
rDiv <- readRDS('PhyDiv_FuncDiv/rDiv.rds') %>% 
  select(-treatment)
  

#merge all data on diversity metrics (phylogenetic, functional, species), experimental treatments, and site characteristics
allDivTrt <- pDiv %>% #phylogenetic metrics
  full_join(fDiv) %>% #functional metrics
  full_join(rDiv) %>% #species metrics
  left_join(trt) %>% #treatments
  filter(treatment_year>0) %>% 
  full_join(read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteBiotic.csv')) %>% #site anpp and regional richness
  full_join(read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteLocationClimate.csv')) %>% #site MAP and MAT
  mutate(site_proj_comm=paste(site_code,  project_name, community_type, sep='::')) %>%
  mutate(site_proj_comm_trt=paste(site_proj_comm, treatment, sep='::')) %>% 
  na.omit() 

hist(allDivTrt$MPD_phylo)
hist(allDivTrt$MPD_phylo_ses)
# hist(allDivTrt$FDis)
hist(allDivTrt$RaoQ)
hist(allDivTrt$RaoQ_ses)
hist(allDivTrt$richness)
hist(allDivTrt$hill)
hist(allDivTrt$Evar)
  
# write.csv(allDivTrt, 'paper 2_PD and FD responses\\data\\CoRRE_allDiversityMetrics_phyFunAnalysis.csv', row.names=F)


##### calculating response ratios #####
#filter control plots
control <- allDivTrt %>%
  filter(trt_type2=='control') %>%
  rename(mpd.raw_ctl=MPD_phylo,
         mpd.ses_ctl=MPD_phylo_ses,
         # mntd.raw_ctl=mntd.raw,
         # mntd.ses_ctl=mntd.ses,
         RaoQ_raw_ctl=RaoQ,
         RaoQ_ses_ctl=RaoQ_ses,
         # MNTD_traits_raw_ctl=MNTD_traits_raw,
         # MNTD_traits_ses_ctl=MNTD_traits_ses,
         # FDis_ctl=FDis,
         # RaoQ_ctl=RaoQ,
         richness_ctl=richness,
         hill_ctl=hill) %>%
  group_by(site_code, project_name, community_type, treatment_year) %>%
  summarize_at(vars(mpd.raw_ctl,
                    mpd.ses_ctl,
                    # mntd.raw_ctl,
                    # mntd.ses_ctl,
                    RaoQ_raw_ctl,
                    RaoQ_ses_ctl,
                    # MNTD_traits_raw_ctl,
                    # MNTD_traits_ses_ctl,
                    # FDis_ctl,
                    # RaoQ_ctl,
                    richness_ctl,
                    hill_ctl),
               list(mean=mean), na.rm=T) %>% #average across plots and years
  ungroup()

#merge on site characteristics
controlEnv <- control %>%
  group_by(site_code, project_name, community_type) %>%
  summarize_at(vars(mpd.raw_ctl_mean,
                    mpd.ses_ctl_mean,
                    # mntd.raw_ctl_mean,
                    # mntd.ses_ctl_mean,
                    RaoQ_raw_ctl_mean,
                    RaoQ_ses_ctl_mean,
                    # MNTD_traits_raw_ctl_mean,
                    # MNTD_traits_ses_ctl_mean,
                    # FDis_ctl_mean,
                    # RaoQ_ctl_mean,
                    richness_ctl_mean,
                    hill_ctl_mean),
               list(mean=mean), na.rm=T) %>% #average across plots and years
  ungroup() %>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteBiotic.csv')) %>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteLocationClimate.csv')) %>%
  gather(key='env_variable', value='env_value', rrich, anpp, MAP, MAT, aridityValues)


allDivRR <- allDivTrt %>%
  filter(trt_type2!='control') %>%
  left_join(control) %>%
  filter(!is.na(RaoQ_raw_ctl_mean)) %>%  #remove lines where there was no control to compare to due to lack of spp cover for traits; lose 3586 data points __ CHECK THIS!!
  mutate(P_mpd_raw_RR=log(MPD_phylo/mpd.raw_ctl_mean),
         P_mpd_ses_RR=(MPD_phylo-mpd.ses_ctl_mean/mpd.ses_ctl_mean), #percent difference for ses due to neg values
         # P_mntd_raw_RR=log(mntd.raw/mntd.raw_ctl_mean),
         # P_mntd_ses_RR=(mntd.ses-mntd.ses_ctl_mean/mntd.ses_ctl_mean), #percent difference for ses due to neg values
         F_mpd_raw_RR=log(RaoQ/RaoQ_raw_ctl_mean),
         F_mpd_ses_RR=(RaoQ_ses-RaoQ_ses_ctl_mean/RaoQ_ses_ctl_mean), #percent difference for ses due to neg values
         # F_mntd_raw_RR=log(MNTD_traits_raw/MNTD_traits_raw_ctl_mean),
         # F_mntd_ses_RR=(MNTD_traits_ses-MNTD_traits_ses_ctl_mean/MNTD_traits_ses_ctl_mean), #percent difference for ses due to neg values
         # FDis_RR=log(FDis/FDis_ctl_mean),
         # RaoQ_RR=log(RaoQ/RaoQ_ctl_mean),
         richness_RR=log(richness/richness_ctl_mean),
         hill_RR=log(hill/hill_ctl_mean)) %>%
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::')) %>%
  select(site_proj_comm, site_code, project_name, community_type, treatment_year, calendar_year, treatment, trt_type2, plot_mani, plot_id, rrich, anpp, MAP, MAT, aridityValues, P_mpd_raw_RR:hill_RR) #%>%
  # filter(P_mpd_raw_RR<1.9) #filter outliers, removes 3 points

allDivRRmean <- allDivRR %>%
  group_by(site_proj_comm, site_code, project_name, community_type, treatment, trt_type2, plot_id) %>%
  summarise_at(vars(P_mpd_raw_RR,
                    P_mpd_ses_RR,
                    # P_mntd_raw_RR,
                    # P_mntd_ses_RR,
                    F_mpd_raw_RR,
                    F_mpd_ses_RR,
                    # F_mntd_raw_RR,
                    # F_mntd_ses_RR,
                    # FDis_RR,
                    # RaoQ_RR,
                    richness_RR,
                    hill_RR),
               list(mean=mean), na.rm=T) %>%
  ungroup()

hist(allDivRR$P_mpd_raw_RR)
hist(allDivRR$P_mpd_ses_RR)
# hist(allDivRR$P_mntd_raw_RR)
# hist(allDivRR$P_mntd_ses_RR)
hist(allDivRR$F_mpd_raw_RR)
hist(allDivRR$F_mpd_ses_RR)
# hist(allDivRR$F_mntd_raw_RR)
# hist(allDivRR$F_mntd_ses_RR)
# hist(allDivRR$FDis_RR)
# hist(allDivRR$RaoQ_RR)
hist(allDivRR$richness_RR)
hist(allDivRR$hill_RR)

##### check if richness difference is correlated with other diversity metrics #####

chart.Correlation(allDivTrt[c(6:12)]) #raw values

chart.Correlation(allDivRR[16:21]) #response ratios

chart.Correlation(allDivRRmean[8:18]) #mean response ratios over all years


##### mixed effects model #####
#NOTE: these models do not account for biotic or abiotic env drivers at a site or for trt magnitude (but do include a random effect of site)

options(contrasts=c('contr.sum','contr.poly')) 

summary(richModel <- lme(richness_RR_mean ~ as.factor(trt_type2),
                         data=na.omit(subset(allDivRRmean)),
                         random=~1|site_proj_comm))
anova.lme(richModel, type='sequential') #significant effect of trt
meansRichModel <- emmeans(richModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansRichModelOutput <- as.data.frame(meansRichModel$emmeans)

ggplot(data=meansRichModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0, size=3) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('Species Richness') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'herb_removal', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), 
                   breaks=c('multiple trts', 'herb_removal', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), 
                   labels=c('Multiple Trts', 'Herbivore Rem.', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'Phosphorus', 'Nitrogen')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'orange', 'blue', 'orange', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


summary(HillModel <- lme(hill_RR_mean ~ as.factor(trt_type2),
                         data=na.omit(subset(allDivRRmean)),
                         random=~1|site_proj_comm))
anova.lme(HillModel, type='sequential') #significant trt effect
meansHillModel <- emmeans(HillModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansHillModelOutput <- as.data.frame(meansHillModel$emmeans)

hillFig <- ggplot(data=meansHillModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0, size=3) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('Hill Number\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) +
  scale_x_discrete(limits=c('multiple trts', 'herb_removal', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), 
                   breaks=c('multiple trts', 'herb_removal', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), 
                   labels=c('Multiple Trts', 'Herbivore Rem.', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'Phosphorus', 'Nitrogen')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'orange', 'blue', 'orange', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


summary(PmpdModel <- lme(P_mpd_raw_RR_mean ~ as.factor(trt_type2),
                         data=allDivRRmean,
                         random=~1|site_code/richness_RR_mean))
anova.lme(PmpdModel, type='sequential') #sig diff among trts
meansPMPDModel <- emmeans(PmpdModel, ~as.factor(trt_type2), adjust="tukey")
meansPMPDModelOutput <- as.data.frame(meansPMPDModel)

mpdFig <- ggplot(data=meansPMPDModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0, size=3) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('Phylogenetic MPD\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'herb_removal', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), 
                   breaks=c('multiple trts', 'herb_removal', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), 
                   labels=c('Multiple Trts', 'Herbivore Rem.', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'Phosphorus', 'Nitrogen')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'orange', 'blue', 'orange', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


summary(RaoQModel <- lme(F_mpd_raw_RR_mean ~ as.factor(trt_type2),
                         data=allDivRRmean,
                         random=~1|site_code/richness_RR_mean))
anova.lme(RaoQModel, type='sequential') #significant effect of trt
meansRaoQModel <- emmeans(RaoQModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansRaoQModelOutput <- as.data.frame(meansRaoQModel$emmeans)

RaoQFig <- ggplot(data=meansRaoQModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0, size=3) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('Functional Rao Q\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'herb_removal', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), 
                   breaks=c('multiple trts', 'herb_removal', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), 
                   labels=c('Multiple Trts', 'Herbivore Rem.', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'Phosphorus', 'Nitrogen')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'orange', 'blue', 'orange', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


pushViewport(viewport(layout=grid.layout(1,3)))
print(hillFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(mpdFig, vp=viewport(layout.pos.row=1, layout.pos.col=3))
print(RaoQFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
#export at 2000x1000








# summary(FDisModel <- lme(FDis_RR_mean ~ as.factor(trt_type2),
#                          data=allDivRRmean,
#                          random=~1|site_code/richness_RR_mean))
# anova.lme(FDisModel, type='sequential') #significant effect of trt
# meansFDisModel <- emmeans(FDisModel, pairwise~as.factor(trt_type2), adjust="tukey")
# meansFDisModelOutput <- as.data.frame(meansFDisModel$emmeans)
# 
# FDisFig <- ggplot(data=meansFDisModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
#   geom_point(size=5) +
#   geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
#   geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0, size=3) +
#   geom_hline(yintercept=0) +
#   coord_flip() +
#   ylab('Rao Q\nEffect Size') + xlab('') +
#   scale_x_discrete(limits=c('multiple trts', 'herb_removal', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), 
#                    breaks=c('multiple trts', 'herb_removal', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), 
#                    labels=c('Multiple Trts', 'Herbivore Rem.', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'Phosphorus', 'Nitrogen')) + 
#   scale_color_manual(values=c('blue', 'orange', 'orange', 'orange', 'blue', 'darkgrey', 'blue', 'blue', 'orange')) +
#   theme(legend.position='none')






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