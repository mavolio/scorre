################################################################################
##  sCoRRE_phylogenetic_diversity_analysis.R: Examining the effects of global change experiments on phylogenetic diversity within the CoRRE database.
##
##  Author: Kimberly Komatsu
##  Date created: Sept 22, 2020
################################################################################

library(lme4)
library(lmerTest)
library(ggeffects)
library(grid)
library(performance)
library(tidyverse)

#laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data')


#functions
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

#theme set
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###read in data
trt <- read.csv('CoRRE_raw_abundance_Nov2019.csv')%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id)%>%
  unique()%>%
  left_join(read.csv('CoRRE_project_summary.csv'))%>%
  left_join(read.csv('CoRRE_treatment_summary.csv'))

pDiv <- read.csv('CoRRE_pd_metrics.csv')%>%
  separate(plot_id2, c('site_code', 'project_name', 'community_type', 'treatment_year', 'plot_id'), sep='::')%>%
  mutate(treatment_year=as.integer(treatment_year))%>%
  left_join(trt)

# ggplot(data=subset(pDivMaxYr, trt_type %in% c('control', 'CO2', 'drought', 'irr', 'N', 'mult_nutrient')), aes(x=trt_type, y=pd.ses)) +
#   geom_boxplot()

control <- pDiv%>%
  filter(trt_type=='control')%>%
  rename(pd.ses_ctl=pd.ses, pd.pval_ctl=pd.pval, mpd.ses_ctl=mpd.ses, mpd.pval_ctl=mpd.pval, mntd.ses_ctl=mntd.ses, mntd.pval_ctl=mntd.pval)%>%
  select(site_code, project_name, community_type, treatment_year, pd.ses_ctl, pd.pval_ctl, mpd.ses_ctl, mpd.pval_ctl, mntd.ses_ctl, mntd.pval_ctl)

pDivRR <- pDiv%>%
  filter(trt_type!='control')%>%
  left_join(control)%>%
  mutate(pd_diff=(pd.ses-pd.ses_ctl), mpd_diff=(mpd.ses-mpd.ses_ctl), mntd_diff=(mntd.ses-mntd.ses_ctl))

pDivMaxYr <- pDivRR%>%
  group_by(site_code, project_name, community_type)%>%
  summarise(max_year=max(treatment_year))%>%
  ungroup()%>%
  left_join(pDivRR)%>%
  filter(treatment_year==max_year)



###figures

#by trt type
PDfig <- ggplot(data=subset(pDivMaxYr, trt_type %in% c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient')), aes(x=trt_type, y=pd_diff)) +
  geom_boxplot() +
  xlab('') + ylab('PD.SES difference (trt-ctl)') +
  geom_hline(yintercept=0) +
  scale_x_discrete(limits=c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient'),
                   labels=c('CO2', 'drought', 'irrigation', 'N', 'N*P', 'many nuts.')) +
  ylim(-10,8)
MPDfig <- ggplot(data=subset(pDivMaxYr, trt_type %in% c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient')), aes(x=trt_type, y=mpd_diff)) +
  geom_boxplot() +
  xlab('') + ylab('MPD.SES difference (trt-ctl)') +
  geom_hline(yintercept=0) +
  scale_x_discrete(limits=c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient'),
                   labels=c('CO2', 'drought', 'irrigation', 'N', 'N*P', 'many nuts.')) +
  ylim(-10,8)
MNTDfig <- ggplot(data=subset(pDivMaxYr, trt_type %in% c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient')), aes(x=trt_type, y=mntd_diff)) +
  geom_boxplot() +
  xlab('') + ylab('MNTD.SES difference (trt-ctl)') +
  geom_hline(yintercept=0) +
  scale_x_discrete(limits=c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient'),
                   labels=c('CO2', 'drought', 'irrigation', 'N', 'N*P', 'many nuts.')) +
  ylim(-10,8)

pushViewport(viewport(layout=grid.layout(1,3)))
print(PDfig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(MPDfig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(MNTDfig, vp=viewport(layout.pos.row=1, layout.pos.col=3))
#export at 1800x600

#by site characteristics (for PD only)
lengthfig <- ggplot(data=subset(pDivMaxYr, trt_type %in% c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient')), aes(x=max_year, y=pd_diff)) +
  geom_point() +
  geom_smooth(method='lm', color='grey', size=3) +
  xlab('Experiment Length') + ylab('PD.SES difference (trt-ctl)') +
  ylim(-10,8)
gammafig <- ggplot(data=subset(pDivMaxYr, trt_type %in% c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient')), aes(x=rrich, y=pd_diff)) +
  geom_point() +
  geom_smooth(method='lm', color='grey', size=3) +
  xlab('Gamma Diversity') + ylab('PD.SES difference (trt-ctl)') +
  ylim(-10,8)
anppfig <- ggplot(data=subset(pDivMaxYr, trt_type %in% c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient')), aes(x=anpp, y=pd_diff)) +
  geom_point() +
  geom_smooth(method='lm', color='grey', size=3) +
  xlab('Site Productivity') + ylab('PD.SES difference (trt-ctl)') +
  ylim(-10,8)
MAPfig <- ggplot(data=subset(pDivMaxYr, trt_type %in% c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient')), aes(x=MAP, y=pd_diff)) +
  geom_point() +
  geom_smooth(method='lm', color='grey', size=3) +
  xlab('Mean Annual Precipitation (mm)') + ylab('PD.SES difference (trt-ctl)') +
  ylim(-10,8)
MATfig <- ggplot(data=subset(pDivMaxYr, trt_type %in% c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient')), aes(x=MAT, y=pd_diff)) +
  geom_point() +
  geom_smooth(method='lm', color='grey', size=3) +
  xlab('Mean Annual Temperature (C)') + ylab('PD.SES difference (trt-ctl)') +
  ylim(-10,8)

pushViewport(viewport(layout=grid.layout(2,3)))
print(lengthfig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(gammafig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(anppfig, vp=viewport(layout.pos.row=1, layout.pos.col=3))
print(MAPfig, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(MATfig, vp=viewport(layout.pos.row=2, layout.pos.col=3))
#export at 1800x1200

#by magnitude of treatment
Nfig <- ggplot(data=subset(pDivMaxYr, trt_type=='N'), aes(x=n, y=pd_diff)) +
  geom_point() +
  geom_smooth(method='lm', color='grey', size=3) +
  xlab('N added (g/m2)') + ylab('PD.SES difference (trt-ctl)') +
  ylim(-10,8)
droughtfig <- ggplot(data=subset(pDivMaxYr, trt_type=='drought'), aes(x=precip, y=pd_diff)) +
  geom_point() +
  geom_smooth(method='lm', color='grey', size=3) +
  xlab('Precipitation Reduced (%)') + ylab('PD.SES difference (trt-ctl)') +
  ylim(-10,8)
irrfig <- ggplot(data=subset(pDivMaxYr, trt_type=='irr'), aes(x=precip, y=pd_diff)) +
  geom_point() +
  geom_smooth(method='lm', color='grey', size=3) +
  xlab('Precipitation Increased (%)') + ylab('PD.SES difference (trt-ctl)') +
  ylim(-10,8)

pushViewport(viewport(layout=grid.layout(1,3)))
print(Nfig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(droughtfig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(irrfig, vp=viewport(layout.pos.row=1, layout.pos.col=3))
#export at 1800x600

#by initial phylogenetic diversity at a site
PDfig <- ggplot(data=subset(pDivMaxYr, trt_type %in% c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient')), aes(x=pd.ses_ctl, y=pd_diff)) +
  geom_point() +
  geom_smooth(method='lm', color='grey', size=3) +
  xlab('PD.SES in Control') + ylab('PD.SES difference (trt-ctl)') +
  ylim(-10,8)
MPDfig <- ggplot(data=subset(pDivMaxYr, trt_type %in% c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient')), aes(x=mpd.ses_ctl, y=mpd_diff)) +
  geom_point() +
  geom_smooth(method='lm', color='grey', size=3) +
  xlab('MPD.SES in Control') + ylab('MPD.SES difference (trt-ctl)') +
  ylim(-10,8)
MNTDfig <- ggplot(data=subset(pDivMaxYr, trt_type %in% c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient')), aes(x=mntd.ses_ctl, y=mntd_diff)) +
  geom_point() +
  geom_smooth(method='lm', color='grey', size=3) +
  xlab('MNTD.SES in Control') + ylab('MNTD.SES difference (trt-ctl)') +
  ylim(-10,8)

pushViewport(viewport(layout=grid.layout(1,3)))
print(PDfig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(MPDfig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(MNTDfig, vp=viewport(layout.pos.row=1, layout.pos.col=3))
#export at 1800x600