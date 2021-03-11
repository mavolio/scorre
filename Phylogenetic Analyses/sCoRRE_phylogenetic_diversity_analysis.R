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
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared')
setwd('C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared')


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
trt <- read.csv('CoRRE data\\CoRRE data\\community composition\\CoRRE_RawAbundanceMar2021.csv')%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id)%>%
  unique()%>%
  left_join(read.csv('CoRRE data\\CoRRE data\\CoRRE_project_summary.csv'))%>%
  left_join(read.csv('CoRRE data\\CoRRE data\\basic dataset info\\ExperimentInfo.csv'))

pDiv <- read.csv('paper 2_PD and FD responses\\data\\CoRRE_pd_metrics_weighted_kjk.csv')%>%
  separate(plot_id2, c('site_code', 'project_name', 'community_type', 'treatment_year', 'plot_id'), sep='::')%>%
  mutate(treatment_year=as.integer(treatment_year))%>%
  left_join(trt)

# ggplot(data=subset(pDivMaxYr, trt_type %in% c('control', 'CO2', 'drought', 'irr', 'N', 'mult_nutrient')), aes(x=trt_type, y=pd.ses)) +
#   geom_boxplot()

control <- pDiv%>%
  filter(trt_type=='control')%>%
  rename(pd.ses_ctl=pd.ses, mpd.ses_ctl=mpd.ses, mntd.ses_ctl=mntd.ses)%>%
  select(site_code, project_name, community_type, treatment_year, pd.ses_ctl, mpd.ses_ctl, mntd.ses_ctl)

pDivRR <- pDiv%>%
  filter(trt_type!='control')%>%
  left_join(control)%>%
  mutate(pd_diff=(pd.ses-pd.ses_ctl), mpd_diff=(mpd.ses-mpd.ses_ctl), mntd_diff=(mntd.ses-mntd.ses_ctl))

pDivAvg <- pDivRR%>%
  group_by(site_code, project_name, community_type, treatment, MAP, MAT, rrich, experiment_length, anpp, trt_type, n, p, k, CO2, precip, temp, herb_removal)%>%
  summarise(pd_diff_avg=mean(pd_diff), mpd_diff_avg=mean(mpd_diff), mntd_diff_avg=mean(mntd_diff))%>%
  ungroup()


trt_analysis<-trt%>%
  mutate(alltrts=ifelse(trt_type %in% c("CO2","CO2*temp", "mow_clip","burn","burn*graze","disturbance","burn*mow_clip","drought","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N","mult_nutrient","N*P","P","N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp","temp*mow_clip","drought*temp","irr*temp","irr"),1,0))%>%
  filter(alltrts==1)%>%
  mutate(dist=ifelse(trt_type %in% c("mow_clip","burn","burn*graze","disturbance","burn*mow_clip"), 1, 0),
         tCO2=ifelse(trt_type %in% c("CO2"), 1, 0),
         drought=ifelse(trt_type %in% c("drought"), 1, 0),
         therb_removal=ifelse(trt_type %in% c("herb_removal"), 1, 0),
         irg=ifelse(trt_type %in% c("irr"), 1, 0),
         ttemp=ifelse(trt_type %in% c("temp"), 1, 0),
         tn=ifelse(trt_type %in% c("N"), 1, 0),
         tp=ifelse(trt_type %in% c("P"), 1, 0),
         multtrts=ifelse(trt_type %in% c("CO2*temp", "burn*graze","burn*mow_clip","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp*mow_clip","drought*temp","irr*temp","mult_nutrient","N*P"),1,0))%>%
  select(site_code, project_name, community_type, treatment, alltrts, dist, tCO2, drought, therb_removal, irg, ttemp, tn, tp, multtrts)%>%
  unique()


#pick treatments here.
pDivAvgTrt <- pDivAvg%>%
  right_join(trt_analysis)

# write.csv(pDivAvgTrt, 'paper 2_PD and FD responses\\data\\PD_mean_effectsizes.csv')

subco2<-pDivAvgTrt%>%
  filter(tCO2==1)%>%
  mutate(trt_type2="co2")
subdrought<-pDivAvgTrt%>%
  filter(drought==1)%>%
  mutate(trt_type2="drought")
subirg<-pDivAvgTrt%>%
  filter(irg==1)%>%
  mutate(trt_type2="irg")
subtemp<-pDivAvgTrt%>%
  filter(ttemp==1)%>%
  mutate(trt_type2="temp")
subn<-pDivAvgTrt%>%
  filter(tn==1)%>%
  mutate(trt_type2="n")
subp<-pDivAvgTrt%>%
  filter(tp==1)%>%
  mutate(trt_type2="p")
subdist<-pDivAvgTrt%>%
  filter(dist==1)%>%
  mutate(trt_type2="dist")
subherbr<-pDivAvgTrt%>%
  filter(therb_removal==1)%>%
  mutate(trt_type2="herbr")
submult<-pDivAvgTrt%>%
  filter(multtrts==1)%>%
  mutate(trt_type2="mult_trts")

alldat<-submult%>%
  bind_rows(subco2, subtemp, subdist, subdrought, subirg, subn, subp, subherbr)

###figures

#by trt type
ggplot(data=alldat, aes(x=trt_type2, y=mpd_diff_avg)) +
  geom_boxplot() +
  xlab('') + ylab('MPD.SES difference (trt-ctl)') +
  geom_hline(yintercept=0) +
  scale_x_discrete(limits=c('n', 'p', 'co2', "drought", 'irg', 'dist', 'herbr', 'mult_trts'),
                   labels=c('N', 'P', 'CO2', "Drought", 'Irrigation', 'Disturbance', 'Herb_removal','mult_trts'))+
  theme(axis.text.x = element_text(angle = 90))



MPDfig <- ggplot(data=subset(pDivAvgTrt, trt_type %in% c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient')), aes(x=trt_type, y=mpd_diff)) +
  geom_boxplot() +
  xlab('') + ylab('MPD.SES difference (trt-ctl)') +
  geom_hline(yintercept=0) +
  scale_x_discrete(limits=c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient'),
                   labels=c('CO2', 'drought', 'irrigation', 'N', 'N*P', 'many nuts.')) +
  ylim(-10,8)
MNTDfig <- ggplot(data=subset(pDivAvgTrt, trt_type %in% c('CO2', 'drought', 'irr', 'N', 'N*P', 'mult_nutrient')), aes(x=trt_type, y=mntd_diff)) +
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