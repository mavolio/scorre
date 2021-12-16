################################################################################
##  sCoRRE_phylogeneticDiversity_causalModels.R: Examining differences in phylogenetic and functional diversity with causal modeling of the CoRRE database.
##
##  Author: Kimberly Komatsu
##  Date created: December 13, 2021
################################################################################

library(data.table)
library(codyn)
library(fixest)
library(lme4)
library(tidyverse)


setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\paper 2_PD and FD responses\\data\\')  #kim's laptop

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

#theme set
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))



##### data #####
trt <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\RawAbundance.csv')%>%
select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id)%>%
  unique()%>%
  left_join(read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\ExperimentInfo.csv'))%>%
  group_by(site_code, project_name, community_type)%>%
  mutate(experiment_length=max(treatment_year))%>%
  ungroup()


pDiv <- read.csv('CoRRE_pd_metrics_non_weighted.csv')%>%
  separate(identifier, into=c("site_code", "project_name", "community_type", "treatment_year", "plot_id"), sep="::")%>%
  mutate(treatment_year=as.integer(treatment_year))%>%
  left_join(trt)%>%
  full_join(read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteBiotic.csv'))%>%
  full_join(read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteLocationClimate.csv'))%>%
  mutate(site_proj_comm=paste(site_code,  project_name, community_type, sep='_'))

fDiv <- read.csv('CoRRE_functionalDiversity_20211215.csv')

relCover <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\RelativeCover.csv')%>%
  mutate(replicate=paste(site_code, project_name, community_type, plot_id, sep='::'))

rDiv <- community_structure(relCover, time.var="treatment_year", abundance.var="relcov", replicate.var="replicate")%>%
  separate(replicate, into=c("site_code", "project_name", "community_type", "plot_id"), sep='::')

allDiv <- pDiv%>%
  full_join(fDiv)%>%
  full_join(rDiv)%>%
  filter(site_code!='DCGS')%>%
  filter(treatment_year>0)

# control <- pDiv%>%
#   filter(trt_type=='control')%>%
#   rename(pd.ses_ctl=pd.ses, mpd.ses_ctl=mpd.ses, mntd.ses_ctl=mntd.ses)%>%
#   select(site_code, project_name, community_type, treatment_year, pd.ses_ctl, mpd.ses_ctl, mntd.ses_ctl)
# 
# pDivRR <- pDiv%>%
#   filter(trt_type!='control')%>%
#   left_join(control)%>%
#   mutate(pd_diff=(pd.ses-pd.ses_ctl), mpd_diff=(mpd.ses-mpd.ses_ctl), mntd_diff=(mntd.ses-mntd.ses_ctl))
# 
# pDivAvg <- pDivRR%>%
#   group_by(site_code, project_name, community_type, treatment_year, treatment, MAP, MAT, rrich, experiment_length, anpp, trt_type, n, p, k, CO2, precip, temp, herb_removal)%>%
#   summarise(pd_diff_avg=mean(pd_diff), mpd_diff_avg=mean(mpd_diff), mntd_diff_avg=mean(mntd_diff))%>%
#   ungroup()

trt_analysis<-trt%>%
  mutate(alltrts=ifelse(trt_type %in% c("control", "CO2","CO2*temp", "mow_clip","burn","burn*graze","disturbance","burn*mow_clip","drought","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N","mult_nutrient","N*P","P","N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp","temp*mow_clip","drought*temp","irr*temp","irr"),1,0))%>%
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
  mutate(trt_type2=ifelse(dist==1, 'disturbance', ifelse(multtrts==1, 'multiple trts', trt_type)))%>%
  select(site_code, project_name, community_type, treatment, alltrts, dist, tCO2, drought, therb_removal, irg, ttemp, tn, tp, multtrts, trt_type2)%>%
  unique()


#pick treatments here.
allDivTrt <- allDiv%>%
  right_join(trt_analysis)%>%
  mutate(trt_binary=ifelse(plot_mani>0, 1, 0))%>%
  select(-FDiv)

# write.csv(allDivTrt, 'CoRRE_allDiversityMetrics_phyFunAnalysis.csv', row.names=F)


##### figures over time by site #####
allDivYearMeans <- allDivTrt %>%
  group_by(site_proj_comm, site_code, project_name, community_type, treatment_year, treatment, trt_type2)%>%
  summarize_at(vars(mpd.ses, mntd.ses, FDis, richness), list(mean=mean, se=se), na.rm=T)%>%
  ungroup()

site_proj_comm_vector <- unique(allDivTrt$site_proj_comm)

#MPD over time by project
for(PROJ in 1:length(site_proj_comm_vector)){
  ggplot(data=filter(allDivYearMeans, site_proj_comm == site_proj_comm_vector[PROJ]),
         aes(x=treatment_year, y=mpd.ses_mean, col=treatment,
             ymin=mpd.ses_mean-mpd.ses_se, ymax=mpd.ses_mean+mpd.ses_se)) +
    geom_errorbar(width=0.1) +
    geom_point() +
    geom_path() +
    ggtitle(site_proj_comm_vector[PROJ]) +
    theme_bw()
  ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\pd figs\\",
                         site_proj_comm_vector[PROJ], "_mpd.png"))
}

#MNTD over time by project
for(PROJ in 1:length(site_proj_comm_vector)){
  ggplot(data=filter(allDivYearMeans, site_proj_comm == site_proj_comm_vector[PROJ]),
         aes(x=treatment_year, y=mntd.ses_mean, col=treatment,
             ymin=mntd.ses_mean-mntd.ses_se, ymax=mntd.ses_mean+mntd.ses_se)) +
    geom_errorbar(width=0.1) +
    geom_point() +
    geom_path() +
    ggtitle(site_proj_comm_vector[PROJ]) +
    theme_bw()
  ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\pd figs\\",
                         site_proj_comm_vector[PROJ], "_mntd.png"))
}


##### figures by plot (avg over years) #####
allDivPlotMeans <- allDivTrt %>%
  group_by(site_proj_comm, site_code, project_name, community_type, plot_id, treatment, trt_type2)%>%
  summarize_at(vars(mpd.ses, mntd.ses, FDis, richness), list(mean=mean, se=se), na.rm=T)%>%
  ungroup()


# #MPD vs trait functional diversity
# for(PROJ in 1:length(site_proj_comm_vector)){
#   ggplot(data=filter(allDivPlotMeans, site_proj_comm == site_proj_comm_vector[PROJ]),
#          aes(x=mpd.ses_mean, y=FDis_mean, col=treatment)) +
#     geom_point() +
#     geom_smooth(method='lm', se=F) +
#     ggtitle(site_proj_comm_vector[PROJ]) +
#     theme_bw()
#   ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\pd figs\\",
#                          site_proj_comm_vector[PROJ], "_mpd_Fdis.png"))
# }
# 
# #MNTD vs trait functional diversity
# for(PROJ in 1:length(site_proj_comm_vector)){
#   ggplot(data=filter(allDivPlotMeans, site_proj_comm == site_proj_comm_vector[PROJ]),
#          aes(x=mntd.ses_mean, y=FDis_mean, col=treatment)) +
#     geom_point() +
#     geom_smooth(method='lm', se=F) +
#     ggtitle(site_proj_comm_vector[PROJ]) +
#     theme_bw()
#   ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\pd figs\\",
#                          site_proj_comm_vector[PROJ], "_mntd_Fdis.png"))
# }
# 
# #richness vs trait functional diversity
# for(PROJ in 1:length(site_proj_comm_vector)){
#   ggplot(data=filter(allDivPlotMeans, site_proj_comm == site_proj_comm_vector[PROJ]),
#          aes(x=richness_mean, y=FDis_mean, col=treatment)) +
#     geom_point() +
#     geom_smooth(method='lm', se=F) +
#     ggtitle(site_proj_comm_vector[PROJ]) +
#     theme_bw()
#   ggsave(filename=paste0("C:\\Users\\lapie\\Desktop\\pd figs\\",
#                          site_proj_comm_vector[PROJ], "_richness_Fdis.png"))
# }


##### global trends - comparing diversity metrics #####
allDivGlobal <- allDivPlotMeans %>%
  group_by(site_proj_comm, site_code, project_name, community_type, treatment, trt_type2)%>%
  summarize_at(vars(mpd.ses_mean, mntd.ses_mean, FDis_mean, richness_mean), list(mean=mean, se=se), na.rm=T)%>%
  ungroup()

# ggplot(data=allDivGlobal, aes(x=mpd.ses_mean_mean, y=FDis_mean_mean)) +
#   geom_point() +
#   geom_smooth(method='lm', se=F, formula = y ~ x + I(x^2)) +
#   facet_wrap(~trt_type2)
# 
# ggplot(data=allDivGlobal, aes(x=mntd.ses_mean_mean, y=FDis_mean_mean)) +
#   geom_point() +
#   geom_smooth(method='lm', se=F, formula = y ~ x + I(x^2)) +
#   facet_wrap(~trt_type2)
# 
# ggplot(data=allDivGlobal, aes(x=richness_mean_mean, y=FDis_mean_mean)) +
#   geom_point() +
#   geom_smooth(method='lm', se=F, formula = y ~ x + I(x^2)) +
#   facet_wrap(~trt_type2)


### TO CONSIDER: should we look at each treatment within project and see if the slopes differ from 1 (more or less Fdiv per change in phy div), and does that differ by which treatment it is?
#can we identify trajectories through time for each plot?

#plot trajectories through time, this is not very straight forward
ggplot(data=subset(allDivTrt, project_name=='Fert1'&treatment %in% c('full_nut', 'control')), aes(x=mntd.ses, y=FDis, color=treatment, label=as.numeric(treatment_year))) +
  geom_text() +
  geom_path() +
  facet_wrap(~plot_id)




##### causal modeling #####
allDivTrtFixed <- allDivTrt%>%
  mutate(site_year=paste(site_code, project_name, community_type, treatment_year, sep="::"))%>% #creates dummy site-by-year variable (time variant site effects)
  mutate(site_trt=paste(site_code, project_name, community_type, trt_type2, sep='::'))%>% #creates dummy site-by-plot variable (time invariant plot effects)
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::')) #creates dummy site-experiment variable

#individual treatments
mpdDivCausalModel <- feols(mpd.ses ~ n + p + CO2 + drought + irg + temp + dist | MAP + MAT + rrich + anpp, data=allDivTrtFixed, panel.id=~site_project_comm + treatment_year)
etable(mpdDivCausalModel, cluster="site_year") #trt effects
fixedEffects = fixef(mpdDivCausalModel) #env effects
plot(fixedEffects)

mntdDivCausalModel <- feols(mntd.ses ~ n + p + CO2 + drought + irg + temp + dist | MAP + MAT + rrich + anpp, data=allDivTrtFixed, panel.id=~site_project_comm + treatment_year)
etable(mntdDivCausalModel, cluster="site_year") #trt effects
fixedEffects = fixef(mntdDivCausalModel) #env effects
plot(fixedEffects)

fdisDivCausalModel <- feols(FDis ~ n + p + CO2 + drought + irg + temp + dist | MAP + MAT + rrich + anpp, data=allDivTrtFixed, panel.id=~site_project_comm + treatment_year)
etable(fdisDivCausalModel, cluster="site_year") #trt effects
fixedEffects = fixef(fdisDivCausalModel) #env effects
plot(fixedEffects)

# #all together
# pDivCausalModel <- feols(mpd_diff_avg ~ trt_type2 | site_year, data=pDivAvgTrtFixed)
# etable(pDivCausalModel,
#        cluster="site_year")

##### exploratory figures #####
control <- allDivTrt%>%
  filter(trt_type=='control')%>%
  rename(mpd.ses_ctl=mpd.ses, mntd.ses_ctl=mntd.ses, FDis_ctl=FDis, richness_ctl=richness)%>%
  group_by(site_code, project_name, community_type, treatment_year)%>%
  summarize_at(vars(mpd.ses_ctl, mntd.ses_ctl, FDis_ctl, richness_ctl), list(mean=mean), na.rm=T)%>%
  ungroup()

allDivRR <- allDivTrt%>%
  filter(trt_type!='control')%>%
  left_join(control)%>%
  mutate(mpd_diff=(mpd.ses-mpd.ses_ctl_mean), mntd_diff=(mntd.ses-mntd.ses_ctl_mean), FDis_RR=((FDis-FDis_ctl_mean)/FDis_ctl_mean), richness_RR=((richness-richness_ctl_mean)/richness_ctl_mean))%>%
  gather(key=env_variable, value=env_value, c("MAP", "MAT", "rrich", "anpp"))%>%
  group_by(site_code, project_name, community_type, treatment, trt_type2, env_variable)%>%
  summarise_at(vars(mpd_diff, mntd_diff, FDis_RR, richness_RR, env_value), list(mean=mean), na.rm=T)%>%
  ungroup%>%
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  full_join(read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteBiotic.csv'))%>%
  full_join(read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteLocationClimate.csv'))

ggplot(data=allDivRR, aes(x=env_value_mean, y=FDis_RR_mean)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  geom_hline(yintercept=0) +
  facet_grid(trt_type2~env_variable, scales='free')
  
  
test <- subset(allDivRR, project_name=='pplots'&treatment_year==10&treatment=='N2P3')


##### repeated measures model #####
library(nlme)
library(emmeans)
library(performance)

options(contrasts=c('contr.sum','contr.poly')) 

summary(FDisModel <- lme(FDis_RR_mean ~ as.factor(trt_type2),
                         data=na.omit(allDivRR),
                         random=~1|site_proj_comm))
anova.lme(FDisModel, type='sequential')
emmeans(FDisModel, pairwise~as.factor(trt_type2), adjust="tukey")

summary(glm(FDis ~ trt_binary*trt_type2*treatment_year, data=allDivTrt))
