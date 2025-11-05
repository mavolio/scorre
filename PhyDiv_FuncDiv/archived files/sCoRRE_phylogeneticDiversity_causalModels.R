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


setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\paper 2_PD and FD responses\\data\\')  #kim's laptop

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
trt <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\RawAbundance.csv')%>%
select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id)%>%
  unique()%>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\ExperimentInfo.csv'))%>%
  group_by(site_code, project_name, community_type)%>%
  mutate(experiment_length=max(treatment_year))%>%
  ungroup()%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, trt_type, experiment_length, plot_mani, n, p, CO2, precip, temp)


pDiv <- read.csv('CoRRE_pd_metrics_non_weighted.csv')%>%
  separate(identifier, into=c("site_code", "project_name", "community_type", "treatment_year", "plot_id"), sep="::")%>%
  mutate(treatment_year=as.integer(treatment_year))

fDiv <- read.csv('CoRRE_functionalDiversity_2022-07-14.csv')

relCover <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeCover_Dec2021.csv')%>%
  mutate(replicate=paste(site_code, project_name, community_type, plot_id, sep='::'))

rDiv <- community_structure(relCover, time.var="treatment_year", abundance.var="relcov", replicate.var="replicate")%>%
  separate(replicate, into=c("site_code", "project_name", "community_type", "plot_id"), sep='::')

allDiv <- pDiv%>%
  full_join(fDiv)%>%
  full_join(rDiv)%>%
  filter(site_code!='DCGS')%>%
  filter(treatment_year>0)%>%
  mutate(pd.raw=ifelse(richness==1, 0, pd.raw),
         pd.ses=ifelse(richness==1, 0, pd.ses),
         mpd.raw=ifelse(richness==1, 0, mpd.raw),
         mpd.ses=ifelse(richness==1, 0, mpd.ses),
         mntd.raw=ifelse(richness==1, 0, mntd.raw),
         mntd.ses=ifelse(richness==1, 0, pd.raw),
         FDis=ifelse(richness==1, 0, FDis),
         RaoQ=ifelse(richness==1, 0, RaoQ))%>%
  left_join(trt)%>%
  full_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteBiotic.csv'))%>%
  full_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteLocationClimate.csv'))%>%
  mutate(site_proj_comm=paste(site_code,  project_name, community_type, sep='_'))%>%
  select(site_code, project_name, community_type, site_proj_comm, treatment_year, calendar_year, treatment, plot_id, mpd.ses, mntd.ses, FDis, RaoQ, richness, Evar, trt_type, experiment_length, plot_mani, rrich, anpp, MAP, MAT, n, p, CO2, precip, temp)%>%
  filter(!(site_proj_comm %in% c('DL_NSFC_0', 'Naiman_Nprecip_0', 'DCGS_gap_0'))) #remove problem expts until they are fixed

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
         # tCO2=ifelse(trt_type %in% c("CO2"), 1, 0),
         drought=ifelse(trt_type %in% c("drought"), 1, 0),
         # therb_removal=ifelse(trt_type %in% c("herb_removal"), 1, 0),
         irg=ifelse(trt_type %in% c("irr"), 1, 0),
         # ttemp=ifelse(trt_type %in% c("temp"), 1, 0),
         # tn=ifelse(trt_type %in% c("N"), 1, 0),
         # tp=ifelse(trt_type %in% c("P"), 1, 0),
         multtrts=ifelse(trt_type %in% c("CO2*temp", "burn*graze","burn*mow_clip","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp*mow_clip","drought*temp","irr*temp","mult_nutrient","N*P"),1,0))%>%
  mutate(trt_type2=ifelse(dist==1, 'disturbance', ifelse(multtrts==1, 'multiple trts', trt_type)))%>%
  select(site_code, project_name, community_type, treatment, alltrts, dist, drought, irg, multtrts, trt_type2)%>%
  unique()


#pick treatments here.
allDivTrt <- allDiv%>%
  right_join(trt_analysis)%>%
  mutate(trt_binary=ifelse(plot_mani>0, 1, 0))%>%
  filter(!(site_proj_comm %in% c('DL_NSFC_0', 'Naiman_Nprecip_0', 'DCGS_gap_0')))%>% #remove problem expts until they are fixed
  mutate_at(c('MAP', 'MAT', 'rrich', 'anpp'), funs(c(scale(.))))%>% #scale site characteristics
  filter(!(is.na(site_proj_comm)))



# write.csv(allDivTrt, 'CoRRE_allDiversityMetrics_phyFunAnalysis.csv', row.names=F)


##### figures over time by site #####
allDivYearMeans <- allDivTrt %>%
  group_by(site_proj_comm, site_code, project_name, community_type, treatment_year, treatment, trt_type2)%>%
  summarize_at(vars(mpd.ses, mntd.ses, FDis, richness), list(mean=mean, se=se), na.rm=T)%>%
  ungroup()

site_proj_comm_vector <- unique(allDivTrt$site_proj_comm)

# #MPD over time by project
# for(PROJ in 1:length(site_proj_comm_vector)){
#   ggplot(data=filter(allDivYearMeans, site_proj_comm == site_proj_comm_vector[PROJ]),
#          aes(x=treatment_year, y=mpd.ses_mean, col=treatment,
#              ymin=mpd.ses_mean-mpd.ses_se, ymax=mpd.ses_mean+mpd.ses_se)) +
#     geom_errorbar(width=0.1) +
#     geom_point() +
#     geom_path() +
#     ggtitle(site_proj_comm_vector[PROJ]) +
#     theme_bw()
#   ggsave(filename=paste0("C:\\Users\\kjkomatsu\\Desktop\\pd figs\\",
#                          site_proj_comm_vector[PROJ], "_mpd.png"))
# }
# 
# #MNTD over time by project
# for(PROJ in 1:length(site_proj_comm_vector)){
#   ggplot(data=filter(allDivYearMeans, site_proj_comm == site_proj_comm_vector[PROJ]),
#          aes(x=treatment_year, y=mntd.ses_mean, col=treatment,
#              ymin=mntd.ses_mean-mntd.ses_se, ymax=mntd.ses_mean+mntd.ses_se)) +
#     geom_errorbar(width=0.1) +
#     geom_point() +
#     geom_path() +
#     ggtitle(site_proj_comm_vector[PROJ]) +
#     theme_bw()
#   ggsave(filename=paste0("C:\\Users\\kjkomatsu\\Desktop\\pd figs\\",
#                          site_proj_comm_vector[PROJ], "_mntd.png"))
# }


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
#   ggsave(filename=paste0("C:\\Users\\kjkomatsu\\Desktop\\pd figs\\",
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
#   ggsave(filename=paste0("C:\\Users\\kjkomatsu\\Desktop\\pd figs\\",
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
#   ggsave(filename=paste0("C:\\Users\\kjkomatsu\\Desktop\\pd figs\\",
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


# ### TO CONSIDER: should we look at each treatment within project and see if the slopes differ from 1 (more or less Fdiv per change in phy div), and does that differ by which treatment it is?
# #can we identify trajectories through time for each plot?
# 
# #plot trajectories through time, this is not very straight forward
# ggplot(data=subset(allDivTrt, project_name=='Fert1'&treatment %in% c('full_nut', 'control')), aes(x=mntd.ses, y=FDis, color=treatment, label=as.numeric(treatment_year))) +
#   geom_text() +
#   geom_path() +
#   facet_wrap(~plot_id)




##### causal modeling #####
allDivTrtFixed <- allDivTrt%>%
  mutate(site_year=paste(site_code, project_name, community_type, treatment_year, sep="::"))%>% #creates dummy site-by-year variable (time variant site effects)
  mutate(site_trt=paste(site_code, project_name, community_type, trt_type2, sep='::'))%>% #creates dummy site-by-plot variable (time invariant plot effects)
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::')) #creates dummy site-experiment variable

#individual treatments
#NOTE: these do account for env variables and trt magnitude
mpdDivCausalModel <- feols(mpd.ses ~ n + p + CO2 + drought + irg + temp + dist + multtrts | MAP + MAT + rrich + anpp, data=allDivTrtFixed, panel.id=~site_project_comm + treatment_year)
summary(mpdDivCausalModel, cluster="site_year") #trt effects
fixedEffects = fixef(mpdDivCausalModel) #env effects
plot(fixedEffects)

mntdDivCausalModel <- feols(mntd.ses ~ n + p + CO2 + drought + irg + temp + dist + multtrts | MAP + MAT + rrich + anpp, data=allDivTrtFixed, panel.id=~site_project_comm + treatment_year)
summary(mntdDivCausalModel, cluster="site_year") #trt effects
fixedEffects = fixef(mntdDivCausalModel) #env effects
plot(fixedEffects)

fdisDivCausalModel <- feols(FDis ~ n + p + CO2 + drought + irg + temp + dist + multtrts | MAP + MAT + rrich + anpp, data=allDivTrtFixed, panel.id=~site_project_comm + treatment_year)
summary(fdisDivCausalModel, cluster="site_year") #trt effects
fixedEffects = fixef(fdisDivCausalModel) #env effects
plot(fixedEffects)

etable(mpdDivCausalModel, mntdDivCausalModel, fdisDivCausalModel, cluster='site_year')

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

controlEnv <- control%>%
  group_by(site_code, project_name, community_type)%>%
  summarize_at(vars(mpd.ses_ctl_mean, mntd.ses_ctl_mean, FDis_ctl_mean, richness_ctl_mean), list(mean=mean), na.rm=T)%>%
  ungroup()%>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteBiotic.csv'))%>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteLocationClimate.csv'))%>%
  gather(key='env_variable', value='env_value', rrich, anpp, MAP, MAT, aridityValues)

#env driver effect in ctl plots
ggplot(data=controlEnv, aes(x=env_value, y=FDis_ctl_mean_mean)) +
  geom_point() +
  # geom_smooth(method='lm', se=F) +
  # geom_hline(yintercept=0) +
  xlab('') + ylab('Functional Dispersion') +
  facet_wrap(~env_variable, scales='free')

summary(lm(FDis_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='anpp'))) #no effect
summary(lm(FDis_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='rrich'))) #no effect
summary(lm(FDis_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='MAP'))) #no effect
summary(lm(FDis_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='MAT'))) #no effect
summary(lm(FDis_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='aridityValues'))) #no effect

ggplot(data=controlEnv, aes(x=env_value, y=mntd.ses_ctl_mean_mean)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  geom_hline(yintercept=0) +
  xlab('') + ylab('MNTD') +
  facet_wrap(~env_variable, scales='free')

summary(lm(mntd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='anpp'))) #significant negative relationship
summary(lm(mntd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='rrich'))) #significant positive relationship
summary(lm(mntd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='MAP'))) #no effect
summary(lm(mntd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='MAT'))) #significant negative relationship
summary(lm(mntd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='aridityValues'))) #marginal

# ggplot(data=controlEnv, aes(x=env_value, y=mpd.ses_ctl_mean_mean)) +
#   geom_point() +
#   geom_smooth(method='lm', se=F) +
#   geom_hline(yintercept=0) +
#   xlab('') + ylab('MPD') +
#   facet_wrap(~env_variable, scales='free')


allDivRR <- allDivTrt%>%
  filter(trt_type!='control')%>%
  left_join(control)%>%
  mutate(mpd_diff=(mpd.ses-mpd.ses_ctl_mean), mntd_diff=(mntd.ses-mntd.ses_ctl_mean), FDis_RR=((FDis-FDis_ctl_mean)/FDis_ctl_mean), richness_RR=((richness-richness_ctl_mean)/richness_ctl_mean))%>%
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteBiotic.csv'))%>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteLocationClimate.csv'))%>%
  gather(key=env_variable, value=env_value, c("MAP", "MAT", "rrich", "anpp", "aridityValues"))%>%
  group_by(site_proj_comm, site_code, project_name, community_type, treatment, trt_type2, env_variable)%>%
  summarise_at(vars(mpd_diff, mntd_diff, FDis_RR, richness_RR, env_value), list(mean=mean), na.rm=T)%>%
  ungroup()

#treatment responses by env drivers
ggplot(data=allDivRR, aes(x=env_value_mean, y=FDis_RR_mean)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  geom_hline(yintercept=0) +
  facet_grid(trt_type2~env_variable, scales='free')

ggplot(data=allDivRR, aes(x=env_value_mean, y=mntd_diff_mean)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  geom_hline(yintercept=0) +
  facet_grid(trt_type2~env_variable, scales='free')

##### mixed effects model #####
#NOTE: these models do not account for biotic or abiotic env drivers at a site or for trt magnitude
library(nlme)
library(emmeans)
library(performance)
library(grid)

options(contrasts=c('contr.sum','contr.poly')) 

summary(FDisModel <- lme(FDis_RR_mean ~ as.factor(trt_type2),
                         data=na.omit(subset(allDivRR, FDis_RR_mean<5 & trt_type2!='herb_removal')),
                         random=~1|site_proj_comm))
anova.lme(FDisModel, type='sequential')
meansFDisModel <- emmeans(FDisModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansFDisModelOutput <- as.data.frame(meansFDisModel$emmeans)

FDisFig <- ggplot(data=meansFDisModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.2) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('Functional Dispersion\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


summary(MNTDModel <- lme(mntd_diff_mean ~ as.factor(trt_type2),
                         data=na.omit(subset(allDivRR, trt_type2!='herb_removal')),
                         random=~1|site_proj_comm))
anova.lme(MNTDModel, type='sequential')
meansMNTDModel <- emmeans(MNTDModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansMNTDModelOutput <- as.data.frame(meansMNTDModel$emmeans)

MNTDFig <- ggplot(data=meansMNTDModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.2) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('Phylogenetic Distance\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


pushViewport(viewport(layout=grid.layout(1,2)))
print(FDisFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(MNTDFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))

### see sCoRRE_dCCA_traits script to find PCA of case study examples of extreme responses for each trt type
# N example of FDis and MNTD decreasing effect size: 	KUFS::E2::0::N1S0H0 or YMN::NitAdd::0::N80 (not CUL::Culardoch::0 N50 or maerc::fireplots::0 unuu)
# irrigation example of increasing FDis and MNTD effect size: KNZ::IRG::l i or SEV::WENNDEx::0 P or MNR::watfer::0::W
# drought example of increasing FDis and MNTD effect size: SFREC::GrazePrecip::G4 D or HAYS::Precip::0 reduction


