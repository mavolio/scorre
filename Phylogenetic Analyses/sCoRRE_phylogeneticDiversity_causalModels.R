################################################################################
##  sCoRRE_phylogeneticDiversity_causalModels.R: Examining differences in phylogenetic and functional diversity with causal modeling of the CoRRE database.
##
##  Author: Kimberly Komatsu
##  Date created: December 13, 2021
################################################################################

library(data.table)
library(fixest)
library(lme4)
library(tidyverse)


setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\paper 2_PD and FD responses\\data\\')  #kim's laptop

##### functions and themes #####

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


pDiv <- read.csv('CoRRE_pd_metrics_weighted_des2021_kjkedit.csv')%>%
  left_join(trt)%>%
  full_join(read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteBiotic.csv'))%>%
  full_join(read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteLocationClimate.csv'))

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
  mutate(trt_type2=ifelse(dist==1, 'disturbance', ifelse(multtrts==1, 'multiple trts', trt_type)))%>%
  select(site_code, project_name, community_type, treatment, alltrts, dist, tCO2, drought, therb_removal, irg, ttemp, tn, tp, multtrts, trt_type2)%>%
  unique()


#pick treatments here.
pDivTrt <- pDiv%>%
  right_join(trt_analysis)

# write.csv(pDivAvgTrt, 'PD_mean_effectsizes_20211214.csv')




##### causal modeling - difference #####
pDivTrtFixed <- pDivTrt%>%
  mutate(site_year=paste(site_code, project_name, community_type, treatment_year, sep="::"))%>% #creates dummy site-by-year variable (time variant site effects)
  mutate(site_trt=paste(site_code, project_name, community_type, trt_type2, sep='::'))%>% #creates dummy site-by-plot variable (time invariant plot effects)
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='::')) #creates dummy site-experiment variable

#individual treatments
pDivCausalModel <- feols(mpd.ses ~ n + p + k+  CO2 + drought + irg + temp + dist | MAP + MAT + rrich + anpp, data=pDivTrtFixed, panel.id=~site_project_comm + treatment_year)
etable(pDivCausalModel, cluster="site_year") #trt effects
fixedEffects = fixef(pDivCausalModel) #env effects
plot(fixedEffects)

# #all together
# pDivCausalModel <- feols(mpd_diff_avg ~ trt_type2 | site_year, data=pDivAvgTrtFixed)
# etable(pDivCausalModel,
#        cluster="site_year")

##### figures #####
ggplot(data=pDivAvgTrtFixed, aes(x=trt_type2, y=mpd_diff_avg)) +
  geom_boxplot() +
  geom_hline(yintercept=0)

ggplot(data=subset(pDivAvgTrtFixed, n<150), aes(x=n, y=mpd_diff_avg, color=rrich)) +
  geom_point() +
  geom_hline(yintercept=0)
