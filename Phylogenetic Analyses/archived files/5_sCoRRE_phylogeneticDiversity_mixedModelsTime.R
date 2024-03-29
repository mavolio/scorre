################################################################################
##  sCoRRE_phylogeneticDiversity_causalModels.R: Examining differences in phylogenetic and functional diversity with causal modeling of the CoRRE database.
##
##  Author: Kimberly Komatsu
##  Date created: December 4, 2022
################################################################################

library(data.table)
library(codyn)
library(fixest)
library(lme4)
library(tidyverse)


setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\paper 2_PD and FD responses\\data\\')  #kim's computer

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
trt <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\RawAbundance.csv') %>%
select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id) %>%
  unique() %>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\ExperimentInfo.csv')) %>%
  group_by(site_code, project_name, community_type) %>%
  mutate(experiment_length=max(treatment_year)) %>%
  ungroup() %>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, trt_type, experiment_length, plot_mani, n, p, CO2, precip, temp)

#phylogenetic diveristy data
pDiv <- read.csv('CoRRE_pd_metrics_non_weighted.csv') %>%
  separate(identifier, into=c("site_code", "project_name", "community_type", "treatment_year", "plot_id"), sep="::") %>%
  mutate(treatment_year=as.integer(treatment_year))

#functional diversity data
fDiv <- read.csv('CoRRE_functionalDiversity_2022-07-14.csv')

#species relative cover data
relCover <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeCover_Dec2021.csv') %>%
  mutate(replicate=paste(site_code, project_name, community_type, plot_id, sep='::')) #creating identifying column of each plot

#getting community diversity metrics for each plot
rDiv <- community_structure(relCover, time.var="treatment_year", abundance.var="relcov", replicate.var="replicate") %>%
  separate(replicate, into=c("site_code", "project_name", "community_type", "plot_id"), sep='::')

#merge all data on diversity metrics (phylogenetic, functional, species), experimental treatments, and site characteristics
allDiv <- pDiv %>% #phylogenetic metrics
  full_join(fDiv) %>% #functional metrics
  full_join(rDiv) %>% #species metrics
  filter(treatment_year>0) %>% #removing pre-treatment data
  mutate(pd.raw=ifelse(richness==1, 0, pd.raw), #mutating diversity metrics to be 0 if species richness is 1
         pd.ses=ifelse(richness==1, 0, pd.ses),
         mpd.raw=ifelse(richness==1, 0, mpd.raw),
         mpd.ses=ifelse(richness==1, 0, mpd.ses),
         mntd.raw=ifelse(richness==1, 0, mntd.raw),
         mntd.ses=ifelse(richness==1, 0, pd.raw),
         FDis=ifelse(richness==1, 0, FDis),
         RaoQ=ifelse(richness==1, 0, RaoQ)) %>%
  left_join(trt) %>% #treatments
  full_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteBiotic.csv')) %>% #site anpp and regional richness
  full_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteLocationClimate.csv')) %>% #site MAP and MAT
  mutate(site_proj_comm=paste(site_code,  project_name, community_type, sep='::')) %>%
  mutate(site_proj_comm_trt=paste(site_proj_comm, treatment, sep='::')) %>% 
  select(site_code, project_name, community_type, site_proj_comm, site_proj_comm_trt, treatment_year, calendar_year, treatment, plot_id, mpd.ses, mntd.ses, FDis, RaoQ, richness, Evar, trt_type, experiment_length, plot_mani, rrich, anpp, MAP, MAT, n, p, CO2, precip, temp) %>%
  filter(!(site_proj_comm %in% c('DL::NSFC::0', 'Naiman::Nprecip::0', 'DCGS::gap::0', 'Sil::NA::NA', 'SORBAS::NA::NA'))) %>%  #remove problem expts until they are fixed
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
  filter(!(site_proj_comm %in% c('DL::NSFC::0', 'Naiman::Nprecip::0', 'DCGS::gap::0'))) %>% #remove problem expts until they are fixed
  mutate_at(c('MAP', 'MAT', 'rrich', 'anpp'), funs(c(scale(.)))) %>% #scale site characteristics
  filter(!(is.na(site_proj_comm)))



# write.csv(allDivTrt, 'CoRRE_allDiversityMetrics_phyFunAnalysis.csv', row.names=F)



##### calculating response ratios #####
#filter control plots
control <- allDivTrt %>% 
  filter(trt_type=='control') %>%
  rename(mpd.ses_ctl=mpd.ses, mntd.ses_ctl=mntd.ses, FDis_ctl=FDis, richness_ctl=richness) %>%
  group_by(site_code, project_name, community_type, treatment_year) %>%
  summarize_at(vars(mpd.ses_ctl, mntd.ses_ctl, FDis_ctl, richness_ctl), list(mean=mean), na.rm=T) %>% #average aross plots and years
  ungroup()

#merge on site characteristics
controlEnv <- control %>%
  group_by(site_code, project_name, community_type) %>%
  summarize_at(vars(mpd.ses_ctl_mean, mntd.ses_ctl_mean, FDis_ctl_mean, richness_ctl_mean), list(mean=mean), na.rm=T) %>% #average aross plots and years
  ungroup() %>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteBiotic.csv')) %>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\siteLocationClimate.csv')) %>%
  gather(key='env_variable', value='env_value', rrich, anpp, MAP, MAT, aridityValues)

allDivRR <- allDivTrt %>%
  filter(trt_type!='control') %>%
  left_join(control) %>%
  mutate(mpd_diff=((mpd.ses-mpd.ses_ctl_mean)/mpd.ses_ctl_mean), mntd_diff=((mntd.ses-mntd.ses_ctl_mean)/mntd.ses_ctl_mean), FDis_RR=((FDis-FDis_ctl_mean)/FDis_ctl_mean), richness_RR=((richness-richness_ctl_mean)/richness_ctl_mean)) %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::'))  %>% 
  group_by(site_proj_comm, site_code, project_name, community_type, treatment, trt_type2, treatment_year) %>%
  summarise_at(vars(mpd_diff, mntd_diff, FDis_RR, richness_RR), list(mean=mean), na.rm=T)  %>%
  ungroup()


##### mixed effects model #####
#NOTE: these models do not account for biotic or abiotic env drivers at a site or for trt magnitude (but do include a random effect of site)
library(lme4)
library(car)
library(ggeffects)
# library(emmeans)
# library(performance)
library(grid)

options(contrasts=c('contr.sum','contr.poly')) 

summary(FDisModel <- lmer(FDis_RR_mean ~ poly(treatment_year,2)*as.factor(trt_type2) + (1 + trt_type2|site_proj_comm),
                         data=na.omit(subset(allDivRR, FDis_RR_mean<5 & trt_type2!='herb_removal' & treatment_year<50))))
Anova(FDisModel, type='III')
plot(FDisModel)
qqnorm(resid(FDisModel)) #not a great line


summary(MNTDModel <- lmer(mntd_diff_mean ~ poly(treatment_year,2)*as.factor(trt_type2) + (1 + trt_type2|site_proj_comm),
                         data=na.omit(subset(allDivRR, FDis_RR_mean<5 & trt_type2!='herb_removal' & treatment_year<50))))