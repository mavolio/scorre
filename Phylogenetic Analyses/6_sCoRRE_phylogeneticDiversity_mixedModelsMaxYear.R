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


# ##### figures over time by site #####
# allDivYearMeans <- allDivTrt  %>%
#   group_by(site_proj_comm, site_code, project_name, community_type, treatment_year, treatment, trt_type2) %>%
#   summarize_at(vars(mpd.ses, mntd.ses, FDis, richness), list(mean=mean, se=se), na.rm=T) %>%
#   ungroup()
# 
# site_proj_comm_vector <- unique(allDivTrt$site_proj_comm)
# 
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
#   ggsave(filename=paste0("C:\\Users\\kjkomatsu\\Desktop\\pd figs\\mpd\\",
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
#   ggsave(filename=paste0("C:\\Users\\kjkomatsu\\Desktop\\pd figs\\mntd\\",
#                          site_proj_comm_vector[PROJ], "_mntd.png"))
# }



##### global trends - do environmental conditions impact FD or PD under ambient conditions? #####
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

# ##### mapping functional and phylogenetic diversity globally #####
# library('wesanderson')
# library("rnaturalearth")
# library("rnaturalearthdata")
# library("rgeos")
# world <- ne_countries(scale = "medium", returnclass = "sf")
# 
# z.pal<-wes_palette("Zissou1", 100, type = "continuous")
# 
# #map of functional dispersion
# ggplot(data=world)+
#   theme(panel.background=element_rect(fill="aliceblue", color="aliceblue"))+
#   theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
#         axis.text.y=element_text(size=20, colour="black"))+
#   geom_sf(color="black", fill="antiquewhite")+
#   geom_point(data=subset(controlEnv, env_variable=='rrich'), mapping=aes(x=Longitude,y=Latitude,fill=FDis_ctl_mean_mean),size=4.5,shape=21)+
#   scale_fill_gradientn(colours = z.pal) +
#   labs(fill="Phylogenetic Diversity (MNTD)")+
#   theme(legend.position = "top")+
#   ylab(expression("Latitude "*degree*""))+
#   xlab(expression("Longitude "*degree*""))
# 
# #map of phylogenetic diversity (MNTD)
# ggplot(data=world)+
#   theme(panel.background=element_rect(fill="aliceblue", color="aliceblue"))+
#   theme(text=element_text(size=20, colour="black"),axis.text.x=element_text(size=20, colour="black"),
#         axis.text.y=element_text(size=20, colour="black"))+
#   geom_sf(color="black", fill="antiquewhite")+
#   geom_point(data=subset(controlEnv, env_variable=='rrich'), mapping=aes(x=Longitude,y=Latitude,fill=mntd.ses_ctl_mean_mean),size=4.5,shape=21)+
#   scale_fill_gradientn(colours = z.pal) +
#   labs(fill="Phylogenetic Diversity (MNTD)")+
#   theme(legend.position = "top")+
#   ylab(expression("Latitude "*degree*""))+
#   xlab(expression("Longitude "*degree*""))

# ##### env drivers of FD #####
# ggplot(data=controlEnv, aes(x=env_value, y=FDis_ctl_mean_mean)) +
#   geom_point() +
#   # geom_smooth(method='lm', se=F) +
#   # geom_hline(yintercept=0) +
#   xlab('') + ylab('Functional Dispersion') +
#   facet_wrap(~env_variable, scales='free')
# 
# summary(lm(FDis_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='anpp'))) #no effect
# summary(lm(FDis_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='rrich'))) #no effect
# summary(lm(FDis_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='MAP'))) #no effect
# summary(lm(FDis_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='MAT'))) #no effect
# summary(lm(FDis_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='aridityValues'))) #no effect
# 
# #env drivers of MNTD
# ggplot(data=controlEnv, aes(x=env_value, y=mntd.ses_ctl_mean_mean)) +
#   geom_point() +
#   geom_smooth(method='lm', se=F) +
#   geom_hline(yintercept=0) +
#   xlab('') + ylab('MNTD') +
#   facet_wrap(~env_variable, scales='free')
# 
# summary(lm(mntd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='anpp'))) #significant negative relationship
# summary(lm(mntd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='rrich'))) #significant positive relationship
# summary(lm(mntd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='MAP'))) #no effect
# summary(lm(mntd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='MAT'))) #significant negative relationship
# summary(lm(mntd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='aridityValues'))) #marginal
# 
# #env drivers of MPD
# ggplot(data=controlEnv, aes(x=env_value, y=mpd.ses_ctl_mean_mean)) +
#   geom_point() +
#   geom_hline(yintercept=0) +
#   xlab('') + ylab('MPD') +
#   facet_wrap(~env_variable, scales='free')
# 
# summary(lm(mpd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='anpp'))) #no effect
# summary(lm(mpd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='rrich'))) #marginal
# summary(lm(mpd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='MAP'))) #no effect
# summary(lm(mpd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='MAT'))) #no effect
# summary(lm(mpd.ses_ctl_mean_mean~env_value, data=subset(controlEnv, env_variable=='aridityValues'))) #no effect



##### calculating response ratios #####
allDivRR <- allDivTrt %>%
  filter(trt_type!='control') %>%
  left_join(control) %>%
  mutate(mpd_diff=((mpd.ses-mpd.ses_ctl_mean)/mpd.ses_ctl_mean), mntd_diff=((mntd.ses-mntd.ses_ctl_mean)/mntd.ses_ctl_mean), FDis_RR=((FDis-FDis_ctl_mean)/FDis_ctl_mean), richness_RR=((richness-richness_ctl_mean)/richness_ctl_mean)) %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::'))  %>% 
  group_by(site_proj_comm, site_code, project_name, community_type, treatment, trt_type2) %>%
  mutate(max_year=ifelse(treatment_year==max(treatment_year), 1, 0)) %>% 
  ungroup() %>% 
  filter(max_year==1) #%>% 
  # group_by(site_proj_comm, site_code, project_name, community_type, treatment, trt_type2) %>% 
  # summarise_at(vars(mpd_diff, mntd_diff, FDis_RR, richness_RR), list(mean=mean), na.rm=T)  %>%
  # ungroup()
  

# ##### treatment responses by env drivers - no strong trends #####
# allDivRRenv <- allDivRR  %>% 
#   gather(key=env_variable, value=env_value, c("MAP", "MAT", "rrich", "anpp")) %>%
#   group_by(site_proj_comm, site_code, project_name, community_type, treatment, trt_type2, env_variable) %>%
#   ungroup()
# 
# ggplot(data=allDivRRenv, aes(x=env_value_mean, y=FDis_RR_mean)) +
#   geom_point() +
#   geom_smooth(method='lm', se=F) +
#   geom_hline(yintercept=0) +
#   facet_grid(trt_type2~env_variable, scales='free')
# 
# ggplot(data=allDivRRenv, aes(x=env_value_mean, y=mntd_diff_mean)) +
#   geom_point() +
#   geom_smooth(method='lm', se=F) +
#   geom_hline(yintercept=0) +
#   facet_grid(trt_type2~env_variable, scales='free')



##### mixed effects model #####
#NOTE: these models do not account for biotic or abiotic env drivers at a site or for trt magnitude (but do include a random effect of site)
library(nlme)
library(emmeans)
library(performance)
library(grid)

options(contrasts=c('contr.sum','contr.poly')) 

summary(FDisModel <- lme(FDis_RR ~ as.factor(trt_type2),
                         data=na.omit(subset(allDivRR, FDis_RR<5 & trt_type2!='herb_removal')),
                         random=~1|site_proj_comm))
anova.lme(FDisModel, type='sequential')
meansFDisModel <- emmeans(FDisModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansFDisModelOutput <- as.data.frame(meansFDisModel$emmeans)

FDisFig <- ggplot(data=meansFDisModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('Functional Dispersion\nEffect Size') + xlab('') +
  scale_x_discrete(limits=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), breaks=c('multiple trts', 'disturbance', 'temp', 'drought', 'CO2', 'irr', 'P', 'N'), labels=c('Multiple Trts', 'Disturbance', 'Temperature', 'Drought', 'CO2','Irrigation', 'P', 'N')) + 
  scale_color_manual(values=c('blue', 'orange', 'orange', 'blue', 'dark grey', 'blue', 'blue', 'orange')) +
  theme(legend.position='none')


summary(MNTDModel <- lme(mntd_diff ~ as.factor(trt_type2),
                         data=na.omit(subset(allDivRR, trt_type2!='herb_removal')),
                         random=~1|site_proj_comm))
anova.lme(MNTDModel, type='sequential')
meansMNTDModel <- emmeans(MNTDModel, pairwise~as.factor(trt_type2), adjust="tukey")
meansMNTDModelOutput <- as.data.frame(meansMNTDModel$emmeans)

MNTDFig <- ggplot(data=meansMNTDModelOutput, aes(x=trt_type2, y=emmean, color=trt_type2)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=emmean-SE*1.96, ymax=emmean+SE*1.96), width=0.2) +
  geom_hline(yintercept=0) +
  coord_flip() +
  ylab('Phylogenetic Diversity (MNTD)\nEffect Size') + xlab('') +
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


##### treatment magnitude #####
allDivRRtrt <- allDivRR %>% 
  left_join(trt) %>% 
  select(-treatment_year, -calendar_year, -plot_id) %>% 
  unique()

#N additions
nDivRR <- allDivRRtrt %>% 
  filter(trt_type2=='N')

summary(nFDisModel <- lme(FDis_RR_mean ~ n, #tried polynomial, but linear was better
                          data=na.omit(subset(nDivRR, FDis_RR_mean<5)),
                          random=~1|site_proj_comm))
anova.lme(nFDisModel, type='sequential')
coef(nFDisModel)
summary.tablefunc(nFDisModel)

nFDisFig <- ggplot(data=nDivRR, aes(x=n, y=FDis_RR_mean)) +
  geom_point(size=2)+
  # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
  geom_smooth(method='lm', formula=y~x, color='black') +
  geom_hline(yintercept=0) +
  coord_cartesian(ylim=c(-0.75,1.25)) +
  ylab('Functional Dispersion\nEffect Size') + xlab(bquote('N added '(gm^-2)))



summary(nMNTDModel <- lme(mntd_diff_mean ~ n, #tried polynomial, but linear was better
                         data=nDivRR,
                         random=~1|site_proj_comm))
anova.lme(nMNTDModel, type='sequential')
coef(nMNTDModel)
summary.tablefunc(nMNTDModel)

nMNTDFig <- ggplot(data=nDivRR, aes(x=n, y=mntd_diff_mean)) +
  geom_point(size=2)+
  # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
  geom_smooth(method='lm', formula=y~x, color='black') +
  geom_hline(yintercept=0) +
  coord_cartesian(ylim=c(-0.75,1.25)) +
  ylab('Phylogenetic Diversity (MNTD)\nEffect Size') + xlab(bquote('N added '(gm^-2)))


#precip
precipDivRR <- allDivRRtrt %>% 
  filter(precip!=0)

summary(precipFDisModel <- lme(FDis_RR_mean ~ precip, #tried a polynomial model, but linear was better fit
                          data=na.omit(subset(precipDivRR, FDis_RR_mean<5)),
                          random=~1|site_proj_comm))
anova.lme(precipFDisModel, type='sequential')
coef(precipFDisModel)
summary.tablefunc(precipFDisModel)

precipFDisFig <- ggplot(data=precipDivRR, aes(x=precip, y=FDis_RR_mean)) +
  geom_point(size=2)+
  # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
  geom_smooth(method='lm', formula=y~x, color='black') +
  geom_hline(yintercept=0) +
  coord_cartesian(ylim=c(-0.5,1.25)) +
  ylab('Functional Dispersion\nEffect Size') + xlab('Precipitation Manipulation (%)')



summary(precipMNTDModel <- lme(mntd_diff_mean ~ precip,
                          data=precipDivRR,
                          random=~1|site_proj_comm))
anova.lme(precipMNTDModel, type='sequential')
coef(precipMNTDModel)
summary.tablefunc(precipMNTDModel)

precipMNTDFig <- ggplot(data=precipDivRR, aes(x=precip, y=mntd_diff_mean)) +
  geom_point(size=2)+
  # geom_abline(linewidth=2, aes(intercept=`(Intercept)`, slope=`poly(n, 2)`, as.data.frame(t(fixef(nMNTDModel))))) +
  # geom_smooth(method='lm', formula=y~poly(x,2), color='black') + #no significant effect
  geom_hline(yintercept=0) +
  coord_cartesian(ylim=c(-0.5,1.25)) +
  ylab('Phylogenetic Diversity (MNTD)\nEffect Size') + xlab('Precipitation Manipulation (%)')


#combined N and precip magnitude figure
pushViewport(viewport(layout=grid.layout(2,2)))
print(nFDisFig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(nMNTDFig, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(precipFDisFig, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(precipMNTDFig, vp=viewport(layout.pos.row=2, layout.pos.col=2))


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



# ##### trends in PD and FD - vote counting across experiments #####
# 
# anovaResults=data.frame(row.names=1) #creates dataframe to put results into
# 
# forANOVA <- allDiv %>% 
#   select(site_proj_comm) %>% 
#   unique()
# 
# 
# #running an anova for each treatment independently
# #includes one for functional dispersion response and one for phylogenetic diversity
# for(i in 1:length(forANOVA$site_proj_comm)) {
#   
#   #creates a dataset for each unique year, trt, exp combo
#   subset <- allDiv[allDiv$site_proj_comm==as.character(forANOVA$site_proj_comm[i]),] %>%
#     select(site_proj_comm_trt, calendar_year, treatment, plot_mani, plot_id, mntd.ses, FDis) %>%
#     mutate(treatment2=ifelse(plot_mani==0, 'TRUECONTROL', as.character(treatment)))
#   
#   #control data only
#   controls <- subset %>% 
#     filter(plot_mani==0)
#   
#   #treatments only
#   treatments <- subset %>%
#     filter(plot_mani!=0) %>% 
#     select(treatment) %>%
#     unique()
#   
#   #splitting into separate anovas for each treatment
#   for(j in 1:length(treatments$treatment)) {
#     
#     #bind control data to treatment data
#     subset2 <- subset[subset$treatment==as.character(treatments$treatment[j]),] %>% 
#       rbind(controls)
#     
#     #treatment label
#     trt_label=as.character(treatments$treatment[j])
#     site_proj_comm_label=as.character(forANOVA$site_proj_comm[i])
#     
#     
#     #ANOVAs
#     anovaFDis <- as.data.frame(summary(aov(FDis~treatment, data=subset2))[[1]]) %>% #generate anova statistics
#       mutate(treatment=trt_label, site_proj_comm=site_proj_comm_label) %>% #add treatment and site_proj_comm labels
#       mutate(variable='FDis')
#     anovaFDis <- anovaFDis[rownames(anovaFDis) %in% c('treatment  '), ]
#     
#     anovaMNTD <- as.data.frame(summary(aov(mntd.ses~treatment, data=subset2))[[1]]) %>% #generate anova statistics
#       mutate(treatment=trt_label, site_proj_comm=site_proj_comm_label) %>% #add treatment and site_proj_comm labels
#       mutate(variable='MNTD')
#     anovaMNTD <- anovaMNTD[rownames(anovaMNTD) %in% c('treatment  '), ]
#       
#     #bind data together
#     anovaResults <- rbind(anovaResults, anovaFDis, anovaMNTD)
#   }
# }
# 
# #generate dataframe of the response ratios
# allDivRRpivot <- allDivRR %>% 
#   select(site_proj_comm, site_code, project_name, community_type, treatment, trt_type2, mntd_diff_mean, FDis_RR_mean) %>%
#   pivot_longer(cols=mntd_diff_mean:FDis_RR_mean, names_to='var', values_to='RR') %>% 
#   mutate(variable=ifelse(var=='mntd_diff_mean', 'MNTD', 'FDis')) %>% 
#   select(-var)
# 
# names(anovaResults) <- c('df', 'sum_sq', 'mean_sq', 'F_value', 'p_value', 'treatment', 'site_proj_comm', 'variable')
# 
# anovaResultsCorrected <- anovaResults %>%  #need to figure out the NAs (were these where there is no RR because richness is 1?)
#   left_join(allDivRRpivot) %>% 
#   mutate(effect=ifelse(p_value<0.00008 & RR>0, 'increasing',
#                        ifelse(p_value<0.00008 & RR<0, 'decreasing', 'no effect'))) %>%  
#   #bonferonni correction: p<0.00008 for 622 comparisons (for each variable type)
#   select(treatment, site_proj_comm, variable, site_code, project_name, community_type, trt_type2, effect) %>% 
#   pivot_wider(names_from='variable', values_from='effect') %>% 
#   filter(!is.na(FDis), !is.na(MNTD)) %>%  #lose 34 treatments to not having complete cases
#   mutate(effect=paste(FDis, MNTD, sep='::'))
# 
# #get numbers of each kind of response
# summary(as.factor(anovaResultsCorrected$effect))
# # decreasing::decreasing  decreasing::no effect increasing::decreasing  increasing::no effect  no effect::decreasing  no effect::increasing 
# # 34                     15                      3                     25                     44                      5 
# # no effect::no effect 
# # 462






































