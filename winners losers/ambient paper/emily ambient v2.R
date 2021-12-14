#investigate winners losers

library(tidyverse)
library(gridExtra)

###read in data

my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "E:/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "/Users/egrman/Dropbox/sDiv_sCoRRE_shared/"

#read in the data

#raw abundance data
dat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RawAbundance_Dec2021.csv", sep=""))

#relative abundance data
reldat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Dec2021.csv", sep=""))

#info on treatments
trts<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfo_Dec2021.csv", sep=""))%>%
  select(site_code, project_name, community_type, treatment, trt_type, pulse, resource_mani)%>%
  unique()

#DCi

#combine relative abundance data with treatment and cleaned species names
allreldat<-reldat%>%
  left_join(trts)%>%
  left_join(sp)%>% #this drops the unknowns
  na.omit()

#get average relative cover for each species in a treatment, over all plots
relave<-allreldat%>%
  group_by(site_code, project_name, community_type, treatment, trt_type, species_matched, calendar_year, treatment_year)%>%
  summarize(mean.relabund=mean(relcov))

#getting frequency of each plot type
allplots<-allreldat%>%
  select(site_code, project_name, community_type, treatment, trt_type, plot_id, calendar_year, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment, trt_type, calendar_year, treatment_year)%>%
  summarize(ntplots=length(plot_id))

freq<-allreldat%>%
  select(site_code, project_name, community_type, treatment, trt_type, species_matched, plot_id, calendar_year, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment, trt_type, species_matched, calendar_year, treatment_year)%>%
  summarize(nplots=length(plot_id))%>%
  left_join(allplots)%>%
  mutate(freq=nplots/ntplots)

DCi.through.time<-freq%>%
  left_join(relave)%>%
  mutate(DCi=(mean.relabund+freq)/2)%>%
  select(site_code, project_name, community_type, treatment, trt_type, species_matched, calendar_year, treatment_year, mean.relabund, freq, DCi)

#add handy labels
DCi.through.time$site_project_comm=as.factor(paste(DCi.through.time$site_code, DCi.through.time$project_name, DCi.through.time$community_type, sep="::"))
DCi.through.time$year=as.factor(paste(DCi.through.time$calendar_year, DCi.through.time$treatment_year, sep="::")) #combine calendar year and treatment year columns into one to save this info to split out later
DCi.through.time$treatment_year=NULL
DCi.through.time$calendar_year=NULL


#assign DCi of zero when a species was absent from all replicate plots of a treatment or control (within a site/project/community/treatment)
spc=unique(DCi.through.time$site_project_comm)
DCi.through.time_filled=NULL

for (j in 1:length(spc)) {
  dat.filled=DCi.through.time[DCi.through.time$site_project_comm==as.character(spc[j]),] %>% 
    select(site_code, project_name, community_type, treatment, trt_type, site_project_comm, species_matched, year, DCi) %>% 
    pivot_wider(names_from="year", values_from="DCi", values_fill=0)
  dat.keep=dat.filled %>% 
    pivot_longer(!c("site_code", "project_name", "community_type", "treatment", "trt_type", "site_project_comm_trt", "species_matched"), names_to="year", values_to="DCi") %>% 
    separate(year, c("calendar_year", "treatment_year"))
  DCi.through.time_filled=rbind(DCi.through.time_filled, dat.keep)
}

#add handy labels

DCi.through.time_filled$site_project_comm=as.factor(paste(DCi.through.time_filled$site_code, DCi.through.time_filled$project_name, DCi.through.time_filled$community_type, sep="::"))
DCi.through.time_filled$site_project_comm_sp=as.factor(paste(DCi.through.time_filled$site_project_comm, DCi.through.time_filled$species_matched, sep="::"))


#plot
for (i in 1:length(spcs)) {
  qplot(treatment_year, DCi, data=DCi.through.time_filled[DCi.through.time_filled$site_project_comm_sp==as.character(spcs[i]),], color=treatment, shape=TorC, main=spcs[i]) + geom_smooth(method="lm", se=F, aes(group=treatment))
  
  
}
