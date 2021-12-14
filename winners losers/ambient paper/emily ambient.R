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
dat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RawAbundance_Dec2021.csv",sep=""))

#relative abundance data
reldat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Dec2021.csv",sep=""))

#info on treatments
trts<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfo_Dec2021.csv", sep=""))%>%
  select(site_code, project_name, community_type, treatment, trt_type, pulse, plot_mani, resource_mani)%>%
  unique()

{#for N and P other, do I include multiple nutrients or not?

trt_analysis<-trts%>%
  mutate(alltrts=ifelse(trt_type %in% c("CO2","CO2*temp", "mow_clip","burn","burn*graze","disturbance","burn*mow_clip","drought","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N","mult_nutrient","N*P","P","N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp","temp*mow_clip","drought*temp","irr*temp","irr"),1,0))%>%
  filter(alltrts==1)%>%
  mutate(dist=ifelse(trt_type %in% c("mow_clip","burn","burn*graze","disturbance","burn*mow_clip"), 1, 0), 
         dist_other=ifelse(trt_type %in% c("drought*mow_clip
","drought*temp*mow_clip", "irr*mow_clip","irr*temp*mow_clip","N*irr*mow_clip","N*P*burn*graze","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip", "mult_nutrient*mow_clip","N*burn*mow_clip", "N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp*mow_clip"), 1, 0),
         CO2=ifelse(trt_type %in% c("CO2"), 1, 0),
         CO2_other=ifelse(trt_type %in% c("CO2*temp","drought*CO2*temp","irr*CO2","irr*CO2*temp","N*CO2*temp","N*irr*CO2","N*irr*CO2*temp","N*CO2"), 1, 0),
         drought=ifelse(trt_type %in% c("drought"), 1, 0),
         drought_other=ifelse(trt_type %in% c("drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","N*drought","drought*temp"), 1, 0),
         herb_removal=ifelse(trt_type %in% c("herb_removal"), 1, 0),
         herb_removal_other=ifelse(trt_type %in% c("herb_removal*mow_clip","irr*herb_removal","N*herb_removal","P*herb_removal","mult_nutrient*herb_removal*mow_clip","mult_nutrient*herb_removal"), 1, 0),
         irg=ifelse(trt_type %in% c("irr"), 1, 0),
         irg_other=ifelse(trt_type %in% c("irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*irr*CO2","N*irr*mow_clip","mult_nutrient*irr","N*irr*CO2*temp","N*irr","N*irr*temp","irr*temp"), 1, 0),
         temp=ifelse(trt_type %in% c("temp"), 1, 0),
         temp_other=ifelse(trt_type %in% c("control","CO2*temp","drought*CO2*temp","drought*temp*mow_clip","irr*CO2*temp","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2*temp","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","temp*mow_clip","drought*temp","irr*temp"), 1, 0),
         nuts_other=ifelse(trt_type %in% c("control","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze","mult_nutrient*irr","N*irr*CO2*tempN*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip"), 1, 0),
         n=ifelse(trt_type=="N", 1, 0),
         n_other=ifelse(trt_type %in% c("N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze","N*irr*CO2*tempN*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","N*drought","N*herb_removal","N*irr","N*irr*temp","N*temp","N*P*temp","N*burn*mow_clip","N*P*burn","N*P*mow_clip"), 1, 0),
         p=ifelse(trt_type=="P", 1, 0),
         p_other=ifelse(trt_type %in% c("N*P*burn*graze","mult_nutrient*irr","P*burn*graze","P*burn*mow_clip","P*herb_removal","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal*mow_clip"), 1, 0),
         multtrts=ifelse(trt_type %in% c("mult_nutrient","N*P","CO2*temp", "burn*graze","burn*mow_clip","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp*mow_clip","drought*temp","irr*temp"),1,0))
  }

#cleaned species names
sp <-read.csv(paste(my.wd,"CoRRE data/CoRRE data/trait data/CoRRE2trykey_2021.csv", sep=""))%>%
  select(genus_species, species_matched)%>%
  unique

##DCi 

#combine relative abundance data with treatment and cleaned species names
allreldat<-reldat%>%
  left_join(trts)%>%
  left_join(sp)%>%# this drops the unknowns
  na.omit()

#get average relative cover for each species in a treatment, over all plots
relave<-allreldat%>%
  group_by(site_code, project_name, community_type, treatment, trt_type, plot_mani, species_matched, calendar_year, treatment_year)%>%
  summarize(mean=mean(relcov))


#work on control plots separately:

#subset out control plots, plot_mani==0
Crelave<-relave%>%
  filter(plot_mani==0)%>%
  ungroup()%>%
  select(-plot_mani)

#to get relative frequency, determine number of control plots
controlplots<-allreldat%>%
  filter(plot_mani==0)%>%
  select(site_code, project_name, community_type, plot_id, calendar_year, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, calendar_year, treatment_year)%>%
  summarize(ncplots=length(plot_id))

#to get relative frequency, determine number of control plots a species is found in, merge in total number of plots and calculate relative frequency  
control_freq<-allreldat%>%
  filter(plot_mani==0)%>%
  select(site_code, project_name, community_type, species_matched, treatment, trt_type, plot_id, calendar_year, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, species_matched, treatment, trt_type, calendar_year, treatment_year)%>%
  summarize(nplots=length(plot_id))%>%
  left_join(controlplots)%>%
  mutate(freq=nplots/ncplots)

#calculate DCi
control_dom<-control_freq%>%
  left_join(Crelave)%>%
  mutate(DCi=(mean+freq)/2)%>%
  select(site_code, project_name, community_type, treatment, trt_type, species_matched, calendar_year, treatment_year, mean, freq, DCi)%>%
  mutate(TorC="C")


#work on treatment plots separately:

#subset out treated plots, plot_mani!=0
Trelave<-relave%>%
  filter(plot_mani!=0)%>%
  ungroup()%>%
  select(-plot_mani)


#getting frequency of treated plots, same code as above but for treated plots
treatplots<-allreldat%>%
  filter(plot_mani!=0)%>%
  select(site_code, project_name, community_type, treatment, trt_type, plot_id, calendar_year, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment, trt_type, calendar_year, treatment_year)%>%
  summarize(ntplots=length(plot_id))

treat_freq<-allreldat%>%
  filter(plot_mani!=0)%>%
  select(site_code, project_name, community_type, treatment, trt_type, species_matched, plot_id, calendar_year, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment, trt_type, species_matched, calendar_year, treatment_year)%>%
  summarize(nplots=length(plot_id))%>%
  left_join(treatplots)%>%
  mutate(freq=nplots/ntplots)

treat_dom<-treat_freq%>%
  left_join(Trelave)%>%
  mutate(DCi=(mean+freq)/2)%>%
  select(site_code, project_name, community_type, treatment, trt_type, species_matched, calendar_year, treatment_year, mean, freq, DCi)%>%
  mutate(TorC="T")


#combine treatment and control plots:
DCi.through.time=rbind(treat_dom, control_dom)
DCi.through.time$site_project_comm_trt=as.factor(paste(DCi.through.time$site_code, DCi.through.time$project_name, DCi.through.time$community_type, DCi.through.time$trt_type, sep="::"))
DCi.through.time$year=as.factor(paste(DCi.through.time$calendar_year, DCi.through.time$treatment_year, sep="::")) #combine calendar year and treatment year columns into one to save this info to split out later
DCi.through.time$treatment_year=NULL
DCi.through.time$calendar_year=NULL


#assign DCi of zero when a species was absent from all replicate plots of a treatment or control (within a site/project/community/treatment)
spct=unique(DCi.through.time$site_project_comm_trt)
DCi.through.time_filled=NULL

for (j in 1:length(spct)) {
  dat=DCi.through.time[DCi.through.time$site_project_comm_trt==as.character(spct[j]),]
  dat.filled=dat %>% 
    select(site_code, project_name, community_type, treatment, trt_type, site_project_comm_trt, species_matched, year, DCi, TorC) %>% 
    pivot_wider(names_from="year", values_from="DCi", values_fill=0)
  dat.keep=dat.filled %>% 
    pivot_longer(!c("site_code", "project_name", "community_type", "treatment", "trt_type", "site_project_comm_trt", "species_matched", "TorC"), names_to="year", values_to="DCi") %>% 
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
