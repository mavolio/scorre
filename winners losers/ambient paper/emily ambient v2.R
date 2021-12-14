rm(list=ls())

#investigate winners losers

library(tidyverse)
library(gridExtra)

###read in data

my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "E:/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "/Users/egrman/Dropbox/sDiv_sCoRRE_shared/"

#raw abundance data
dat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RawAbundance_Dec2021.csv",sep=""))

#relative abundance data
reldat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Dec2021.csv", sep=""))

#info on treatments
trts<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfo_Dec2021.csv", sep=""))%>%
  select(site_code, project_name, community_type, treatment, trt_type, pulse, resource_mani)%>%
  unique()

#cleaned species names
sp <-read.csv(paste(my.wd,"CoRRE data/CoRRE data/trait data/CoRRE2trykey_2021.csv", sep=""))%>%
  select(genus_species, species_matched)%>%
  unique

#combine relative abundance data with treatment and cleaned species names, drop mosses
myreldat<-reldat%>%
  left_join(trts)%>%
  left_join(sp)%>% #this drops the unknowns
  na.omit() %>% 
  mutate(drop=ifelse(species_matched %in% c("Andreaea obovata", "Anthelia juratzkana", "Aulacomnium turgidum", "Barbilophozia hatcheri", "Barbilophozia kunzeana", "Blepharostoma trichophyllum", "Brachythecium albicans", "Bryum arcticum", "Bryum pseudotriquetrum",  "Campylium stellatum", "Cyrtomnium hymenophyllum", "Dicranoweisia crispula", "Dicranum brevifolium", "Dicranum elongatum", "Dicranum fuscescens", "Dicranum groenlandicum",  "Dicranum scoparium", "Distichium capillaceum", "Ditrichum flexicaule", "Gymnomitrion concinnatum", "Hamatocaulis vernicosus", "Homalothecium pinnatifidum", "Hylocomium splendens", "Hypnum cupressiforme", "Hypnum hamulosum", "Isopterygiopsis pulchella", "Kiaeria starkei", "Leiocolea heterocolpos", "Marchantia polymorpha", "Marsupella brevissima", "Meesia uliginosa", "Myurella tenerrima", "Oncophorus virens", "Oncophorus wahlenbergii", "Pleurozium schreberi", "Pogonatum urnigerum", "Pohlia cruda", "Pohlia nutans", "Polytrichastrum alpinum", "Polytrichum juniperinum", "Polytrichum piliferum", "Polytrichum strictum", "Preissia quadrata", "Ptilidium ciliare", "Racomitrium lanuginosum", "Rhytidium rugosum", "Saelania glaucescens", "Sanionia uncinata",  "Schistidium apocarpum", "Syntrichia ruralis","Tomentypnum nitens", "Tortella tortuosa", "Tritomaria quinquedentata", "Nephroma arcticum", "Unknown NA", "Campylopus flexuosus", "Hypnum jutlandicum", "Plagiothecium undulatum", "Polytrichum commune", "Pseudoscleropodium purum", "Rhytidiadelphus loreus", "Rhytidiadelphus triquetrus", "Thuidium tamariscinum"), 1, 0)) %>% 
  filter(drop==0) %>% 
  group_by(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block, plot_id, trt_type, species_matched) %>% 
  summarize(relcov=sum(relcov)) %>% 
  mutate(site_project_comm=as.factor(paste(site_code, project_name, community_type, sep="::")))


#adding in zeros for species that were absent from a plot

spc=unique(myreldat$site_project_comm)
myreldat_filled=NULL

for (j in 1:length(spc)) {
  dat.filled=myreldat[myreldat$site_project_comm==as.character(spc[j]),] %>% 
    select(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block, plot_id, data_type, version, genus_species, relcov, trt_type, pulse, resource_mani, species_matched, site_project_comm) %>% 
    pivot_wider(names_from="species_matched", values_from="relcov", values_fill=0)
  dat.keep=dat.filled %>% 
    pivot_longer(!c("site_code", "project_name", "community_type", "calendar_year", "treatment_year", "treatment", "block", "plot_id", "data_type", "version", "genus_species", "trt_type", "pulse", "resource_mani", "site_project_comm"), names_to="species.matched", values_to="relcov")
  myreldat_filled=rbind(myreldat_filled, dat.keep)
}


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
