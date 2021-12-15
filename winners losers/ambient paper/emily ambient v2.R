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
    select(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block, plot_id, relcov, trt_type, species_matched, site_project_comm) %>% 
    pivot_wider(names_from="species_matched", values_from="relcov", values_fill=0)
  dat.keep=dat.filled %>% 
    pivot_longer(!c("site_code", "project_name", "community_type", "calendar_year", "treatment_year", "treatment", "block", "plot_id", "trt_type", "site_project_comm"), names_to="species_matched", values_to="relcov")
  myreldat_filled=rbind(myreldat_filled, dat.keep)
}

#get average relative cover for each species in a treatment, over all plots
relave<-myreldat_filled%>%
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, species_matched, calendar_year, treatment_year)%>%
  summarize(mean.relabund=mean(relcov))

#getting frequency of each plot type
myplots<-myreldat_filled%>%
  select(site_code, project_name, community_type, site_project_comm, treatment, block, trt_type, plot_id, calendar_year, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, calendar_year, treatment_year)%>%
  summarize(ntotplots=length(plot_id))

#getting number of plots of each type in which a species was present
freq<-myreldat_filled %>%
  select(site_code, project_name, community_type, site_project_comm, treatment, trt_type, species_matched, block, plot_id, calendar_year, treatment_year, relcov) %>%
  unique() %>%
  filter(relcov>0) %>% 
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, species_matched, calendar_year, treatment_year) %>%
  summarize(nplots=length(plot_id)) %>%
  left_join(myplots) %>%
  mutate(freq=nplots/ntotplots)

DCi.through.time<-relave %>%
  left_join(freq) %>%
  mutate(freq=replace_na(freq, 0)) %>%
  mutate(nplots=replace_na(nplots, 0)) %>%
  mutate(DCi=(mean.relabund+freq)/2) %>%
  select(site_code, project_name, community_type, site_project_comm, treatment, trt_type, species_matched, calendar_year, treatment_year, mean.relabund, nplots, freq, DCi)

filename=(paste(my.wd, "/WinnersLosers paper/DCi trends through time/ DCi trends through time.csv", sep=""))
write.csv(DCi.through.time, filename, row.names=F)

#add handy labels
DCi.through.time$site_project_comm=as.factor(paste(DCi.through.time$site_code, DCi.through.time$project_name, DCi.through.time$community_type, sep="_")); spc=unique(DCi.through.time$site_project_comm); length(spc)

DCi.through.time$site_project_comm_trt=as.factor(paste(DCi.through.time$site_project_comm, DCi.through.time$treatment, sep="::")); spct=unique(DCi.through.time$site_project_comm_trt); length(spct)

DCi.through.time$site_project_comm_sp=as.factor(paste(DCi.through.time$site_project_comm, DCi.through.time$species_matched, sep="::")); spcs=unique(DCi.through.time$site_project_comm_sp); length(spcs)

DCi.through.time$site_project_comm_sp_yr=as.factor(paste(DCi.through.time$site_project_comm, DCi.through.time$species_matched, DCi.through.time$calendar_year, sep="::")); spcsy=unique(DCi.through.time$site_project_comm_sp_yr); length(spcsy)

DCi.through.time$site_project_comm_sp_trt_yr=as.factor(paste(DCi.through.time$site_project_comm, DCi.through.time$species_matched, DCi.through.time$treatment, DCi.through.time$calendar_year, sep="::")); spcsty=unique(DCi.through.time$site_project_comm_sp_trt_yr); length(spcsty)

DCi.through.time$site_project_comm_sp_trt=as.factor(paste(DCi.through.time$site_project_comm, DCi.through.time$species_matched, DCi.through.time$treatment, sep="::")); spcst=unique(DCi.through.time$site_project_comm_sp_trt); length(spct)


#plot
for (i in 1:length(spc)) {
  qplot(treatment_year, DCi, data=DCi.through.time[DCi.through.time$site_project_comm==as.character(spc[i]),], color=treatment, main=spc[i]) + geom_smooth(method="lm", se=F, aes(group=treatment)) + facet_wrap(~species_matched)
  filename=paste(my.wd, "/WinnersLosers paper/DCi trends through time/", as.character(spc[i]), ".pdf", sep="")
  ggsave(filename, width=20, height=20)
}


####  how do we decide which species are rare enough to drop?  ###


for (i in 1:length(spc)) {
  ggplot(DCi.through.time[DCi.through.time$site_project_comm==as.character(spc[i]),], aes(treatment_year, DCi)) + geom_point() + geom_smooth(method="lm", se=F, aes(group=treatment)) + facet_wrap(~species_matched) + ggtitle(spc[i])
  filename=paste(my.wd, "/WinnersLosers paper/DCi trends through time/", as.character(spc[i]), ".pdf", sep="")
  ggsave(filename, width=20, height=20)
}

