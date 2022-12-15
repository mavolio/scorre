#### Emily Grman modified emily DCi trends through time.R
#### Dec 13, 2022 sCoRRE meeting 4 at iDiv



rm(list=ls())

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
sp <-read.csv(paste(my.wd,"CoRRE data/trait data/CoRRE2trykey_2021.csv", sep=""))%>%
  select(genus_species, species_matched)%>%
  unique

#reading in categorical traits
my_cat=read.csv(paste(my.wd, "CoRRE data/trait data/sCoRRE categorical trait data_11302021.csv", sep="")) %>%
  select(species_matched, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal, n_fixation)

#combine relative abundance data with treatment, cleaned species names, categorical traits, of life forms keep only graminoids, forbs, vines, woody.
myreldat<-reldat%>%
  left_join(trts)%>%
  left_join(sp)%>% #this drops the unknowns??
  left_join(my_cat)%>% # this adds ~5k rows. why?
  na.omit()%>% 
  mutate(drop=ifelse(species_matched %in% c("Andreaea obovata", "Anthelia juratzkana", "Aulacomnium turgidum", "Barbilophozia hatcheri", "Barbilophozia kunzeana", "Blepharostoma trichophyllum", "Brachythecium albicans", "Bryum arcticum", "Bryum pseudotriquetrum", "Campylium stellatum", "Cyrtomnium hymenophyllum", "Dicranoweisia crispula", "Dicranum brevifolium", "Dicranum elongatum", "Dicranum fuscescens", "Dicranum groenlandicum",  "Dicranum scoparium", "Distichium capillaceum", "Ditrichum flexicaule", "Gymnomitrion concinnatum", "Hamatocaulis vernicosus", "Homalothecium pinnatifidum", "Hylocomium splendens", "Hypnum cupressiforme", "Hypnum hamulosum", "Isopterygiopsis pulchella", "Kiaeria starkei", "Leiocolea heterocolpos", "Marchantia polymorpha", "Marsupella brevissima", "Meesia uliginosa", "Myurella tenerrima", "Oncophorus virens", "Oncophorus wahlenbergii", "Pleurozium schreberi", "Pogonatum urnigerum", "Pohlia cruda", "Pohlia nutans", "Polytrichastrum alpinum", "Polytrichum juniperinum", "Polytrichum piliferum", "Polytrichum strictum", "Preissia quadrata", "Ptilidium ciliare", "Racomitrium lanuginosum", "Rhytidium rugosum", "Saelania glaucescens", "Sanionia uncinata",  "Schistidium apocarpum", "Syntrichia ruralis","Tomentypnum nitens", "Tortella tortuosa", "Tritomaria quinquedentata", "Nephroma arcticum", "Unknown NA", "Campylopus flexuosus", "Hypnum jutlandicum", "Plagiothecium undulatum", "Polytrichum commune", "Pseudoscleropodium purum", "Rhytidiadelphus loreus", "Rhytidiadelphus triquetrus", "Thuidium tamariscinum"), 1, ifelse(growth_form %in% c("forb", "graminoid", "vine"), 0, 1))) %>% 
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
  summarize(mean.relabund=mean(relcov)) %>% 
  ungroup()

#getting frequency of each plot type
myplots<-myreldat_filled%>%
  select(site_code, project_name, community_type, site_project_comm, treatment, block, trt_type, plot_id, calendar_year, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, calendar_year, treatment_year)%>%
  summarize(ntotplots=length(plot_id)) %>% 
  ungroup()

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

filename=(paste(my.wd, "/WinnersLosers paper/DCi/ DCi trends through time.csv", sep="")); write.csv(DCi.through.time, filename, row.names=F)

#DCi.through.time=read.csv(paste(my.wd, "/WinnersLosers paper/DCi/ DCi trends through time.csv", sep=""))

#add handy labels
DCi.through.time$site_project_comm=as.factor(paste(DCi.through.time$site_code, DCi.through.time$project_name, DCi.through.time$community_type, sep="_")); spc=unique(DCi.through.time$site_project_comm); length(spc)

#plot species trends through time, omitting pre-treatment data

for (i in 1:length(spc)) {
  qplot(treatment_year, DCi, data=DCi.through.time[DCi.through.time$site_project_comm==as.character(spc[i]) & DCi.through.time$treatment_year>0,], color=treatment, main=spc[i]) + geom_smooth(method="lm", se=F, aes(group=treatment)) + facet_wrap(~species_matched)
  filename=paste(my.wd, "/WinnersLosers paper/DCi/species trends through time/", as.character(spc[i]), ".pdf", sep="")
  ggsave(filename, width=20, height=20)
}


#averaging over time (omitting pre-treatment data) and adding column for length of study:

study.length=DCi.through.time %>% 
  select(site_code, project_name, community_type, site_project_comm, treatment, trt_type, treatment_year) %>% 
  filter(treatment_year>0) %>% 
  unique %>% 
  group_by(site_code, project_name, community_type, treatment, trt_type) %>% 
  summarize(study_length=max(treatment_year), number_years=length(treatment_year)) %>% 
  ungroup()

DCi.averaged=DCi.through.time %>% 
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, species_matched) %>% 
  filter(treatment_year>0) %>% 
  summarize(across(c(DCi, freq, mean.relabund, nplots), mean)) %>% 
  ungroup() %>% 
  right_join(study.length)

filename=(paste(my.wd, "/WinnersLosers paper/DCi/ DCi averaged through time.csv", sep=""))
write.csv(DCi.averaged, filename, row.names=F)

#making list of studies that have at least 5 years of data (and should be kept in future analysis):
long_enough_studies=DCi.averaged %>% 
  select(site_code, project_name, community_type, site_project_comm, treatment, trt_type, study_length, number_years) %>% 
  unique() %>% 
  filter(!number_years<5) 

# how do we decide which species are rare enough to drop? 

for (i in 1:length(spc)) {
  ggplot(DCi.averaged[DCi.averaged$site_project_comm==as.character(spc[i]),], aes(DCi)) + geom_histogram(aes(fill=treatment)) + ggtitle(spc[i])
  filename=paste(my.wd, "/WinnersLosers paper/DCi/histograms to look at rarity/", as.character(spc[i]), ".pdf", sep="")
  ggsave(filename, width=10, height=10)
}

#exploring DCi cutoffs: keep DCi>0.1

DCi_0.1=DCi.averaged %>% 
  select(site_code, project_name, community_type, site_project_comm, treatment, trt_type, species_matched, DCi, freq, mean.relabund, nplots) %>% 
  filter(DCi<0.1)

ggplot(DCi_0.1, aes(freq)) + geom_histogram()
filename=paste(my.wd, "/WinnersLosers paper/DCi/freq histogram with DCi less than 0.1", ".pdf", sep=""); ggsave(filename, width=6, height=6)

ggplot(DCi_0.1, aes(mean.relabund)) + geom_histogram() 
filename=paste(my.wd, "/WinnersLosers paper/DCi/relabund histogram with DCi less than 0.1", ".pdf", sep=""); ggsave(filename, width=6, height=6)

DCi_0.1=DCi.averaged %>% 
  select(site_code, project_name, community_type, site_project_comm, treatment, trt_type, species_matched, DCi, freq, mean.relabund, nplots) %>% 
  filter(!DCi<0.1)

ggplot(DCi_0.1, aes(freq)) + geom_histogram()
filename=paste(my.wd, "/WinnersLosers paper/DCi/freq histogram with DCi greater than 0.1", ".pdf", sep=""); ggsave(filename, width=6, height=6)

ggplot(DCi_0.1, aes(mean.relabund)) + geom_histogram() 
filename=paste(my.wd, "/WinnersLosers paper/DCi/relabund histogram with DCi greater than 0.1", ".pdf", sep=""); ggsave(filename, width=6, height=6)


#cutting dataset to only species with average (across all plots all years) DCi >=0.1 and to studies with 5 or more years of data (omitting pre-treatment data)

DCi_dropped=DCi_0.1 %>% 
  select(site_code, project_name, community_type, site_project_comm, species_matched) %>% 
  unique %>% 
  left_join(DCi.through.time) %>% 
  right_join(long_enough_studies)

filename=(paste(my.wd, "/WinnersLosers paper/DCi/ DCi through time dropping rare species and short studies.csv", sep="")); write.csv(DCi_dropped, filename, row.names=F)

#plotting species trends through time (omitting pre-treatment data)

spc.dropped=unique(DCi_dropped$site_project_comm)
for (i in 1:length(spc.dropped)) {
  qplot(treatment_year, DCi, data=DCi_dropped[DCi_dropped$site_project_comm==as.character(spc.dropped[i]) & DCi_dropped$treatment_year>0,], color=treatment, main=spc.dropped[i]) + geom_smooth(method="lm", se=F, aes(group=treatment)) + facet_wrap(~species_matched)
  filename=paste(my.wd, "/WinnersLosers paper/DCi/species trends through time, dropping rare and short studies/", as.character(spc.dropped[i]), ".pdf", sep="")
  ggsave(filename, width=10, height=10)
}

#running regressions for each species each treatment and extracting slopes (omitting pre-treatment data)

DCi_dropped$site_project_comm_trt_sp=as.factor(paste(DCi_dropped$site_project_comm, DCi_dropped$treatment, DCi_dropped$species_matched, sep="::"))

spcts=unique(DCi_dropped$site_project_comm_trt_sp)
spcts_regressions=NULL

for (h in 1:length(spcts)) {
  mod=lm(DCi~treatment_year, data=DCi_dropped[DCi_dropped$site_project_comm_trt_sp==as.character(spcts[h]) & DCi_dropped$treatment_year>0,])
  mod2=lm(mean.relabund~treatment_year, data=DCi_dropped[DCi_dropped$site_project_comm_trt_sp==as.character(spcts[h]) & DCi_dropped$treatment_year>0,])
  keep=cbind(as.character(spcts[h]), mod$coefficients[[2]], anova(mod)$Pr[1], mod2$coefficients[[2]], anova(mod2)$Pr[1])
  spcts_regressions=rbind(spcts_regressions, keep)
}

spcts_regressions=data.frame(spcts_regressions)
names(spcts_regressions)=c("site_project_comm_trt_sp", "DCi.slope", "DCi.p", "relabund.slope", "relabund.p")

DCi_regressions=DCi_dropped %>% 
  select(site_code, project_name, community_type, site_project_comm, species_matched, treatment, trt_type, site_project_comm_trt_sp, study_length, number_years) %>% 
  unique() %>% 
  left_join(spcts_regressions)

filename=(paste(my.wd, "/WinnersLosers paper/DCi/ regression results.csv", sep="")); write.csv(DCi_regressions, filename, row.names=F)




#making some more figures to check out our DCi cutoff (to compare with figures Meghan had previously made)

#averaging again across treatments, leaving control plots out:

DCi.av.control=DCi.averaged[DCi.averaged$trt_type=="control",]
names(DCi.av.control)[names(DCi.av.control)=="DCi"]<-"DCi.control"
names(DCi.av.control)[names(DCi.av.control)=="mean.relabund"]<-"mean.relabund.control"
names(DCi.av.control)[names(DCi.av.control)=="freq"]<-"freq.control"

DCi.av.trt=DCi.averaged[!DCi.averaged$trt_type=="control",]
names(DCi.av.trt)[names(DCi.av.trt)=="DCi"]<-"DCi.trt"
names(DCi.av.trt)[names(DCi.av.trt)=="mean.relabund"]<-"mean.relabund.trt"
names(DCi.av.trt)[names(DCi.av.trt)=="freq"]<-"freq.trt"

DCi.av=merge(DCi.av.trt, DCi.av.control, by=c("site_code", "project_name", "community_type", "site_project_comm", "site_project_comm_sp", "species_matched"))
DCi.av$sp=ifelse(DCi.av$species_matched=="Andropogon gerardii", "Andropogon gerardii", ifelse(DCi.av$species_matched=="Ambrosia psilostachya", "Ambrosia psilostachya", ifelse(DCi.av$species_matched=="Triodanis perfoliata", "Triodanis perfoliata", "a")))

ggplot(data=DCi.av, aes(DCi.control, DCi.trt)) + geom_point(color="lightgray") + geom_point(data=DCi.av[DCi.av$sp %in% c("Ambrosia psilostachya", "Andropogon gerardii", "Triodanis perfoliata"),], aes(color=sp)) + geom_abline(intercept=0, slope=1)
filename=paste(my.wd, "/WinnersLosers paper/DCi/DCi correlation", ".pdf", sep=""); ggsave(filename, width=10, height=10)

ggplot(data=DCi.av, aes(freq.control, freq.trt)) + geom_point(color="lightgray") + geom_point(data=DCi.av[DCi.av$sp %in% c("Ambrosia psilostachya", "Andropogon gerardii", "Triodanis perfoliata"),], aes(color=sp)) + geom_abline(intercept=0, slope=1)
filename=paste(my.wd, "/WinnersLosers paper/DCi/freq correlation", ".pdf", sep=""); ggsave(filename, width=10, height=10)

ggplot(data=DCi.av, aes(mean.relabund.control, mean.relabund.trt)) + geom_point(color="lightgray") + geom_point(data=DCi.av[DCi.av$sp %in% c("Ambrosia psilostachya", "Andropogon gerardii", "Triodanis perfoliata"),], aes(color=sp)) + geom_abline(intercept=0, slope=1)
filename=paste(my.wd, "/WinnersLosers paper/DCi/mean relabund correlation", ".pdf", sep=""); ggsave(filename, width=10, height=10)


#averaging again across treatments (including controls):

DCi.averaged2=DCi.averaged %>% 
  group_by(site_code, project_name, community_type, site_project_comm, site_project_comm_sp, species_matched) %>% 
  summarize(across(c(DCi, freq, mean.relabund, nplots), mean)) %>% 
  ungroup()

#how many site_project_comms is each species present in?

sp.widespread=DCi.averaged2 %>% 
  group_by(species_matched) %>% 
  summarize(number.site_project_comm=length(site_project_comm)) %>% 
  ungroup()
sp.widespread[sp.widespread$number.site_project_comm>20,]

ggplot(DCi.averaged2, aes(site_project_comm, DCi)) + geom_violin() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
filename=paste(my.wd, "/WinnersLosers paper/DCi/violin plot of species mean DCi across all treatment all years", ".pdf", sep=""); ggsave(filename, width=40, height=5)

ggplot(DCi.averaged2, aes(site_project_comm, freq)) + geom_violin() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
filename=paste(my.wd, "/WinnersLosers paper/DCi/violin plot of species mean frequencies across all treatment all years", ".pdf", sep=""); ggsave(filename, width=40, height=5)




