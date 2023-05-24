#### Emily Grman modified emily DCi trends through time.R
#### Dec 13, 2022 sCoRRE meeting 4 at iDiv
#### May 10, 2023 sCoRRE meeting 4.1 at UNCG


rm(list=ls())

library(tidyverse)
library(gridExtra)
library(lme4)
library(car)
library(ggpubr)

#emily's plotting preferences:
theme_set(theme_bw())
theme_eg=theme_update(panel.grid.minor=element_blank(), strip.background=element_blank())

###   READ IN TREATMENT, SPECIES NAMES, TRAIT DATA    ###

my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "E:/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "/Users/egrman/Library/CloudStorage/Dropbox/sDiv_sCoRRE_shared/"

#info on treatments, remove pre-treatment data, remove experiments in successional environments or where seeds were added, remove experiments where we have less than 6 years of species comp data, remove treatments we don't want
trts<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfo_Dec2021.csv", sep="")) %>%
  mutate(site_project_comm=as.factor(paste(site_code, project_name, community_type, sep="::"))) %>% 
  filter(treatment_year>0) %>% 
  filter(successional==0) %>%
  filter(plant_mani==0) %>%
  mutate(my_trt=ifelse(trt_type %in% c("control", "temp", "CO2", "N", "irr", "drought", "P", "mult_nutrient"), trt_type, ifelse(trt_type %in% c("drought*temp", "irr*temp", "N*P", "N*temp", "mult_nutrient*drought", "N*CO2", "CO2*temp", "drought*CO2*temp", "N*irr", "mult_nutrient*temp", "N*drought", "N*CO2*temp", "irr*CO2*temp", "N*irr*CO2*temp", "irr*CO2", "N*irr*CO2", "N*irr*temp", "N*P*temp", "mult_nutrient*irr"), "GCD", "other"))) %>%
  filter(!my_trt=="other") %>%
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, my_trt) %>%
  summarize(expt_length=max(treatment_year), n_yrs_data=length(unique(calendar_year)), start_year=min(calendar_year), end_year=max(calendar_year)) %>%
  mutate(year_span=end_year-start_year+1) %>% 
  ungroup() %>%
  filter(n_yrs_data>6) %>% 
  mutate(start_year_minus5=start_year-5, start_year_minus10=start_year-10)

#cleaned species names
sp <-read.csv(paste(my.wd,"CoRRE data/trait data/FullList_Nov2021.csv", sep=""))%>%
  select(genus_species, species_matched)%>%
  unique

#reading in categorical traits
my_cat=read.csv(paste(my.wd, "CoRRE data/trait data/sCoRRE categorical trait data_12142022.csv", sep="")) %>%
  select(species_matched, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal, n_fixation) %>%
  mutate(drop=ifelse(species_matched %in% c("Andreaea obovata", "Anthelia juratzkana", "Aulacomnium turgidum", "Barbilophozia hatcheri", "Barbilophozia kunzeana", "Blepharostoma trichophyllum", "Brachythecium albicans", "Bryum arcticum", "Bryum pseudotriquetrum", "Campylium stellatum", "Cyrtomnium hymenophyllum", "Dicranoweisia crispula", "Dicranum brevifolium", "Dicranum elongatum", "Dicranum fuscescens", "Dicranum groenlandicum",  "Dicranum scoparium", "Distichium capillaceum", "Ditrichum flexicaule", "Gymnomitrion concinnatum", "Hamatocaulis vernicosus", "Homalothecium pinnatifidum", "Hylocomium splendens", "Hypnum cupressiforme", "Hypnum hamulosum", "Isopterygiopsis pulchella", "Kiaeria starkei", "Leiocolea heterocolpos", "Marchantia polymorpha", "Marsupella brevissima", "Meesia uliginosa", "Myurella tenerrima", "Oncophorus virens", "Oncophorus wahlenbergii", "Pleurozium schreberi", "Pogonatum urnigerum", "Pohlia cruda", "Pohlia nutans", "Polytrichastrum alpinum", "Polytrichum juniperinum", "Polytrichum piliferum", "Polytrichum strictum", "Preissia quadrata", "Ptilidium ciliare", "Racomitrium lanuginosum", "Rhytidium rugosum", "Saelania glaucescens", "Sanionia uncinata",  "Schistidium apocarpum", "Syntrichia ruralis","Tomentypnum nitens", "Tortella tortuosa", "Tritomaria quinquedentata", "Nephroma arcticum", "Unknown NA", "Campylopus flexuosus", "Hypnum jutlandicum", "Plagiothecium undulatum", "Polytrichum commune", "Pseudoscleropodium purum", "Rhytidiadelphus loreus", "Rhytidiadelphus triquetrus", "Thuidium tamariscinum"), 1, 0)) %>% 
  filter(drop==0) 


###   READ IN SPECIES RAW AND RELATIVE ABUNDANCE DATA, CALCULATE MEANS AND DCI   ###


#raw abundance data
rawdat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RawAbundance_Jan2023.csv",sep=""))

#relative abundance data
reldat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Jan2023.csv", sep=""))

#combine relative and raw abundance data with treatment, cleaned species names
mydat<-reldat%>%
  left_join(rawdat) %>% 
  left_join(trts) %>%
  left_join(sp) %>% #this drops the unknowns??
  na.omit() %>% 
  mutate(drop=ifelse(species_matched %in% c("Andreaea obovata", "Anthelia juratzkana", "Aulacomnium turgidum", "Barbilophozia hatcheri", "Barbilophozia kunzeana", "Blepharostoma trichophyllum", "Brachythecium albicans", "Bryum arcticum", "Bryum pseudotriquetrum", "Campylium stellatum", "Cyrtomnium hymenophyllum", "Dicranoweisia crispula", "Dicranum brevifolium", "Dicranum elongatum", "Dicranum fuscescens", "Dicranum groenlandicum",  "Dicranum scoparium", "Distichium capillaceum", "Ditrichum flexicaule", "Gymnomitrion concinnatum", "Hamatocaulis vernicosus", "Homalothecium pinnatifidum", "Hylocomium splendens", "Hypnum cupressiforme", "Hypnum hamulosum", "Isopterygiopsis pulchella", "Kiaeria starkei", "Leiocolea heterocolpos", "Marchantia polymorpha", "Marsupella brevissima", "Meesia uliginosa", "Myurella tenerrima", "Oncophorus virens", "Oncophorus wahlenbergii", "Pleurozium schreberi", "Pogonatum urnigerum", "Pohlia cruda", "Pohlia nutans", "Polytrichastrum alpinum", "Polytrichum juniperinum", "Polytrichum piliferum", "Polytrichum strictum", "Preissia quadrata", "Ptilidium ciliare", "Racomitrium lanuginosum", "Rhytidium rugosum", "Saelania glaucescens", "Sanionia uncinata",  "Schistidium apocarpum", "Syntrichia ruralis","Tomentypnum nitens", "Tortella tortuosa", "Tritomaria quinquedentata", "Nephroma arcticum", "Unknown NA", "Campylopus flexuosus", "Hypnum jutlandicum", "Plagiothecium undulatum", "Polytrichum commune", "Pseudoscleropodium purum", "Rhytidiadelphus loreus", "Rhytidiadelphus triquetrus", "Thuidium tamariscinum"), 1, 0)) %>% 
  filter(drop==0) %>% 
  group_by(site_code, project_name, community_type, site_project_comm, calendar_year, treatment_year, treatment, expt_length, block, plot_id, trt_type, my_trt, species_matched) %>% 
  summarize(relcov=sum(relcov), abundance=sum(abundance)) %>% 
  ungroup() 

#adding in zeros for species that were absent from a plot
spc=unique(mydat$site_project_comm)
myreldat_filled=NULL

for (j in 1:length(spc)) {
  dat.filled=mydat[mydat$site_project_comm==as.character(spc[j]),] %>% 
    select(site_code, project_name, community_type, calendar_year, treatment_year, treatment, expt_length, block, plot_id, relcov, trt_type, my_trt, species_matched, site_project_comm) %>% 
    pivot_wider(names_from="species_matched", values_from="relcov", values_fill=0)
  dat.keep=dat.filled %>% 
    pivot_longer(!c("site_code", "project_name", "community_type", "calendar_year", "treatment_year", "treatment", "block", "plot_id", "trt_type", "my_trt", "site_project_comm", "expt_length"), names_to="species_matched", values_to="relcov")
  myreldat_filled=rbind(myreldat_filled, dat.keep)
}

#adding raw abundance measures, filling zeros for species absent from plots
mydat_filled <- myreldat_filled %>% 
  left_join(mydat) %>% 
  mutate(abundance=replace_na(abundance, 0))
  
#averaging over replicate plots and calculating DCi:
#get average relative cover for each species in a treatment, over all plots
relave<-mydat_filled %>%
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, my_trt, species_matched, calendar_year, treatment_year, expt_length) %>%
  summarize(mean.relabund=mean(relcov), mean.rawabund=mean(abundance)) %>% 
  ungroup()
#getting frequency of each plot type
myplots<-mydat_filled %>%
  select(site_code, project_name, community_type, site_project_comm, treatment, block, trt_type, my_trt, plot_id, calendar_year, treatment_year, expt_length) %>%
  unique() %>%
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, my_trt, calendar_year, treatment_year, expt_length) %>%
  summarize(ntotplots=length(plot_id)) %>% 
  ungroup()
#getting number of plots of each type in which a species was present
freq<-mydat_filled %>%
  select(site_code, project_name, community_type, site_project_comm, treatment, trt_type, my_trt, species_matched, block, plot_id, calendar_year, treatment_year, expt_length, relcov) %>%
  unique() %>%
  filter(relcov>0) %>% 
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, my_trt, species_matched, calendar_year, treatment_year, expt_length) %>%
  summarize(nplots=length(plot_id)) %>%
  left_join(myplots) %>%
  mutate(freq=nplots/ntotplots)
#calculate DCi per species
DCi.species.per.year<-relave %>%
  left_join(freq) %>%
  mutate(freq=replace_na(freq, 0)) %>%
  mutate(nplots=replace_na(nplots, 0)) %>%
  mutate(DCi=(mean.relabund+freq)/2) %>%
  select(site_code, project_name, community_type, site_project_comm, treatment, trt_type, my_trt, species_matched, calendar_year, treatment_year, expt_length, mean.relabund, mean.rawabund, nplots, freq, DCi)

#summarize across trait groups (lumping annuals and biennials; selecting only the traits we want)
DCi.cat.per.year<-DCi.species.per.year %>%
  left_join(my_cat) %>% 
  mutate(lifespan=ifelse(lifespan=="annual", "ann.bien", ifelse(lifespan=="biennial", "ann.bien", lifespan))) %>%
  mutate(clonal=ifelse(clonal=="yes", "clonal", ifelse(clonal=="no", "nonclonal", clonal))) %>%
  mutate(mycorrhizal=ifelse(mycorrhizal=="yes", "mycorrhizal", ifelse(mycorrhizal=="no", "nonmycorrhizal", mycorrhizal))) %>%
  mutate(n_fixation=ifelse(n_fixation=="yes", "Nfixer", ifelse(n_fixation=="no", "nonNfixer", n_fixation))) %>%
  pivot_longer(growth_form:n_fixation, names_to="prop", values_to="trait") %>%
  mutate(trait=as.factor(trait)) %>%
  filter(trait %in% c("forb", "C3", "perennial", "clonal", "mycorrhizal", "ann.bien", "nonclonal", "graminoid", "C4", "nonmycorrhizal", "Nfixer", "nonNfixer")) %>%
  mutate(property=as.factor(ifelse(trait %in% c("forb", "graminoid"), "Growth form", ifelse(trait %in% c("C3", "C4"), "Photosynthetic pathway", ifelse(trait %in% c("perennial", "ann.bien"), "Life span", ifelse(trait %in% c("clonal", "nonclonal"), "Clonality", ifelse(trait %in% c("Nfixer", "nonNfixer"), "N fixation", ifelse(trait %in% c("mycorrhizal", "nonmycorrhizal"), "Mycorrhizal", my_trt)))))))) %>%
  mutate(property=factor(property, levels=c("Life span", "Photosynthetic pathway", "Clonality", "N fixation", "Mycorrhizal", "Growth form"))) %>%
  mutate(trait=factor(trait, levels=c("ann.bien", "perennial", "C3", "C4", "clonal", "nonclonal", "Nfixer", "nonNfixer", "mycorrhizal", "nonmycorrhizal", "graminoid", "forb"))) %>%
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, my_trt, calendar_year, treatment_year, expt_length, trait, property) %>%
  summarize(mean.sp.DCi=mean(DCi), sum.sp.relabund=sum(mean.relabund), sum.sp.rawabund=sum(mean.rawabund)) %>% 
  ungroup()

#need to go back and try Padu's advice from Dec 2022 about how to sum up cover

#plotting functional group abundances through time:
ggplot(aes(calendar_year, sum.sp.relabund), data=DCi.cat.per.year[DCi.cat.per.year$trt_type=="control" & DCi.cat.per.year$property=="Life span",]) + geom_point(aes(color=trait)) + facet_wrap(~site_project_comm, scales="free") + geom_smooth(method="lm", se=F, aes(color=trait))
ggsave(paste(my.wd, "ambient change paper/figs may 2023/lifespan sum relabund through time.pdf", sep=""), width=20, height=12)

ggplot(aes(calendar_year, mean.sp.DCi), data=DCi.cat.per.year[DCi.cat.per.year$trt_type=="control" & DCi.cat.per.year$property=="Life span",]) + geom_point(aes(color=trait)) + facet_wrap(~site_project_comm, scales="free") + geom_smooth(method="lm", se=F, aes(color=trait))
ggsave(paste(my.wd, "ambient change paper/figs may 2023/lifespan mean DCi through time.pdf", sep=""), width=20, height=12)

ggplot(aes(calendar_year, sum.sp.rawabund), data=DCi.cat.per.year[DCi.cat.per.year$trt_type=="control" & DCi.cat.per.year$property=="Life span",]) + geom_point(aes(color=trait)) + facet_wrap(~site_project_comm, scales="free") + geom_smooth(method="lm", se=F, aes(color=trait))
ggsave(paste(my.wd, "ambient change paper/figs may 2023/lifespan sum rawabund through time.pdf", sep=""), width=20, height=12)



###   CALCULATING CHANGE OVER TIME IN FUNCTIONAL GROUPS   ###


#by site_project_comm then treatment then trait:

spclist=as.character(unique(DCi.cat.per.year$site_project_comm))
change_over_time=NULL

for(i in 1:length(spclist)) {
  dati=DCi.cat.per.year[DCi.cat.per.year$site_project_comm==as.character(spclist[i]) & DCi.cat.per.year$calendar_year<2017,] #need to drop these for now because we don't have environmental data past 2016
  change_over_timei=NULL
  trtlist=unique(dati$treatment)
  
  for(j in 1:length(trtlist)) {
    datj=dati[dati$treatment==as.character(trtlist[j]),]
    change_over_timej=NULL
    traitlist=as.character(unique(datj$trait))
    
    for(k in 1:length(traitlist)) {
      datk=datj[datj$trait==as.character(traitlist[k]),]
      DCireg=lm(mean.sp.DCi~calendar_year, data=datk)
      relabundreg=lm(sum.sp.relabund~calendar_year, data=datk)
      rawabundreg=lm(sum.sp.rawabund~calendar_year, data=datk)
      change_over_timek=data.frame(row.names=k, site_code=datk[1,"site_code"], site_project_comm=datk[1,"site_project_comm"], property=datk[1, "property"], treatment=datk[1, "treatment"], trait=datk[1, "trait"], DCi.slope=DCireg$coefficients[2], relabund.slope=relabundreg$coefficients[2], rawabund.slope=rawabundreg$coefficients[2])
      change_over_timej=rbind(change_over_timej, change_over_timek)
    }
    change_over_timei=rbind(change_over_timei, change_over_timej)
  }
  change_over_time=rbind(change_over_time, change_over_timei)
}
change_over_time <- change_over_time %>% 
  left_join(trts) 

write.csv(change_over_time, paste(my.wd, "ambient change paper/slopes of 3 metrics of FG abund change.csv", sep=""), row.names=F)
#compare slope against study length to see if there is an obvious cutoff or any trends
hist(change_over_time[change_over_time$my_trt=="control",]$expt_length)
min(change_over_time$expt_length)
ggplot(aes(expt_length, relabund.slope), data=change_over_time[change_over_time$my_trt=="control",]) + geom_point(alpha=I(0.2)) + facet_wrap(~trait, scales="free", ncol=5)
ggsave(paste(my.wd, "ambient change paper/figs may 2023/slopes of summed species relative abundances vs treatment length, control plots.pdf", sep=""), width=9, height=5)


###   DO FUNCTIONAL GROUP ABUNDANCES CHANGE THROUGH TIME?   ###   averaging across replicate studies at a site regardless of their length


sitemean_control_change_over_time <- change_over_time %>%
  filter(my_trt=="control") %>%
  group_by(site_code, trait, property) %>% 
  summarize(control.DCi.slope=mean(DCi.slope), control.relabund.slope=mean(relabund.slope), control.rawabund.slope=mean(rawabund.slope)) %>%
  ungroup() %>% 
  select(site_code, trait, property, control.DCi.slope, control.relabund.slope, control.rawabund.slope) 

#FIG 1:

#do functional groups differ from each other?
#can look into nonparametric tests but for now:

propertylist=unique(sitemean_control_change_over_time$property)
fg_responses=NULL

for(i in 1:length(propertylist)) {
  mod=lm(control.relabund.slope ~ trait, data=sitemean_control_change_over_time[sitemean_control_change_over_time$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), p=paste("p=", round(p$Pr[1], 2), sep=""))
  hist(resid(mod), main=as.character(propertylist[i]))
  qqPlot(resid(mod), main=as.character(propertylist[i]))
  fg_responses=rbind(fg_responses, temp)
}

# can only use 95% CI to test whether slopes differ from zero if data are normally distributed. i think these are better enough--see qqPlots. 
# can still try nonparametric tests. 

ggplot(aes(trait, control.relabund.slope, color=trait), data=sitemean_control_change_over_time) + geom_boxplot() + facet_wrap(~property, scales="free") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4"), name=NULL) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_text(data=fg_responses, aes(label=p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black") + ggtitle("Control plots, averaged across all experiments at a site (regardless of study duration)")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses, global controls boxplots.pdf", sep=""), width=8, height=5)

global_control_change_over_time <- sitemean_control_change_over_time %>%
  group_by(trait, property) %>% 
  summarize(globalavgC=mean(control.relabund.slope), globalsdC=sd(control.relabund.slope), n=length(control.relabund.slope)) %>%
  mutate(globalseC=globalsdC/sqrt(n), trt=as.factor("control"), global95CIC=1.96*globalsdC/sqrt(n)) %>%
  mutate(trt=factor(trt, levels=c("control", "treatment"))) %>% 
  ungroup()

fig1a=ggplot(aes(trait, globalavgC, color=trait), data=global_control_change_over_time) + geom_point(aes(shape=trt)) + facet_wrap(~property, scales="free") + geom_errorbar(aes(ymin=globalavgC-global95CIC, ymax=globalavgC+global95CIC, width=0.1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + ylab("Change in relative abundance over time (+/- 95% CI)") + ggtitle("a) Control plots") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4"), name=NULL) + scale_shape_manual(values=c(1, 16), name=NULL, drop=F) + geom_hline(yintercept=0, color="black") + geom_text(data=fg_responses, aes(label=p, x=Inf, y=Inf), vjust=1.5, hjust=1.1, color="black") + scale_y_continuous(expand = expansion(mult=0.2))


# + scale_y_continuous(expand = expansion(mult = c(0, 0.2)))


###   DO RESPONSES TO DRIVERS PREDICT TREATMENT EFFECTS? ARE CHANGEY SITES (IN CONTROLS) MORE CHANGEY (IN TREATMENTS)?   ###


#first need to select only the controls for which we have corresponding treatment data (at the site_project_comm level)
control_change_over_time <- change_over_time %>%
  filter(my_trt=="control") %>%
  mutate(control.DCi.slope=DCi.slope, control.relabund.slope=relabund.slope, control.rawabund.slope=rawabund.slope) %>%
  select(site_code, site_project_comm, expt_length, start_year, end_year, start_year_minus5, start_year_minus10, trait, property, control.DCi.slope, control.relabund.slope, control.rawabund.slope) 

trt_change_over_time <- change_over_time %>%
  filter(!my_trt=="control") %>%
  mutate(trt.DCi.slope=DCi.slope, trt.relabund.slope=relabund.slope, trt.rawabund.slope=rawabund.slope) %>%
  select(site_code, site_project_comm, expt_length, start_year, end_year, start_year_minus5, start_year_minus10, my_trt, trait, property, trt.DCi.slope, trt.relabund.slope, trt.rawabund.slope) %>%
  left_join(control_change_over_time, multiple="all") %>%
  mutate(DCi.TminusC=trt.DCi.slope-control.DCi.slope, relabund.TminusC=trt.relabund.slope-control.relabund.slope, rawabund.TminusC=trt.rawabund.slope-control.rawabund.slope)

#then averaging across experiments with the same trt at a site
mean_trt_change_over_time <- trt_change_over_time %>% 
  group_by(site_code, my_trt, trait, property) %>% 
  summarize(control.DCi.slope=mean(control.DCi.slope), control.relabund.slope=mean(control.relabund.slope), control.rawabund.slope=mean(control.rawabund.slope), trt.DCi.slope=mean(trt.DCi.slope), trt.relabund.slope=mean(trt.relabund.slope), trt.rawabund.slope=mean(trt.rawabund.slope)) %>% 
  ungroup()

#FIG 2:

ggplot(aes(control.relabund.slope, trt.relabund.slope, color=trait), data=mean_trt_change_over_time) + geom_point(aes(shape=my_trt)) + facet_wrap(~property, scales="free") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + scale_shape_manual(values=c(8, 2, 16, 17, 10, 1, 5, 15)) + geom_smooth(method="lm", se=F) + geom_abline(intercept=0, slope=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund change in controls vs change in trt.pdf", sep=""), width=10, height=7)
#need to do major axis regression or reduced major axis regression--they differ in some details but test deviation from 1:1 line rather than deviation from slope=0?
#adam says orthogonal regression because assumes error in both axes and therefore tries to minimize distance from point to line in both directions and not just y


###   CAN WE INFER WHAT GCD MIGHT BE DRIVING CHANGE OVER TIME IN CONTROLS BY LEARNING FROM FG RESPONSES IN EXPERIMENTS?   ### (relabund changes only)


#first averaging the matched trt-control data above across sites to get global trt responses and global controls (that match those experiments)
global_trt_change_over_time <- mean_trt_change_over_time %>%
  group_by(my_trt, trait, property) %>% 
  summarize(avg.relabundT=mean(trt.relabund.slope), sd.relabundT=sd(trt.relabund.slope), n=length(trt.relabund.slope), avg.relabundC=mean(control.relabund.slope), sd.relabundC=sd(control.relabund.slope)) %>% 
  ungroup() %>% 
  mutate(CI.relabundT=1.96*sd.relabundT/sqrt(n), CI.relabundC=1.96*sd.relabundC/sqrt(n))

#rearrange for plotting:
toplot.trt<-global_trt_change_over_time %>%
  mutate(avg=avg.relabundT, CI=CI.relabundT, trt="treatment") %>%
  select(my_trt, trait, property, avg, CI, trt)
toplot.control<-global_trt_change_over_time %>%
  mutate(avg=avg.relabundC, CI=CI.relabundC, trt="control") %>%
  select(my_trt, trait, property, avg, CI, trt)
toplot=rbind(toplot.trt, toplot.control)

fig1b=ggplot(aes(my_trt, avg, color=trait), data=toplot) + geom_point(aes(shape=trt)) + facet_wrap(~property, scales="free") + geom_errorbar(aes(ymin=avg-CI, ymax=avg+CI, width=0.1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") + xlab("") + ylab("Change in relative abundance over time (+/- 95% CI)") + ggtitle("b) Treatment plots relative to control plots in those experiments") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + scale_shape_manual(values=c(1, 16)) + geom_hline(yintercept=0, color="black")

#wow: http://www.sthda.com/english/articles/32-r-graphics-essentials/126-combine-multiple-ggplots-in-one-graph/
ggarrange(fig1a, fig1b, nrow=2)
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund, global controls and treatment vs control.pdf", sep=""), width=7.5, height=12)



###   GLOBAL CHANGE DRIVERS   ###


#read in CO2 data (downloaded from https://scrippsco2.ucsd.edu/data/atmospheric_co2/mlo.html on 11 may 2023) and calculate yearly mean
co2<-read.csv(paste(my.wd,"CoRRE data/CoRRE data/environmental data/monthly_in_situ_co2_mlo_EmilyforR.csv", sep="")) %>%
  filter(ppm.CO2>0) %>% #they have missing data as -99
  group_by(Year) %>% 
  summarize(ppm.CO2=mean(ppm.CO2)) %>% 
  ungroup()

ggplot(aes(Year, ppm.CO2), data=co2[co2$Year>1980,]) + geom_point() + geom_smooth(method="lm"); ggsave(paste(my.wd, "ambient change paper/figs may 2023/global CO2 over our years.pdf", sep=""), width=8, height=10)

#read in data from Adam's GIS person and calculate the summary variables we want (SummerTmax, WinterTmin, AnnualPrecip, Ndep in 2016)
driverdata=read.csv(paste(my.wd, "CoRRE data/CoRRE data/environmental data/ScorreSitesTmaxTminPrecip 1901-2016.csv", sep="")) %>%
  group_by(site_code, Year)
SummerTmax=driverdata %>%
  filter(Season=="Summer") %>%
  summarize(SummerTmax=mean(Tmax))
WinterTmin=driverdata %>%
  filter(Season=="Winter") %>%
  summarize(WinterTmin=mean(Tmin))
AnnualPrecip=driverdata %>%
  summarize(AnnualPrecip=sum(Precip), Ndep_2016=mean(Ndep_2016))
drivers=AnnualPrecip %>%
  left_join(SummerTmax) %>%
  left_join(WinterTmin) %>% 
  left_join(co2) %>% 
  ungroup()


###   CALCULATING CHANGE OVER TIME IN ENVIRONMENTAL DRIVERS 3 WAYS   ###

#1) for each site_project_comm, get environmental data from the years for which we have veg data (indicated with start and end columns in trts)
drivers_by_vegyear <- trts %>% 
  filter(my_trt=="control") %>% 
  group_by(site_project_comm) %>% 
  mutate(vegdate=map2(start_year, end_year, seq, by=1)) %>% 
  unnest(cols=c(vegdate)) %>% 
  select(site_code, site_project_comm, my_trt, vegdate) %>% 
  left_join(drivers, join_by(site_code, vegdate==Year)) %>% 
  ungroup()

#2) then do the same for the expanded years (5 years before and 10 years before)
drivers_by_lag5year <- trts %>% 
  filter(my_trt=="control") %>% 
  group_by(site_project_comm) %>% 
  mutate(lag5date=map2(start_year_minus5, end_year, seq, by=1)) %>% 
  unnest(cols=c(lag5date)) %>% 
  select(site_code, site_project_comm, my_trt, lag5date) %>% 
  left_join(drivers, join_by(site_code, lag5date==Year)) %>% 
  ungroup()

drivers_by_lag10year <- trts %>% 
  filter(my_trt=="control") %>% 
  group_by(site_project_comm) %>% 
  mutate(lag10date=map2(start_year_minus10, end_year, seq, by=1)) %>% 
  unnest(cols=c(lag10date)) %>% 
  select(site_code, site_project_comm, my_trt, lag10date)  %>% 
  left_join(drivers, join_by(site_code, lag10date==Year)) %>% 
  ungroup()

#3) calculate change over time for each timespan (3 for each site_project_comm), carrying Ndep_2016 through
spclist=unique(drivers_by_year$site_project_comm)
global_change=NULL

for(i in 1:length(spclist)) {
  dat.vegyear=drivers_by_vegyear[drivers_by_vegyear$site_project_comm==as.character(spclist[i]),]
  summerreg.vegyear=lm(SummerTmax~vegdate, data=dat.vegyear)
  winterreg.vegyear=lm(WinterTmin~vegdate, data=dat.vegyear)
  annprecip.vegyear=lm(AnnualPrecip~vegdate, data=dat.vegyear)
  
  dat.lag5year=drivers_by_lag5year[drivers_by_lag5year$site_project_comm==as.character(spclist[i]),]
  summerreg.lag5year=lm(SummerTmax~lag5date, data=dat.lag5year)
  winterreg.lag5year=lm(WinterTmin~lag5date, data=dat.lag5year)
  annprecip.lag5year=lm(AnnualPrecip~lag5date, data=dat.lag5year)
  
  dat.lag10year=drivers_by_lag10year[drivers_by_lag10year$site_project_comm==as.character(spclist[i]),]
  summerreg.lag10year=lm(SummerTmax~lag10date, data=dat.lag10year)
  winterreg.lag10year=lm(WinterTmin~lag10date, data=dat.lag10year)
  annprecip.lag10year=lm(AnnualPrecip~lag10date, data=dat.lag10year)
  
  tempi=data.frame(row.names=i, site_project_comm=as.character(spclist[i]), change_in_SummerTmax.vegyear=summerreg.vegyear$coefficients[2], change_in_WinterTmin.vegyear=winterreg.vegyear$coefficients[2], change_in_AnnualPrecip.vegyear=annprecip.vegyear$coefficients[2], change_in_SummerTmax.lag5year=summerreg.lag5year$coefficients[2], change_in_WinterTmin.lag5year=winterreg.lag5year$coefficients[2], change_in_AnnualPrecip.lag5year=annprecip.lag5year$coefficients[2], change_in_SummerTmax.lag10year=summerreg.lag10year$coefficients[2], change_in_WinterTmin.lag10year=winterreg.lag10year$coefficients[2], change_in_AnnualPrecip.lag10year=annprecip.lag10year$coefficients[2], Ndep_2016=mean(dat.vegyear$Ndep_2016, na.rm=T))

  global_change=rbind(global_change, tempi)
}

#4) merge with functional group changes over time 
global_change_in_fg <- change_over_time %>% 
  filter(my_trt=="control") %>% 
  select(site_code, site_project_comm, property, trait, DCi.slope, relabund.slope, rawabund.slope, expt_length) %>% 
  left_join(global_change) 

#5) trying that with averages across all experiments at a site (ignoring differences in study duration)
mean_global_change_in_fg <- global_change_in_fg %>% 
  group_by(site_code, property, trait) %>% 
  summarize(nstudies=length(DCi.slope), DCi.slope=mean(DCi.slope), relabund.slope=mean(relabund.slope), rawabund.slope=mean(rawabund.slope), change_in_SummerTmax.vegyear=mean(change_in_SummerTmax.vegyear), change_in_WinterTmin.vegyear=mean(change_in_WinterTmin.vegyear), change_in_AnnualPrecip.vegyear=mean(change_in_AnnualPrecip.vegyear), change_in_SummerTmax.lag5year=mean(change_in_SummerTmax.lag5year), change_in_WinterTmin.lag5year=mean(change_in_WinterTmin.lag5year), change_in_AnnualPrecip.lag5year=mean(change_in_AnnualPrecip.lag5year), change_in_SummerTmax.lag10year=mean(change_in_SummerTmax.lag10year), change_in_WinterTmin.lag10year=mean(change_in_WinterTmin.lag10year), change_in_AnnualPrecip.lag10year=mean(change_in_AnnualPrecip.lag10year), Ndep_2016=mean(Ndep_2016)) %>% 
  ungroup()
  

###   DO FUNCTIONAL GROUPS IN THE CONTROL PLOTS RESPOND DIFFERENTLY TO THE DRIVERS OVER THE TIMESCALES OF THE EXPERIMENTS?   ###  sum relabund only


#first, checking for correlations among our predictor variables:

pairs(global_change_in_fg[,c("change_in_SummerTmax.vegyear", "change_in_WinterTmin.vegyear", "change_in_AnnualPrecip.vegyear", "Ndep_2016")]) #2 measures of temp are correlated
pairs(global_change_in_fg[,c("change_in_SummerTmax.lag5year", "change_in_WinterTmin.lag5year", "change_in_AnnualPrecip.lag5year", "Ndep_2016")])
pairs(global_change_in_fg[,c("change_in_SummerTmax.lag10year", "change_in_WinterTmin.lag10year", "change_in_AnnualPrecip.lag10year", "Ndep_2016")])

#how much difference does it make to include the different lags?
pairs(global_change_in_fg[,c("change_in_SummerTmax.vegyear", "change_in_SummerTmax.lag5year", "change_in_SummerTmax.lag10year")])
pairs(global_change_in_fg[,c("change_in_WinterTmin.vegyear", "change_in_WinterTmin.lag5year", "change_in_WinterTmin.lag10year")])
pairs(global_change_in_fg[,c("change_in_AnnualPrecip.vegyear", "change_in_AnnualPrecip.lag5year", "change_in_AnnualPrecip.lag10year")])
#5 and 10 year tend to be correlated, but much weaker correlations between lagged datasets and the unlagged. 


# no lag, no averages

propertylist=unique(global_change_in_fg$property)
fg_relabund.vegdate=NULL

for(i in 1:length(propertylist)) {
  mod=lm(relabund.slope ~ trait + change_in_SummerTmax.vegyear + change_in_WinterTmin.vegyear + change_in_AnnualPrecip.vegyear + Ndep_2016 + trait:change_in_SummerTmax.vegyear + trait:change_in_WinterTmin.vegyear + trait:change_in_AnnualPrecip.vegyear + trait:Ndep_2016, data=global_change_in_fg[global_change_in_fg$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), SummerTmax.p=paste("p=", round(p$Pr[6], 3), sep=""), WinterTmin.p=paste("p=", round(p$Pr[7], 3), sep=""), AnnualPrecip.p=paste("p=", round(p$Pr[8], 3), sep=""), Ndep.p=paste("p=", round(p$Pr[9], 3), sep=""))
  hist(resid(mod), main=as.character(propertylist[i]))
  qqPlot(resid(mod), main=as.character(propertylist[i]))
  fg_relabund.vegdate=rbind(fg_relabund.vegdate, temp)
}

#qqplots look pretty bad

ggplot(aes(change_in_SummerTmax.vegyear, relabund.slope, color=trait), data=global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.vegdate, aes(label=SummerTmax.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in SummerTmax no lag.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_WinterTmin.vegyear, relabund.slope, color=trait), data=global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.vegdate, aes(label=WinterTmin.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in WinterTmin no lag.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_AnnualPrecip.vegyear, relabund.slope, color=trait), data=global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.vegdate, aes(label=AnnualPrecip.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in AnnualPrecip no lag.pdf", sep=""), width=8, height=5)

ggplot(aes(Ndep_2016, relabund.slope, color=trait), data=global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.vegdate, aes(label=Ndep.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to N deposition no lag.pdf", sep=""), width=8, height=5)


# 5 year lag, no averages

propertylist=unique(global_change_in_fg$property)
fg_relabund.lag5date=NULL

for(i in 1:length(propertylist)) {
  mod=lm(relabund.slope ~ trait + change_in_SummerTmax.lag5year + change_in_WinterTmin.lag5year + change_in_AnnualPrecip.lag5year + Ndep_2016 + trait:change_in_SummerTmax.lag5year + trait:change_in_WinterTmin.lag5year + trait:change_in_AnnualPrecip.lag5year + trait:Ndep_2016, data=global_change_in_fg[global_change_in_fg$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), SummerTmax.p=paste("p=", round(p$Pr[6], 3), sep=""), WinterTmin.p=paste("p=", round(p$Pr[7], 3), sep=""), AnnualPrecip.p=paste("p=", round(p$Pr[8], 3), sep=""), Ndep.p=paste("p=", round(p$Pr[9], 3), sep=""))
  hist(resid(mod), main=as.character(propertylist[i]))
  qqPlot(resid(mod), main=as.character(propertylist[i]))
  fg_relabund.lag5date=rbind(fg_relabund.lag5date, temp)
}

#qqplots still pretty bad

ggplot(aes(change_in_SummerTmax.lag5year, relabund.slope, color=trait), data=global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag5date, aes(label=SummerTmax.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in SummerTmax 5 year lag.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_WinterTmin.lag5year, relabund.slope, color=trait), data=global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag5date, aes(label=WinterTmin.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in WinterTmin 5 year lag.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_AnnualPrecip.lag5year, relabund.slope, color=trait), data=global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag5date, aes(label=AnnualPrecip.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in AnnualPrecip 5 year lag.pdf", sep=""), width=8, height=5)

ggplot(aes(Ndep_2016, relabund.slope, color=trait), data=global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag5date, aes(label=Ndep.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to N deposition 5 year lag.pdf", sep=""), width=8, height=5)



# 10 year lag, no averages

propertylist=unique(global_change_in_fg$property)
fg_relabund.lag10date=NULL

for(i in 1:length(propertylist)) {
  mod=lm(relabund.slope ~ trait + change_in_SummerTmax.lag10year + change_in_WinterTmin.lag10year + change_in_AnnualPrecip.lag10year + Ndep_2016 + trait:change_in_SummerTmax.lag10year + trait:change_in_WinterTmin.lag10year + trait:change_in_AnnualPrecip.lag10year + trait:Ndep_2016, data=global_change_in_fg[global_change_in_fg$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), SummerTmax.p=paste("p=", round(p$Pr[6], 3), sep=""), WinterTmin.p=paste("p=", round(p$Pr[7], 3), sep=""), AnnualPrecip.p=paste("p=", round(p$Pr[8], 3), sep=""), Ndep.p=paste("p=", round(p$Pr[9], 3), sep=""))
  hist(resid(mod), main=as.character(propertylist[i]))
  qqPlot(resid(mod), main=as.character(propertylist[i]))
  fg_relabund.lag10date=rbind(fg_relabund.lag10date, temp)
}

#qqplots still pretty bad

ggplot(aes(change_in_SummerTmax.lag10year, relabund.slope, color=trait), data=global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag10date, aes(label=SummerTmax.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in SummerTmax 10 year lag.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_WinterTmin.lag10year, relabund.slope, color=trait), data=global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag10date, aes(label=WinterTmin.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in WinterTmin 10 year lag.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_AnnualPrecip.lag10year, relabund.slope, color=trait), data=global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag10date, aes(label=AnnualPrecip.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in AnnualPrecip 10 year lag.pdf", sep=""), width=8, height=5)

ggplot(aes(Ndep_2016, relabund.slope, color=trait), data=global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag10date, aes(label=Ndep.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to N deposition 10 year lag.pdf", sep=""), width=8, height=5)



# no lag, with site averages

propertylist=unique(mean_global_change_in_fg$property)
fg_relabund.vegdate=NULL

for(i in 1:length(propertylist)) {
  mod=lm(relabund.slope ~ trait + change_in_SummerTmax.vegyear + change_in_WinterTmin.vegyear + change_in_AnnualPrecip.vegyear + Ndep_2016 + trait:change_in_SummerTmax.vegyear + trait:change_in_WinterTmin.vegyear + trait:change_in_AnnualPrecip.vegyear + trait:Ndep_2016, data=mean_global_change_in_fg[mean_global_change_in_fg$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), SummerTmax.p=paste("p=", round(p$Pr[6], 3), sep=""), WinterTmin.p=paste("p=", round(p$Pr[7], 3), sep=""), AnnualPrecip.p=paste("p=", round(p$Pr[8], 3), sep=""), Ndep.p=paste("p=", round(p$Pr[9], 3), sep=""))
  hist(resid(mod), main=as.character(propertylist[i]))
  qqPlot(resid(mod), main=as.character(propertylist[i]))
  fg_relabund.vegdate=rbind(fg_relabund.vegdate, temp)
}

#qqplots look better

ggplot(aes(change_in_SummerTmax.vegyear, relabund.slope, color=trait), data=mean_global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.vegdate, aes(label=SummerTmax.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in SummerTmax no lag, site averages.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_WinterTmin.vegyear, relabund.slope, color=trait), data=mean_global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.vegdate, aes(label=WinterTmin.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in WinterTmin no lag, site averages.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_AnnualPrecip.vegyear, relabund.slope, color=trait), data=mean_global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.vegdate, aes(label=AnnualPrecip.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in AnnualPrecip no lag, site averages.pdf", sep=""), width=8, height=5)

ggplot(aes(Ndep_2016, relabund.slope, color=trait), data=mean_global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.vegdate, aes(label=Ndep.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to N deposition no lag, site averages.pdf", sep=""), width=8, height=5)


# 5 year lag, with site averages

propertylist=unique(mean_global_change_in_fg$property)
fg_relabund.lag5date=NULL

for(i in 1:length(propertylist)) {
  mod=lm(relabund.slope ~ trait + change_in_SummerTmax.lag5year + change_in_WinterTmin.lag5year + change_in_AnnualPrecip.lag5year + Ndep_2016 + trait:change_in_SummerTmax.lag5year + trait:change_in_WinterTmin.lag5year + trait:change_in_AnnualPrecip.lag5year + trait:Ndep_2016, data=mean_global_change_in_fg[mean_global_change_in_fg$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), SummerTmax.p=paste("p=", round(p$Pr[6], 3), sep=""), WinterTmin.p=paste("p=", round(p$Pr[7], 3), sep=""), AnnualPrecip.p=paste("p=", round(p$Pr[8], 3), sep=""), Ndep.p=paste("p=", round(p$Pr[9], 3), sep=""))
  hist(resid(mod), main=as.character(propertylist[i]))
  qqPlot(resid(mod), main=as.character(propertylist[i]))
  fg_relabund.lag5date=rbind(fg_relabund.lag5date, temp)
}

#qqplots maybe OK?

ggplot(aes(change_in_SummerTmax.lag5year, relabund.slope, color=trait), data=mean_global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag5date, aes(label=SummerTmax.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in SummerTmax 5 year lag, site averages.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_WinterTmin.lag5year, relabund.slope, color=trait), data=mean_global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag5date, aes(label=WinterTmin.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in WinterTmin 5 year lag, site averages.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_AnnualPrecip.lag5year, relabund.slope, color=trait), data=mean_global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag5date, aes(label=AnnualPrecip.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in AnnualPrecip 5 year lag, site averages.pdf", sep=""), width=8, height=5)

ggplot(aes(Ndep_2016, relabund.slope, color=trait), data=mean_global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag5date, aes(label=Ndep.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to N deposition 5 year lag, site averages.pdf", sep=""), width=8, height=5)



# 10 year lag, with site averages

propertylist=unique(mean_global_change_in_fg$property)
fg_relabund.lag10date=NULL

for(i in 1:length(propertylist)) {
  mod=lm(relabund.slope ~ trait + change_in_SummerTmax.lag10year + change_in_WinterTmin.lag10year + change_in_AnnualPrecip.lag10year + Ndep_2016 + trait:change_in_SummerTmax.lag10year + trait:change_in_WinterTmin.lag10year + trait:change_in_AnnualPrecip.lag10year + trait:Ndep_2016, data=mean_global_change_in_fg[mean_global_change_in_fg$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), SummerTmax.p=paste("p=", round(p$Pr[6], 3), sep=""), WinterTmin.p=paste("p=", round(p$Pr[7], 3), sep=""), AnnualPrecip.p=paste("p=", round(p$Pr[8], 3), sep=""), Ndep.p=paste("p=", round(p$Pr[9], 3), sep=""))
  hist(resid(mod), main=as.character(propertylist[i]))
  qqPlot(resid(mod), main=as.character(propertylist[i]))
  fg_relabund.lag10date=rbind(fg_relabund.lag10date, temp)
}

#qqplots maybe OK?

ggplot(aes(change_in_SummerTmax.lag10year, relabund.slope, color=trait), data=mean_global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag10date, aes(label=SummerTmax.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in SummerTmax 10 year lag, site averages.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_WinterTmin.lag10year, relabund.slope, color=trait), data=mean_global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag10date, aes(label=WinterTmin.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in WinterTmin 10 year lag, site averages.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_AnnualPrecip.lag10year, relabund.slope, color=trait), data=mean_global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag10date, aes(label=AnnualPrecip.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in AnnualPrecip 10 year lag, site averages.pdf", sep=""), width=8, height=5)

ggplot(aes(Ndep_2016, relabund.slope, color=trait), data=mean_global_change_in_fg) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "plum", "orchid4")) + geom_text(data=fg_relabund.lag10date, aes(label=Ndep.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to N deposition 10 year lag, site averages.pdf", sep=""), width=8, height=5)






















# DOES SITE-LEVEL N DEPOSITION PREDICT HOW STRONGLY N FIXERS DECLINE OVER TIME?

ndep=read.csv(paste(my.wd, "ambient change paper/CoRRE_siteLocationClimateNdep_Dec2021.csv", sep="")) %>%
  left_join(change_over_time) %>%
  filter(my_trt=="control" & trait %in% c("Nfixer", "nonNfixer"))

qplot(ndep_2016, slope, data=ndep, alpha=log(expt_length), color=trait) + geom_smooth(method="lm") + scale_color_manual(values=c("forestgreen", "black")) + guides(alpha=FALSE) + xlab("N deposition in 2016") + ylab("Change in relative abundance\n over time in control plots")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/N deposition N fixers.pdf", sep=""), width=5, height=3)

try=lmer(slope~trait*ndep_2016 + (1|site_code/site_project_comm), data=ndep)
Anova(try) #interaction p=0.002
tryN=lmer(slope~ndep_2016 + (1|site_code), data=ndep[ndep$trait=="Nfixer",]) #ndep p=0.096
trynonN=lmer(slope~ndep_2016 + (1|site_code), data=ndep[ndep$trait=="nonNfixer",]) #ndep p=0.02










#junk




mutate(vegdate=map(data, ~seq(.x$start_year, .x$end_year, by="1 year")), lag5date=map(data, ~seq(.x$start_year_minus5, .x$end_year, by="1 year")), lag10date=map(data, ~seq(.x$start_year_minus10, .x$end_year, by="1 year"))) %>% 
  
  
  
  
  
  #there do not appear to be too many duplicate control plots with the same study length
  ggplot(aes(expt_length, fill=my_trt), data=change_over_time[change_over_time$trait=="C3" & change_over_time$my_trt=="control",]) + geom_histogram() + facet_wrap(~site_code, scales="free")
#this means it doesn't make sense to average across experiments at a site, unless we also want to somehow composite site change over time for all the experiments of different durations
#but averaging across replicate studies with the same my_trt and same year span at a site
mean_change_over_time <- change_over_time %>% 
  group_by(site_code, property, trait, my_trt, expt_length, n_yrs_data, start_year, end_year, year_span, start_year_minus5, start_year_minus10) %>% 
  summarize(n.studies=length(DCi.slope), DCi.slope=mean(DCi.slope), relabund.slope=mean(relabund.slope), rawabund.slope=mean(rawabund.slope)) %>% 
  ungroup()
#so where are all those points stacked up at same x value coming from?




drivers_by_year1 <- change_over_time %>% 
  group_by(site_code, expt_length, start_year, start_year_minus5, start_year_minus10, end_year, property, trait)
left_join(drivers, join_by(site_code, calendar_year==Year), multiple="all") %>% 
  filter(calendar_year<2017) #need to drop these for now because we don't have environmental data past 2016



