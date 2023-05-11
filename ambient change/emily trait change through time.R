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

###read in data

my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "E:/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "/Users/egrman/Library/CloudStorage/Dropbox/sDiv_sCoRRE_shared/"

#raw abundance data
dat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RawAbundance_Jan2023.csv",sep=""))

#relative abundance data
reldat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Jan2023.csv", sep=""))

#info on treatments, remove pre-treatment data, remove experiments where we have less than 6 years of species comp data, remove data in years 31+ from experiments longer than 30 years
trts<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfo_Dec2021.csv", sep="")) %>%
  filter(treatment_year>0) %>% 
  filter(treatment_year<31) %>%
  filter(successional==0) %>%
  filter(plant_mani==0) %>%
  group_by(site_code, project_name, community_type, treatment, trt_type) %>%
  summarize(expt_length=max(treatment_year), n_yrs_data=length(unique(calendar_year))) %>%
  filter(n_yrs_data>6)

#cleaned species names
sp <-read.csv(paste(my.wd,"CoRRE data/trait data/FullList_Nov2021.csv", sep=""))%>%
  select(genus_species, species_matched)%>%
  unique

#reading in categorical traits
my_cat=read.csv(paste(my.wd, "CoRRE data/trait data/sCoRRE categorical trait data_12142022.csv", sep="")) %>%
  select(species_matched, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal, n_fixation) %>%
  mutate(drop=ifelse(species_matched %in% c("Andreaea obovata", "Anthelia juratzkana", "Aulacomnium turgidum", "Barbilophozia hatcheri", "Barbilophozia kunzeana", "Blepharostoma trichophyllum", "Brachythecium albicans", "Bryum arcticum", "Bryum pseudotriquetrum", "Campylium stellatum", "Cyrtomnium hymenophyllum", "Dicranoweisia crispula", "Dicranum brevifolium", "Dicranum elongatum", "Dicranum fuscescens", "Dicranum groenlandicum",  "Dicranum scoparium", "Distichium capillaceum", "Ditrichum flexicaule", "Gymnomitrion concinnatum", "Hamatocaulis vernicosus", "Homalothecium pinnatifidum", "Hylocomium splendens", "Hypnum cupressiforme", "Hypnum hamulosum", "Isopterygiopsis pulchella", "Kiaeria starkei", "Leiocolea heterocolpos", "Marchantia polymorpha", "Marsupella brevissima", "Meesia uliginosa", "Myurella tenerrima", "Oncophorus virens", "Oncophorus wahlenbergii", "Pleurozium schreberi", "Pogonatum urnigerum", "Pohlia cruda", "Pohlia nutans", "Polytrichastrum alpinum", "Polytrichum juniperinum", "Polytrichum piliferum", "Polytrichum strictum", "Preissia quadrata", "Ptilidium ciliare", "Racomitrium lanuginosum", "Rhytidium rugosum", "Saelania glaucescens", "Sanionia uncinata",  "Schistidium apocarpum", "Syntrichia ruralis","Tomentypnum nitens", "Tortella tortuosa", "Tritomaria quinquedentata", "Nephroma arcticum", "Unknown NA", "Campylopus flexuosus", "Hypnum jutlandicum", "Plagiothecium undulatum", "Polytrichum commune", "Pseudoscleropodium purum", "Rhytidiadelphus loreus", "Rhytidiadelphus triquetrus", "Thuidium tamariscinum"), 1, 0)) %>% 
  filter(drop==0) 

#combine relative abundance data with treatment, cleaned species names
myreldat<-reldat%>%
  left_join(trts)%>%
  left_join(sp)%>% #this drops the unknowns??
  na.omit()%>% 
  mutate(drop=ifelse(species_matched %in% c("Andreaea obovata", "Anthelia juratzkana", "Aulacomnium turgidum", "Barbilophozia hatcheri", "Barbilophozia kunzeana", "Blepharostoma trichophyllum", "Brachythecium albicans", "Bryum arcticum", "Bryum pseudotriquetrum", "Campylium stellatum", "Cyrtomnium hymenophyllum", "Dicranoweisia crispula", "Dicranum brevifolium", "Dicranum elongatum", "Dicranum fuscescens", "Dicranum groenlandicum",  "Dicranum scoparium", "Distichium capillaceum", "Ditrichum flexicaule", "Gymnomitrion concinnatum", "Hamatocaulis vernicosus", "Homalothecium pinnatifidum", "Hylocomium splendens", "Hypnum cupressiforme", "Hypnum hamulosum", "Isopterygiopsis pulchella", "Kiaeria starkei", "Leiocolea heterocolpos", "Marchantia polymorpha", "Marsupella brevissima", "Meesia uliginosa", "Myurella tenerrima", "Oncophorus virens", "Oncophorus wahlenbergii", "Pleurozium schreberi", "Pogonatum urnigerum", "Pohlia cruda", "Pohlia nutans", "Polytrichastrum alpinum", "Polytrichum juniperinum", "Polytrichum piliferum", "Polytrichum strictum", "Preissia quadrata", "Ptilidium ciliare", "Racomitrium lanuginosum", "Rhytidium rugosum", "Saelania glaucescens", "Sanionia uncinata",  "Schistidium apocarpum", "Syntrichia ruralis","Tomentypnum nitens", "Tortella tortuosa", "Tritomaria quinquedentata", "Nephroma arcticum", "Unknown NA", "Campylopus flexuosus", "Hypnum jutlandicum", "Plagiothecium undulatum", "Polytrichum commune", "Pseudoscleropodium purum", "Rhytidiadelphus loreus", "Rhytidiadelphus triquetrus", "Thuidium tamariscinum"), 1, 0)) %>% 
  filter(drop==0) %>% 
  group_by(site_code, project_name, community_type, calendar_year, treatment_year, treatment, expt_length, block, plot_id, trt_type, species_matched) %>% 
  summarize(relcov=sum(relcov)) %>% 
  mutate(site_project_comm=as.factor(paste(site_code, project_name, community_type, sep="::")))

#adding in zeros for species that were absent from a plot
spc=unique(myreldat$site_project_comm)
myreldat_filled=NULL

for (j in 1:length(spc)) {
  dat.filled=myreldat[myreldat$site_project_comm==as.character(spc[j]),] %>% 
    select(site_code, project_name, community_type, calendar_year, treatment_year, treatment, expt_length, block, plot_id, relcov, trt_type, species_matched, site_project_comm) %>% 
    pivot_wider(names_from="species_matched", values_from="relcov", values_fill=0)
  dat.keep=dat.filled %>% 
    pivot_longer(!c("site_code", "project_name", "community_type", "calendar_year", "treatment_year", "treatment", "block", "plot_id", "trt_type", "site_project_comm", "expt_length"), names_to="species_matched", values_to="relcov")
  myreldat_filled=rbind(myreldat_filled, dat.keep)
}

#get average relative cover for each species in a treatment, over all plots
relave<-myreldat_filled%>%
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, species_matched, calendar_year, treatment_year, expt_length)%>%
  summarize(mean.relabund=mean(relcov)) %>% 
  ungroup()

#getting frequency of each plot type
myplots<-myreldat_filled%>%
  select(site_code, project_name, community_type, site_project_comm, treatment, block, trt_type, plot_id, calendar_year, treatment_year, expt_length)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, calendar_year, treatment_year, expt_length)%>%
  summarize(ntotplots=length(plot_id)) %>% 
  ungroup()

#getting number of plots of each type in which a species was present
freq<-myreldat_filled %>%
  select(site_code, project_name, community_type, site_project_comm, treatment, trt_type, species_matched, block, plot_id, calendar_year, treatment_year, expt_length, relcov) %>%
  unique() %>%
  filter(relcov>0) %>% 
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, species_matched, calendar_year, treatment_year, expt_length) %>%
  summarize(nplots=length(plot_id)) %>%
  left_join(myplots) %>%
  mutate(freq=nplots/ntotplots)

#calculate DCi per species
DCi.species.per.year<-relave %>%
  left_join(freq) %>%
  mutate(freq=replace_na(freq, 0)) %>%
  mutate(nplots=replace_na(nplots, 0)) %>%
  mutate(DCi=(mean.relabund+freq)/2) %>%
  select(site_code, project_name, community_type, site_project_comm, treatment, trt_type, species_matched, calendar_year, treatment_year, expt_length, mean.relabund, nplots, freq, DCi)

#summarize across trait groups (lumping annuals and biennials; selecting only the traits we want), arranging/selecting treatments we want
DCi.cat.per.year<-DCi.species.per.year %>%
  left_join(my_cat) %>% 
  mutate(lifespan=ifelse(lifespan=="annual", "ann.bien", ifelse(lifespan=="biennial", "ann.bien", lifespan))) %>%
  mutate(clonal=ifelse(clonal=="yes", "clonal", ifelse(clonal=="no", "nonclonal", clonal))) %>%
  mutate(mycorrhizal=ifelse(mycorrhizal=="yes", "mycorrhizal", ifelse(mycorrhizal=="no", "nonmycorrhizal", mycorrhizal))) %>%
  mutate(n_fixation=ifelse(n_fixation=="yes", "Nfixer", ifelse(n_fixation=="no", "nonNfixer", n_fixation))) %>%
  pivot_longer(growth_form:n_fixation, names_to="prop", values_to="trait") %>%
  mutate(trait=as.factor(trait)) %>%
  filter(trait %in% c("forb", "C3", "perennial", "clonal", "mycorrhizal", "ann.bien", "nonclonal", "graminoid", "C4", "nonmycorrhizal", "woody", "Nfixer", "nonNfixer", "vine")) %>%
  mutate(my_trt=ifelse(trt_type %in% c("control", "temp", "CO2", "N", "irr", "drought", "P", "mult_nutrient"), trt_type, ifelse(trt_type %in% c("drought*temp", "irr*temp", "N*P", "N*temp", "mult_nutrient*drought", "N*CO2", "CO2*temp", "drought*CO2*temp", "N*irr", "mult_nutrient*temp", "N*drought", "N*CO2*temp", "irr*CO2*temp", "N*irr*CO2*temp", "irr*CO2", "N*irr*CO2", "N*irr*temp", "N*P*temp", "mult_nutrient*irr"), "GCD", "other"))) %>%
  mutate(property=as.factor(ifelse(trait %in% c("forb", "woody", "vine", "graminoid"), "Growth form", ifelse(trait %in% c("C3", "C4"), "Photosynthetic pathway", ifelse(trait %in% c("perennial", "ann.bien"), "Life span", ifelse(trait %in% c("clonal", "nonclonal"), "Clonality", ifelse(trait %in% c("Nfixer", "nonNfixer"), "N fixation", ifelse(trait %in% c("mycorrhizal", "nonmycorrhizal"), "Mycorrhizal", my_trt)))))))) %>%
  mutate(property=factor(property, levels=c("Life span", "Photosynthetic pathway", "Clonality", "N fixation", "Mycorrhizal", "Growth form"))) %>%
  mutate(trait=factor(trait, levels=c("ann.bien", "perennial", "C3", "C4", "clonal", "nonclonal", "Nfixer", "nonNfixer", "mycorrhizal", "nonmycorrhizal", "graminoid", "forb", "woody", "vine"))) %>%
  filter(!my_trt=="other") %>%
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, my_trt, calendar_year, treatment_year, expt_length, trait, property) %>%
  summarize(mean.sp.DCi=mean(DCi), sum.sp.relabund=sum(mean.relabund))

#plotting functional group abundances through time:
ggplot(aes(treatment_year, sum.sp.relabund), data=DCi.cat.per.year[DCi.cat.per.year$trt_type=="control" & DCi.cat.per.year$property=="Photosynthetic pathway",]) + geom_point(aes(color=trait)) + facet_wrap(~site_project_comm, scales="free") + geom_smooth(method="lm", se=F, aes(color=trait))
ggsave(paste(my.wd, "ambient change paper/figs may 2023/C3 C4 sum relabund through time.pdf", sep=""), width=20, height=12)

ggplot(aes(treatment_year, mean.sp.DCi), data=DCi.cat.per.year[DCi.cat.per.year$trt_type=="control" & DCi.cat.per.year$property=="Photosynthetic pathway",]) + geom_point(aes(color=trait)) + facet_wrap(~site_project_comm, scales="free") + geom_smooth(method="lm", se=F, aes(color=trait))
ggsave(paste(my.wd, "ambient change paper/figs may 2023/C3 C4 mean DCi through time.pdf", sep=""), width=20, height=12)

#which metric to use? 
ggplot(aes(sum.sp.relabund, mean.sp.DCi), data=DCi.cat.per.year[DCi.cat.per.year$trt_type=="control" & DCi.cat.per.year$property=="Photosynthetic pathway",]) + geom_point(aes(color=trait)) + facet_wrap(~site_project_comm) + geom_smooth(method="lm", se=F, aes(color=trait))
ggsave(paste(my.wd, "ambient change paper/figs may 2023/C3 C4 comparing metrics through time.pdf", sep=""), width=20, height=16)
#sum.sp.relabund seems more sensitive, so using that one!

#need to go back and try Padu's advice from Dec 2022 about how to sum up cover


#get slopes of sum.sp.abund through time separately for each treatment
spc_list=as.character(unique(DCi.cat.per.year$site_project_comm))
change_over_time=numeric(0)

for(i in 1:length(spc_list)) {
  dati=DCi.cat.per.year[DCi.cat.per.year$site_project_comm==as.character(spc_list[i]),]
  trt_list=as.character(unique(dati$treatment))
  change_over_timei=numeric(0)
  
  for(j in 1:length(trt_list)) {
    datj=dati[dati$treatment==as.character(trt_list[j]),]
    trait_list=as.character(unique(datj$trait))
    change_over_timej=numeric(0)
    
    for(k in 1:length(trait_list)) {
      datk=datj[datj$trait==as.character(trait_list[k]),]
      reg=lm(sum.sp.relabund~treatment_year, data=datk)
      change_over_timek=data.frame(row.names=k, site_code=datk[1,c("site_code")], project_name=datk[1,c("project_name")], community_type=datk[1,c("community_type")], site_project_comm=as.character(spc_list[i]), expt_length=datk[1, "expt_length"], treatment=as.character(trt_list[j]), my_trt=datk[1,c("my_trt")], trait=datk[1, "trait"], property=datk[1, "property"], slope=reg$coefficients[2], R2=summary(reg)$r.squared)
      change_over_timej=rbind(change_over_timej, change_over_timek)
    }
    change_over_timei=rbind(change_over_timei, change_over_timej)
  }
  change_over_time=rbind(change_over_time, change_over_timei)
}

write.csv(change_over_time, paste(my.wd, "ambient change paper/slopes of summed sp abund.csv", sep=""), row.names=F)


#compare slope against study length to see if there is an obvious cutoff
ggplot(aes(expt_length, R2), data=change_over_time[change_over_time$my_trt=="control",]) + geom_point(alpha=I(0.1)) + facet_wrap(~trait)
ggsave(paste(my.wd, "ambient change paper/figs may 2023/slopes of summed species relative abundances vs treatment length, control plots.pdf", sep=""), width=7, height=7)

hist(change_over_time[change_over_time$my_trt=="control",]$expt_length)


# PLOT AMBIENT SLOPE VS TREATMENT SLOPE FOR EACH TREATMENT AND EACH TRAIT GROUP

#first, identify the sites with each treatment (so we know which control slopes to average)
#then calculate treatment-control slope for each spc

control_change_over_time <- change_over_time %>%
  filter(my_trt=="control") %>%
  mutate(control.slope=slope, abs.control.slope=abs(slope)) %>%
  select(site_code, project_name, community_type, site_project_comm, expt_length, trait, property, control.slope, abs.control.slope) %>%
  group_by(trait, property)

trt_change_over_time <- change_over_time %>%
  filter(!my_trt=="control") %>%
  mutate(trt.slope=slope, abs.trt.slope=abs(slope)) %>%
  select(site_code, project_name, community_type, site_project_comm, expt_length, treatment, my_trt, trait, property, trt.slope, abs.trt.slope) %>%
  left_join(control_change_over_time) %>%
  mutate(TminusC=trt.slope-control.slope, TminusCabs=abs.trt.slope-abs.control.slope) %>%
  filter(!is.na(control.slope)) %>%
  filter(!is.na(trt.slope)) %>%
  group_by(my_trt, trait, property) 


#do responses to drivers predict treatment effects? (or vice versa)

ggplot(aes(control.slope, TminusC, color=trait), data=trt_change_over_time) + geom_point(aes(shape=my_trt)) + facet_wrap(~property) + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + scale_shape_manual(values=c(8, 2, 16, 17, 10, 1, 5, 15)) + geom_smooth(method="lm", se=F) 
ggsave(paste(my.wd, "ambient change paper/figs may 2023/change in controls vs TminusC by FG.pdf", sep=""), width=7, height=7)
#can't do this because of mathematical correlation between the axes (change in control is in both)

#FIG 2:

ggplot(aes(control.slope, trt.slope, color=trait), data=trt_change_over_time) + geom_point(aes(shape=my_trt)) + facet_wrap(~property) + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + scale_shape_manual(values=c(8, 2, 16, 17, 10, 1, 5, 15)) + geom_smooth(method="lm", se=F) + geom_abline(intercept=0, slope=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/change in controls vs change in trt by FG.pdf", sep=""), width=7, height=7)
#need to do major axis regression or reduced major axis regression--they differ in some details but test deviation from 1:1 line rather than deviation from slope=0?

mean_trt_change_over_time <- trt_change_over_time %>%
  summarize(avgT=mean(trt.slope), sdT=sd(trt.slope), avgTminusC=mean(TminusC), sdTminusC=sd(TminusC), n=length(TminusC), avgC=mean(control.slope), sdC=sd(control.slope), avgTminusCabs=mean(TminusCabs), sdTminusCabs=sd(TminusCabs)) %>%
  mutate(TminusCcolor=ifelse(avgC>0 & avgT>0, "consistent increases", ifelse(avgC<0 & avgT<0, "consistent decreases", "disagreement")), seTminusCabs=sdTminusCabs/sqrt(n))

# we are trying to get a single metric that will communicate whether the treatment and controls have similar magnitudes and similar directions. but we are not happy with this figure:
ggplot(aes(trait, avgTminusCabs), data=mean_trt_change_over_time) + geom_point(aes(color=TminusCcolor)) + facet_wrap(~my_trt, scales="free_y") + geom_errorbar(aes(ymin=avgTminusCabs-seTminusCabs, ymax=avgTminusCabs+seTminusCabs, width=0.1, color=TminusCcolor)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + ylab("Difference between treatment and control\nin absolute value of change\nin relative abundance over time") + scale_color_manual(values=c("darkblue", "cornflowerblue", "orange"), name="Trait group changes\nin treated and control plots") + geom_hline(yintercept=0, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/trait groups responses, treatment minus control.pdf", sep=""), width=10, height=6)

#FIG 1:

mean_control_change_over_time <- control_change_over_time %>%
  summarize(globalavgC=mean(control.slope), globalsdC=sd(control.slope), n=length(control.slope)) %>%
  mutate(globalseC=globalsdC/sqrt(n), trt=as.factor("control")) %>%
  mutate(trt=factor(trt, levels=c("control", "treatment")))

fig1a=ggplot(aes(trait, globalavgC, color=trait), data=mean_control_change_over_time) + geom_point(aes(shape=trt)) + facet_wrap(~property, scales="free") + geom_errorbar(aes(ymin=globalavgC-globalseC, ymax=globalavgC+globalseC, width=0.1)) + theme(axis.text.x = element_blank()) + xlab("") + ylab("Change in relative abundance over time") + ggtitle("a) Control plots") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4"), name=NULL) + scale_shape_manual(values=c(1, 16), name=NULL, drop=F) + geom_hline(yintercept=0, color="black")

toplot.trt<-mean_trt_change_over_time %>%
  select(my_trt, trait, property, avgT, sdT, n) %>%
  mutate(se=sdT/sqrt(n), avg=avgT, trt="treatment") %>%
  select(my_trt, trait, property, avg, se, trt)
toplot.TminusC<-mean_trt_change_over_time %>%
  select(my_trt, trait, property, avgTminusC, sdTminusC, n) %>%
  mutate(se=sdTminusC/sqrt(n), avg=avgTminusC, trt="T-C") %>%
  select(my_trt, trait, property, avg, se, trt)
toplot.control<-mean_trt_change_over_time %>%
  select(my_trt, trait, property, avgC, sdC, n) %>%
  mutate(se=sdC/sqrt(n), avg=avgC, trt="control") %>%
  select(my_trt, trait, property, avg, se, trt)
toplot=rbind(toplot.trt, toplot.TminusC, toplot.control)

fig1b=ggplot(aes(my_trt, avg, color=trait), data=toplot[toplot$trt %in% c("control", "treatment"),]) + geom_point(aes(shape=trt)) + facet_wrap(~property, scales="free_y") + geom_errorbar(aes(ymin=avg-se, ymax=avg+se, width=0.1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") + xlab("") + ylab("Change in relative abundance over time") + ggtitle("b) Treatment plots relative to control plots in those experiments") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + scale_shape_manual(values=c(1, 16)) + geom_hline(yintercept=0, color="black")

#wow: http://www.sthda.com/english/articles/32-r-graphics-essentials/126-combine-multiple-ggplots-in-one-graph/
ggarrange(fig1a, fig1b, nrow=2)
ggsave(paste(my.wd, "ambient change paper/figs may 2023/trait groups responses, global controls and treatment vs control.pdf", sep=""), width=8, height=10)


# CAN WE RELATE CHANGE OVER TIME IN CONTROL PLOTS TO GLOBAL CHANGE DRIVERS?

#read in data from Adam's GIS person, subset to last 30 years
driverdata=read.csv(paste(my.wd, "CoRRE data/CoRRE data/environmental data/ScorreSitesTmaxTminPrecip 1901-2016.csv", sep="")) %>%
  filter(Year>1986) %>%
  group_by(site_code)

SummerTmax=driverdata %>%
  filter(Season=="Summer") %>%
  summarize(SummerTmax=mean(Tmax))

WinterTmin=driverdata %>%
  filter(Season=="Winter") %>%
  summarize(WinterTmin=mean(Tmin))

Annualprecip=driverdata %>%
  summarize(Annualprecip=sum(Precip), Ndep2016=mean(Ndep_2016))

drivers=Annualprecip %>%
  left_join(SummerTmax) %>%
  left_join(WinterTmin) %>%
  left_join(control_change_over_time, multiple="all")






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










####### OLD JUNK CODE ##################3
  





  
  pivot_wider(names_from="species_matched", values_from="relcov", values_fill=0)
  summarize(mean.sp.DCi=mean(DCi), sum.sp.relabund=sum(mean.relabund))
  
  group_by(site_code, project_name, community_type, site_project_comm, treatment, trt_type, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal, n_fixation, calendar_year, treatment_year) %>%










filename=(paste(my.wd, "/", sep="")); write.csv(DCi.through.time, filename, row.names=F)

#DCi.through.time=read.csv(paste(my.wd, "/", sep=""))

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




