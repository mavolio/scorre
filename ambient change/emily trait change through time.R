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

#info on treatments, remove pre-treatment data, remove experiments where we have less than 6 years of species comp data, remove treatments we don't want
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
  filter(trait %in% c("forb", "C3", "perennial", "clonal", "mycorrhizal", "ann.bien", "nonclonal", "graminoid", "C4", "nonmycorrhizal", "woody", "Nfixer", "nonNfixer", "vine")) %>%
  mutate(property=as.factor(ifelse(trait %in% c("forb", "woody", "vine", "graminoid"), "Growth form", ifelse(trait %in% c("C3", "C4"), "Photosynthetic pathway", ifelse(trait %in% c("perennial", "ann.bien"), "Life span", ifelse(trait %in% c("clonal", "nonclonal"), "Clonality", ifelse(trait %in% c("Nfixer", "nonNfixer"), "N fixation", ifelse(trait %in% c("mycorrhizal", "nonmycorrhizal"), "Mycorrhizal", my_trt)))))))) %>%
  mutate(property=factor(property, levels=c("Life span", "Photosynthetic pathway", "Clonality", "N fixation", "Mycorrhizal", "Growth form"))) %>%
  mutate(trait=factor(trait, levels=c("ann.bien", "perennial", "C3", "C4", "clonal", "nonclonal", "Nfixer", "nonNfixer", "mycorrhizal", "nonmycorrhizal", "graminoid", "forb", "woody", "vine"))) %>%
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

#need to go back and try Padu's advice from Dec 2022 about how to sum up cover


###   DO FUNCTIONAL GROUP ABUNDANCES CHANGE THROUGH TIME?   ###


control_change_over_time <- change_over_time %>%
  filter(my_trt=="control") %>%
  mutate(control.DCi.slope=DCi.slope, control.relabund.slope=relabund.slope, control.rawabund.slope=rawabund.slope) %>%
  select(site_code, project_name, community_type, site_project_comm, expt_length, trait, property, control.DCi.slope, control.relabund.slope, control.rawabund.slope) 

#FIG 1:

#do functional groups differ from each other?
#can look into nonparametric tests but for now:

propertylist=unique(control_change_over_time$property)
fg_responses=NULL

for(i in 1:length(propertylist)) {
  mod=lm(control.relabund.slope ~ trait, data=control_change_over_time[control_change_over_time$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), p=paste("p=", round(p$Pr[1], 2), sep=""))
  hist(resid(mod))
  qqPlot(resid(mod))
  fg_responses=rbind(fg_responses, temp)
}
# all of these models are singular when use site_code as random effect. p values are identical. which seems suspicious. 
# residuals don't look good either. need to try nonparametric tests. 

ggplot(aes(trait, control.relabund.slope, color=trait), data=control_change_over_time) + geom_boxplot() + facet_wrap(~property, scales="free") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4"), name=NULL) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_text(data=fg_responses, aes(label=p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses, global controls boxplots.pdf", sep=""), width=8, height=5)

#can we use 95% CI to test whether slopes differ from zero? need to see if they are normally distributed. probably not--see qqPlots. 

ggplot(aes(control.relabund.slope), data=control_change_over_time) + geom_histogram() + facet_wrap(~trait, scales="free")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/histograms of control relabund slopes.pdf", sep=""), width=10, height=6)
unique(control_change_over_time$trait)
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="ann.bien"])
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="perennial"])
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="C3"])
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="clonal"])
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="nonclonal"])
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="Nfixer"])
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="nonNfixer"])
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="mycorrhizal"])
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="nonmycorrhizal"])
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="graminoid"])
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="forb"])
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="vine"])
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="woody"])
qqPlot(control_change_over_time$control.relabund.slope[control_change_over_time$trait=="C4"])

mean_control_change_over_time <- control_change_over_time %>%
  group_by(trait, property) %>% 
  summarize(globalavgC=mean(control.relabund.slope), globalsdC=sd(control.relabund.slope), n=length(control.relabund.slope)) %>%
  mutate(globalseC=globalsdC/sqrt(n), trt=as.factor("control"), global95CIC=1.96*globalsdC/sqrt(n)) %>%
  mutate(trt=factor(trt, levels=c("control", "treatment"))) %>% 
  ungroup

fig1a=ggplot(aes(trait, globalavgC, color=trait), data=mean_control_change_over_time) + geom_point(aes(shape=trt)) + facet_wrap(~property, scales="free") + geom_errorbar(aes(ymin=globalavgC-global95CIC, ymax=globalavgC+global95CIC, width=0.1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + ylab("Change in relative abundance over time (+/- 95% CI)") + ggtitle("a) Control plots") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4"), name=NULL) + scale_shape_manual(values=c(1, 16), name=NULL, drop=F) + geom_hline(yintercept=0, color="black") + geom_text(data=fg_responses, aes(label=p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")


###   HOW DO CHANGES IN CONTROLS COMPARE TO CHANGES IN TREATED PLOTS?   ###

#first need to select only the controls for which we have corresponding treatment data

trt_change_over_time <- change_over_time %>%
  filter(!my_trt=="control") %>%
  mutate(trt.DCi.slope=DCi.slope, trt.relabund.slope=relabund.slope, trt.rawabund.slope=rawabund.slope) %>%
  select(site_code, project_name, community_type, site_project_comm, expt_length, treatment, my_trt, trait, property, trt.DCi.slope, trt.relabund.slope, trt.rawabund.slope) %>%
  left_join(control_change_over_time) %>%
  mutate(DCi.TminusC=trt.DCi.slope-control.DCi.slope, relabund.TminusC=trt.relabund.slope-control.relabund.slope, rawabund.TminusC=trt.rawabund.slope-control.rawabund.slope)

#do responses to drivers predict treatment effects? (or vice versa)

ggplot(aes(control.relabund.slope, relabund.TminusC, color=trait), data=trt_change_over_time) + geom_point(aes(shape=my_trt)) + facet_wrap(~property, scales="free") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + scale_shape_manual(values=c(8, 2, 16, 17, 10, 1, 5, 15)) + geom_smooth(method="lm", se=F) 
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund change in controls vs TminusC.pdf", sep=""), width=7, height=7)
#can't do this because of mathematical correlation between the axes (change in control is in both)

#FIG 2:

ggplot(aes(control.relabund.slope, trt.relabund.slope, color=trait), data=trt_change_over_time) + geom_point(aes(shape=my_trt)) + facet_wrap(~property, scales="free") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + scale_shape_manual(values=c(8, 2, 16, 17, 10, 1, 5, 15)) + geom_smooth(method="lm", se=F) + geom_abline(intercept=0, slope=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund change in controls vs change in trt.pdf", sep=""), width=10, height=7)
#need to do major axis regression or reduced major axis regression--they differ in some details but test deviation from 1:1 line rather than deviation from slope=0?

mean_trt_change_over_time <- trt_change_over_time %>%
  group_by(my_trt, trait, property) %>% 
  summarize(avg.relabundT=mean(trt.relabund.slope), sd.relabundT=sd(trt.relabund.slope), n=length(trt.relabund.slope), avg.relabundC=mean(control.relabund.slope), sd.relabundC=sd(control.relabund.slope)) %>% 
  ungroup() %>% 
  mutate(se.relabundT=sd.relabundT/sqrt(n), se.relabundC=sd.relabundC/sqrt(n))

toplot.trt<-mean_trt_change_over_time %>%
  mutate(avg=avg.relabundT, se=se.relabundT, trt="treatment") %>%
  select(my_trt, trait, property, avg, se, trt)
toplot.control<-mean_trt_change_over_time %>%
  mutate(avg=avg.relabundC, se=se.relabundC, trt="control") %>%
  select(my_trt, trait, property, avg, se, trt)
toplot=rbind(toplot.trt, toplot.control)

fig1b=ggplot(aes(my_trt, avg, color=trait), data=toplot) + geom_point(aes(shape=trt)) + facet_wrap(~property, scales="free") + geom_errorbar(aes(ymin=avg-se, ymax=avg+se, width=0.1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") + xlab("") + ylab("Change in relative abundance over time") + ggtitle("b) Treatment plots relative to control plots in those experiments") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + scale_shape_manual(values=c(1, 16)) + geom_hline(yintercept=0, color="black")

#wow: http://www.sthda.com/english/articles/32-r-graphics-essentials/126-combine-multiple-ggplots-in-one-graph/
ggarrange(fig1a, fig1b, nrow=2)
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund, global controls and treatment vs control.pdf", sep=""), width=6, height=10)


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


###   DO FUNCTIONAL GROUPS IN THE CONTROL PLOTS RESPOND DIFFERENTLY TO THE DRIVERS OVER THE TIMESCALES OF THE EXPERIMENTS?   ###

#can't avoid 2-step analysis (first extracting slopes, then seeing what they correlate with) because we only have N deposition at one point in time. 

drivers_by_year <- DCi.cat.per.year %>% 
  left_join(drivers, join_by(site_code, calendar_year==Year), multiple="all") %>% 
  filter(calendar_year<2017) #need to drop these for now because we don't have environmental data past 2016

#extracting slopes over the timeperiods of the study, carrying Ndep_2016 through

#for environmental variables by site_project_comm (because some experiments at a site run for different durations than others):

spclist=unique(drivers_by_year$site_project_comm)
global_change=NULL

for(i in 1:length(spclist)) {
  dati=drivers_by_year[drivers_by_year$site_project_comm==as.character(spclist[i]),]
  tempi=NULL
  trtlist=unique(dati$treatment)
  
  for(j in 1:length(trtlist)) {
    datj=dati[dati$treatment==as.character(trtlist[j]),]
    summerreg=lm(SummerTmax~calendar_year, data=dati)
    winterreg=lm(WinterTmin~calendar_year, data=dati)
    annprecip=lm(AnnualPrecip~calendar_year, data=dati)
    co2reg=lm(ppm.CO2~calendar_year, data=dati)
    tempj=data.frame(row.names=i, site_project_comm=as.character(spclist[i]), treatment=as.character(trtlist[j]), change_in_SummerTmax=summerreg$coefficients[2], change_in_WinterTmin=winterreg$coefficients[2], change_in_AnnualPrecip=annprecip$coefficients[2], change_in_CO2=co2reg$coefficients[2], Ndep_2016=mean(dati$Ndep_2016, na.rm=T))
    tempi=rbind(tempi, tempj)
  }
  global_change=rbind(global_change, tempi)
}

global_change <- global_change %>% 
  left_join(trts)

#for functional groups, by site_project_comm then treatment then trait:

spclist=as.character(unique(drivers_by_year$site_project_comm))
change_over_time=NULL

for(i in 1:length(spc_list)) {
  dati=drivers_by_year[drivers_by_year$site_project_comm==as.character(spclist[i]),]
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
  left_join(global_change) 

write.csv(change_over_time, paste(my.wd, "ambient change paper/slopes of 3 metrics of FG abund change.csv", sep=""), row.names=F)

#compare slope against study length to see if there is an obvious cutoff or any trends
ggplot(aes(expt_length, relabund.slope), data=change_over_time[change_over_time$my_trt=="control",]) + geom_point(alpha=I(0.2)) + facet_wrap(~trait, scales="free", ncol=5)
ggsave(paste(my.wd, "ambient change paper/figs may 2023/slopes of summed species relative abundances vs treatment length, control plots.pdf", sep=""), width=9, height=5)

hist(change_over_time[change_over_time$my_trt=="control",]$expt_length)
min(change_over_time$expt_length)

#DCi without CO2:

propertylist=unique(change_over_time$property)
fg_DCi=NULL

for(i in 1:length(propertylist)) {
  mod=lmer(DCi.slope ~ trait + change_in_SummerTmax + change_in_WinterTmin + change_in_AnnualPrecip + Ndep_2016 + trait:change_in_SummerTmax + trait:change_in_WinterTmin + trait:change_in_AnnualPrecip + trait:Ndep_2016 + (1|site_code), data=change_over_time[change_over_time$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), SummerTmax.p=paste("p=", round(p$Pr[6], 3), sep=""), WinterTmin.p=paste("p=", round(p$Pr[7], 3), sep=""), AnnualPrecip.p=paste("p=", round(p$Pr[8], 3), sep=""), Ndep.p=paste("p=", round(p$Pr[9], 3), sep=""))
  fg_DCi=rbind(fg_DCi, temp)
}
#all of these models are singular. try removing random effect?

ggplot(aes(change_in_SummerTmax, DCi.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_DCi, aes(label=SummerTmax.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG DCi responses to changes in SummerTmax.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_WinterTmin, DCi.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_DCi, aes(label=WinterTmin.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG DCi responses to changes in WinterTmin.pdf", sep=""), width=8, height=5)

ggplot(aes(Ndep_2016, DCi.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_DCi, aes(label=Ndep.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG DCi responses to N deposition.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_AnnualPrecip, DCi.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_DCi, aes(label=AnnualPrecip.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG DCi responses to changes in AnnualPrecip.pdf", sep=""), width=8, height=5)

# sum relabund without CO2

propertylist=unique(change_over_time$property)
fg_relabund=NULL

for(i in 1:length(propertylist)) {
  mod=lmer(relabund.slope ~ trait + change_in_SummerTmax + change_in_WinterTmin + change_in_AnnualPrecip + Ndep_2016 + trait:change_in_SummerTmax + trait:change_in_WinterTmin + trait:change_in_AnnualPrecip + trait:Ndep_2016 + (1|site_code), data=change_over_time[change_over_time$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), SummerTmax.p=paste("p=", round(p$Pr[6], 3), sep=""), WinterTmin.p=paste("p=", round(p$Pr[7], 3), sep=""), AnnualPrecip.p=paste("p=", round(p$Pr[8], 3), sep=""), Ndep.p=paste("p=", round(p$Pr[9], 3), sep=""))
  fg_relabund=rbind(fg_relabund, temp)
}
#all of these models are singular. try removing random effect?

ggplot(aes(change_in_SummerTmax, relabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_relabund, aes(label=SummerTmax.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in SummerTmax.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_WinterTmin, relabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_relabund, aes(label=WinterTmin.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in WinterTmin.pdf", sep=""), width=8, height=5)

ggplot(aes(Ndep_2016, relabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_relabund, aes(label=Ndep.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to N deposition.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_AnnualPrecip, relabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_relabund, aes(label=AnnualPrecip.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in AnnualPrecip.pdf", sep=""), width=8, height=5)

# summed raw abundance without CO2

propertylist=unique(change_over_time$property)
fg_rawabund=NULL

for(i in 1:length(propertylist)) {
  mod=lmer(rawabund.slope ~ trait + change_in_SummerTmax + change_in_WinterTmin + change_in_AnnualPrecip + Ndep_2016 + trait:change_in_SummerTmax + trait:change_in_WinterTmin + trait:change_in_AnnualPrecip + trait:Ndep_2016 + (1|site_code), data=change_over_time[change_over_time$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), SummerTmax.p=paste("p=", round(p$Pr[6], 3), sep=""), WinterTmin.p=paste("p=", round(p$Pr[7], 3), sep=""), AnnualPrecip.p=paste("p=", round(p$Pr[8], 3), sep=""), Ndep.p=paste("p=", round(p$Pr[9], 3), sep=""))
  fg_rawabund=rbind(fg_rawabund, temp)
}
#all of these models are singular. try removing random effect?

ggplot(aes(change_in_SummerTmax, rawabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_rawabund, aes(label=SummerTmax.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG rawabund responses to changes in SummerTmax.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_WinterTmin, rawabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_rawabund, aes(label=WinterTmin.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG rawabund responses to changes in WinterTmin.pdf", sep=""), width=8, height=5)

ggplot(aes(Ndep_2016, rawabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_rawabund, aes(label=Ndep.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG rawabund responses to N deposition.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_AnnualPrecip, rawabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_rawabund, aes(label=AnnualPrecip.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG rawabund responses to changes in AnnualPrecip.pdf", sep=""), width=8, height=5)

#DCi with CO2:

propertylist=unique(change_over_time$property)
fg_DCi=NULL

for(i in 1:length(propertylist)) {
  mod=lmer(DCi.slope ~ trait + change_in_SummerTmax + change_in_WinterTmin + change_in_AnnualPrecip + change_in_CO2 + Ndep_2016 + trait:change_in_SummerTmax + trait:change_in_WinterTmin + trait:change_in_AnnualPrecip + trait:change_in_CO2 + trait:Ndep_2016 + (1|site_code), data=change_over_time[change_over_time$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), SummerTmax.p=paste("p=", round(p$Pr[7], 3), sep=""), WinterTmin.p=paste("p=", round(p$Pr[8], 3), sep=""), AnnualPrecip.p=paste("p=", round(p$Pr[9], 3), sep=""), CO2.p=paste("p=", round(p$Pr[10], 3), sep=""), Ndep.p=paste("p=", round(p$Pr[11], 3), sep=""))
  fg_DCi=rbind(fg_DCi, temp)
}
#all of these models are singular. try removing random effect?

ggplot(aes(change_in_SummerTmax, DCi.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_DCi, aes(label=SummerTmax.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG DCi responses to changes in SummerTmax.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_WinterTmin, DCi.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_DCi, aes(label=WinterTmin.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG DCi responses to changes in WinterTmin.pdf", sep=""), width=8, height=5)

ggplot(aes(Ndep_2016, DCi.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_DCi, aes(label=Ndep.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG DCi responses to N deposition.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_AnnualPrecip, DCi.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_DCi, aes(label=AnnualPrecip.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG DCi responses to changes in AnnualPrecip.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_CO2, DCi.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_DCi, aes(label=CO2.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG DCi responses to changes in CO2.pdf", sep=""), width=8, height=5)

# sum relabund with CO2

propertylist=unique(change_over_time$property)
fg_relabund=NULL

for(i in 1:length(propertylist)) {
  mod=lmer(relabund.slope ~ trait + change_in_SummerTmax + change_in_WinterTmin + change_in_AnnualPrecip + change_in_CO2 + Ndep_2016 + trait:change_in_SummerTmax + trait:change_in_WinterTmin + trait:change_in_AnnualPrecip + trait:change_in_CO2 + trait:Ndep_2016 + (1|site_code), data=change_over_time[change_over_time$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), SummerTmax.p=paste("p=", round(p$Pr[7], 3), sep=""), WinterTmin.p=paste("p=", round(p$Pr[8], 3), sep=""), AnnualPrecip.p=paste("p=", round(p$Pr[9], 3), sep=""), CO2.p=paste("p=", round(p$Pr[10], 3), sep=""), Ndep.p=paste("p=", round(p$Pr[11], 3), sep=""))
  fg_relabund=rbind(fg_relabund, temp)
}
#all of these models are singular. try removing random effect?

ggplot(aes(change_in_SummerTmax, relabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_relabund, aes(label=SummerTmax.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in SummerTmax.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_WinterTmin, relabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_relabund, aes(label=WinterTmin.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in WinterTmin.pdf", sep=""), width=8, height=5)

ggplot(aes(Ndep_2016, relabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_relabund, aes(label=Ndep.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to N deposition.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_AnnualPrecip, relabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_relabund, aes(label=AnnualPrecip.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in AnnualPrecip.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_CO2, relabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_relabund, aes(label=CO2.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG relabund responses to changes in CO2.pdf", sep=""), width=8, height=5)

# summed raw abundance with CO2

propertylist=unique(change_over_time$property)
fg_rawabund=NULL

for(i in 1:length(propertylist)) {
  mod=lmer(rawabund.slope ~ trait + change_in_SummerTmax + change_in_WinterTmin + change_in_AnnualPrecip + change_in_CO2 + Ndep_2016 + trait:change_in_SummerTmax + trait:change_in_WinterTmin + trait:change_in_AnnualPrecip + trait:change_in_CO2 + trait:Ndep_2016 + (1|site_code), data=change_over_time[change_over_time$property==as.character(propertylist[i]),]); p=Anova(mod, test.statistic="F")
  temp=data.frame(row.names=i, property=as.factor(propertylist[i]), SummerTmax.p=paste("p=", round(p$Pr[7], 3), sep=""), WinterTmin.p=paste("p=", round(p$Pr[8], 3), sep=""), AnnualPrecip.p=paste("p=", round(p$Pr[9], 3), sep=""), CO2.p=paste("p=", round(p$Pr[10], 3), sep=""), Ndep.p=paste("p=", round(p$Pr[11], 3), sep=""))
  fg_rawabund=rbind(fg_rawabund, temp)
}
#all of these models are singular. try removing random effect?

ggplot(aes(change_in_SummerTmax, rawabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_rawabund, aes(label=SummerTmax.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG rawabund responses to changes in SummerTmax.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_WinterTmin, rawabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_rawabund, aes(label=WinterTmin.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG rawabund responses to changes in WinterTmin.pdf", sep=""), width=8, height=5)

ggplot(aes(Ndep_2016, rawabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_rawabund, aes(label=Ndep.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG rawabund responses to N deposition.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_AnnualPrecip, rawabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_rawabund, aes(label=AnnualPrecip.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG rawabund responses to changes in AnnualPrecip.pdf", sep=""), width=8, height=5)

ggplot(aes(change_in_CO2, rawabund.slope, color=trait), data=change_over_time) + geom_point(shape=1) + facet_wrap(~property, scale="free") + geom_smooth(method="lm") + scale_color_manual(values=c("darksalmon", "darkred", "orange", "darkorange3", "gold", "darkgoldenrod2", "greenyellow", "green4", "dodgerblue", "dodgerblue4", "thistle", "plum", "orchid4", "darkorchid4")) + geom_text(data=fg_rawabund, aes(label=CO2.p, x=Inf, y=Inf), vjust=1.5, hjust=1, color="black")
ggsave(paste(my.wd, "ambient change paper/figs may 2023/FG rawabund responses to changes in CO2.pdf", sep=""), width=8, height=5)













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
