### Looking at mathematical patterns among DCi, abundance, and frequency
###
### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Created March 8, 2021; last updated March 8, 2021

####investigate winners losers

library(tidyverse)
library(gridExtra)
rm(list=ls())
my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "C:/Users/wilco/Dropbox/shared working groups/sDiv_sCoRRE_shared/"

#raw abundance data
dat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RawAbundanceMar2021.csv",sep=""))

#relative abundance data
reldat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RelativeAbundanceMar2021.csv",sep=""))

#info on treatments
trts<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfoMar2021.csv", sep=""))%>%
  select(-X)

#cleaned species names
sp <-read.csv(paste(my.wd,"CoRRE data/CoRRE data/trait data/CoRRE2trykey.csv", sep=""))%>%
  select(genus_species, species_matched)%>%
  unique

#combine relative abundance data with treatment and cleaned species names
allreldat<-reldat%>%
  left_join(trts)%>%
  left_join(sp)%>%# this drops the unknowns
  na.omit()

#get average realtive cover for each species in a treatment, over all years of the experiment over all plots
relave <- allreldat%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani, species_matched) %>%
  summarize(mean=mean(relcov))

#subset out control plots, plot_mani==0
Crelave<-relave%>%
  filter(plot_mani==0)%>%
  ungroup()%>%
  select(-treatment, -plot_mani)

#subset out treated plots, plot_mani!=0
Trelave<-relave%>%
  filter(plot_mani!=0)%>%
  ungroup()%>%
  select(-plot_mani)

#to get relative frequency, determine number of control plots
controlplots<-allreldat%>%
  filter(plot_mani==0)%>%
  select(site_code, project_name, community_type, plot_id)%>%
  unique()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(ncplots=length(plot_id))

#to get relative frequency, determine number of control plots a species is found in, merge in total number of plots and calculate relative frequency  
control_freq<-allreldat%>%
  filter(plot_mani==0)%>%
  select(site_code, project_name, community_type, species_matched, plot_id)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, species_matched)%>%
  summarize(nplots=length(plot_id))%>%
  left_join(controlplots)%>%
  mutate(freq=nplots/ncplots)

#calculate DCi
control_dom<-control_freq%>%
  left_join(Crelave)%>%
  mutate(DCi=(mean+freq)/2)%>%
  select(site_code, project_name, community_type, species_matched, mean, freq, DCi)%>%
  rename(cmean=mean, cfreq=freq)

#getting frequency of treated plots, same code as above but for treated plots
treatplots<-allreldat%>%
  filter(plot_mani!=0)%>%
  select(site_code, project_name, community_type, treatment, plot_id)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarize(ntplots=length(plot_id))

treat_freq<-allreldat%>%
  filter(plot_mani!=0)%>%
  select(site_code, project_name, community_type, treatment, species_matched, plot_id)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment, species_matched)%>%
  summarize(nplots=length(plot_id))%>%
  left_join(treatplots)%>%
  mutate(freq=nplots/ntplots)

treat_dom<-treat_freq%>%
  left_join(Trelave)%>%
  mutate(treatDCi=(mean+freq)/2)%>%
  select(site_code, project_name, community_type, treatment, species_matched, mean, freq, treatDCi)

#Exploring DCi

#get contol and treated, realtive abundacne (mean), relative frequency (freq), and DCi into one dataset 
abund_freqc<-control_freq%>%
  left_join(Crelave)%>%
  ungroup()%>%
  select(species_matched, mean, freq)%>%
  mutate(dci=(mean+freq)/2)
abund_freqt<-treat_freq%>%
  left_join(Trelave)%>%
  ungroup()%>%
  select(species_matched, mean, freq)%>%
  mutate(dci=(mean+freq)/2)
abund_freq<-abund_freqc%>%
  bind_rows(abund_freqt) 

#are abund and freq related? weakly correlated
with(abund_freq, plot(mean, freq))
with(abund_freq, cor.test(mean, freq))

#how is rel freq related to DCI? strongly correlated, r= 0.95
with(abund_freq, plot(freq, dci))
with(abund_freq, cor.test(dci, freq))
#how is rel abund related to DCI? strongly correlated, but not as strong as frequency r = 0.55
with(abund_freq, plot(mean, dci))
with(abund_freq, cor.test(dci, mean))

with(abund_freq, hist(dci))

#how to determine rare/dominat species
with(abund_freq, quantile(dci, probs=c(0.1, 0.9)))

approach1<-abund_freq%>%
  mutate(sptype=ifelse(dci<0.102, "rare", ifelse(dci>0.554, "dom", "present")))%>%
  group_by(sptype)%>%
  summarise(n=length(dci),
            pct=n/38458)

approach2<-abund_freq%>%
  mutate(sptype=ifelse(dci<0.2, "rare", ifelse(dci>0.56, "dom", "present")))%>%
  group_by(sptype)%>%
  summarise(n=length(dci),
            pct=n/38458)

with(abund_freq, hist(mean))
with(abund_freq, hist(freq))

approach3<-abund_freq%>%
  mutate(sptype=ifelse(dci<0.08, "rare", ifelse(dci>0.63, "dom", "present")))%>%
  group_by(sptype)%>%
  summarise(n=length(dci),
            pct=n/38458)

with(abund_freq, hist(dci))
# 
# toplot_test1<-Trelave%>%
#   full_join(control_dom)
# 
# toplot$cmean[is.na(test2$cmean)] <- 0
# toplot$mean[is.na(test2$mean)] <- 0
# 
# toplot2<-toplot%>%
#   mutate(logrr=log10((mean+1)/(cmean+1)))%>%
#   mutate(spcolor=ifelse(species_matched=="Andropogon gerardii", "ag", ifelse(species_matched=="Ambrosia psilostachya", "ap", ifelse(species_matched=="Triodanis perfoliata", "tp", "a"))))
# 
# abundplot_dom<-
#   
#   ggplot(data=toplot2, aes(x=DCi, y=log10(mean+1), color=spcolor))+
#   geom_point(data=subset(toplot2, spcolor=="a"))+
#   geom_point(data=subset(toplot2, spcolor!="a"))+
#   scale_color_manual(values=c("gray", "red", 'blue', "darkgreen"))+
#   geom_abline(aes(intercept=0, slope=1))
# 
# log_dom<-ggplot(data=toplot2, aes(x=DCi, y=logrr, color=spcolor))+
#   geom_point(data=subset(toplot2, spcolor=="a"))+
#   geom_point(data=subset(toplot2, spcolor!="a"))+
#   scale_color_manual(values=c("gray", "red", 'blue', "darkgreen"))+
#   geom_hline(yintercept=0)
# 
# grid.arrange(abundplot_dom, log_dom, ncol=2)

# Calculating changes in DCi for each species comparing control to treatments.
CT<-treat_dom%>%
  full_join(control_dom)#this introduces NA into the treatment column when it is present in the controls but not the treated plots.

#this is to replace the NAs with 0
CT$DCi[is.na(CT$DCi)] <- 0
CT$treatDCi[is.na(CT$treatDCi)] <- 0

# CT_plot<-CT%>%
#   mutate(spcolor=ifelse(species_matched=="Andropogon gerardii", "ag", ifelse(species_matched=="Ambrosia psilostachya", "ap", ifelse(species_matched=="Linum sulcatum", "ls", "a"))))
# 
# 
# abundplot_dom2<-
#   
#   ggplot(data=CT_plot, aes(x=DCi, y=treatDCi, color=spcolor))+
#   geom_point(data=subset(CT_plot, spcolor=="a"))+
#   geom_point(data=subset(CT_plot, spcolor!="a"))+
#   scale_color_manual(values=c("gray", "red", 'blue', "darkgreen"))+
#   geom_abline(aes(intercept=0, slope=1))
# 
# grid.arrange(abundplot, abundplot_dom2, ncol=2)

# #getting residuals
# residuals<-as.data.frame(with(toplot2_test2,(resid(lm(treatDCi-DCi~0)))))
# 
# #Getting residuals from a 1:1 line is the same as subtracting x-y.
# #postive values, increased abundance in treatment
# #negative values, decreased abudnance in treatment

CT_diff<-CT%>%
  mutate(diff=treatDCi-DCi,
         absdiff=abs(diff))
with(CT_diff, hist(diff))
with(CT_diff, hist(absdiff))

CT_pd<-CT%>%
  filter(treatDCi!=0|DCi!=0)%>%
  mutate(pd=(treatDCi-DCi)/DCi)
with(CT_pd, hist(pd))


#exploring a cut-off for winners/losers/tolerators
#I am basing these on absolute difference..

#how does abundance, and freqency in control pltos affect difference?
with(CT_diff, plot(cmean, absdiff))
with(CT_diff, cor.test(cmean, absdiff))#no relationship

with(CT_diff, plot(cfreq, absdiff))
with(CT_diff, cor.test(cfreq, absdiff))# no relationship

#how does abundance, and freqency in treated pltos affect difference?
with(CT_diff, plot(mean, absdiff))
with(CT_diff, cor.test(mean, absdiff))#no relationship

with(CT_diff, plot(freq, absdiff))
with(CT_diff, cor.test(freq, absdiff))# no relationship

###how does DCi relate to change in rel abund or rel freq?

CT_diff2<-CT_diff%>%
  mutate(dfreq=abs(freq-cfreq),
         dabund=abs(mean-cmean))

with(CT_diff2, plot(dfreq, absdiff))
abline(h=0.1, col="red")

with(CT_diff2, plot(dabund, absdiff))
abline(h=0.1, col="red")


##tolerator cut-off
options<-CT_diff%>%
  mutate(response1=ifelse(absdiff<0.1, "t", "wl"),
         response5=ifelse(absdiff<0.05, "t", "wl"),
         response2=ifelse(absdiff<0.2, "t", "wl"))
opt1<-options%>%
  group_by(response1)%>%
  summarize(n=length(absdiff),
            p=n/34240)

opt2<-options%>%
  group_by(response5)%>%
  summarize(n=length(absdiff),
            p=n/34240)

opt3<-options%>%
  group_by(response2)%>%
  summarize(n=length(absdiff),
            p=n/34240)


#For this analysis I am using only species that are found in treatments and controls, using differneces btwn trt controls, and the more stringent approach to rare/dominant, 
categories<-CT_diff%>%
  na.omit%>%
  filter(treatDCi!=0|DCi!=0)%>%
  mutate(cat=ifelse(DCi>0.63&treatDCi<0.08, "Sup. Lose",
                    ifelse(DCi<0.08&treatDCi<0.08, "Rare Tol.", 
                           ifelse(DCi>0.63&treatDCi>0.63, "Dom. Tol.", 
                                  ifelse(DCi>0.63&treatDCi<0.63, "Loose Dom.",
                                         ifelse(DCi<0.08&treatDCi>0.08, "No Long. Rare",
                                                ifelse(diff>0.1, "Winner", 
                                                       ifelse(absdiff>0.1&diff<0, "Loser", 
                                                              ifelse(absdiff<0.1, "Tol.", "none")))))))),
         DR=ifelse(DCi>0.63, "Dom", 
                   ifelse(DCi<0.08, "Rare", "none")),
         Pres=ifelse(DCi>0.08&DCi<0.63, 1, 0))

#old code
# ifelse(DCi!=0&DCi<0.08&treatDCi>0.6, "Sup. Win.", 
#        ifelse(DCi==0&treatDCi>0.6, "Sup. Dup. Win", 
# ifelse(DCi==0&treatDCi<0.08, "Appear Rare",
#        ifelse(DCi==0&treatDCi>0.08, "Appear", 
# AppDis=ifelse(DCi==0&treatDCi>0, "Appear",
#               ifelse(DCi>0&treatDCi==0, "Disapp.", "none")),

#getting a table of responses
cat_sum<-categories%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/26497)

# dr_sum<-categories%>%
#   group_by(DR)%>%
#   summarize(n=length(DR))
# 
# ad_sum<-categories%>%
#   group_by(AppDis)%>%
#   summarize(n=length(AppDis))
# 
# sum(categories$Pres)


###contrasting different treatments
trts2<-trts%>%
  select(site_code, project_name, community_type, treatment, trt_type)%>%
  unique

cat_trts<-categories%>% 
  left_join(trts2)%>%
  mutate(cat2=ifelse(cat=="Appear"|cat=="Winner"|cat=="No Long. Rare"|cat=="Appear Rare"|cat=="Sup. Dup. Win", "Winner", ifelse(cat=="Loser"|cat=="Loose Dom.", "Loser", ifelse(cat=="Tol."|cat=="Dom. Tol."|cat=="Rare Tol.", "Tol.", "None"))))%>%
  group_by(trt_type)%>%
  summarize(avediff=mean(diff))

ntrt<-cat_trts%>%
  filter(trt_type=="N")
ntrt_sum<-ntrt%>%
  group_by(trt_type, cat2)%>%
  summarize(n=length(cat2))%>%
  mutate(prop=n/2249)

drttrt<-cat_trts%>%
  filter(trt_type=="drought")
drttrt_sum<-drttrt%>%
  group_by(trt_type,cat2)%>%
  summarize(n=length(cat2))%>%
  mutate(prop=n/939)

irrtrt<-cat_trts%>%
  filter(trt_type=="irr")
irrtrt_sum<-irrtrt%>%
  group_by(trt_type,cat2)%>%
  summarize(n=length(cat2))%>%
  mutate(prop=n/798)

multnuttrt<-cat_trts%>%
  filter(trt_type=="mult_nutrient")
multnuttrt_sum<-multnuttrt%>%
  group_by(trt_type,cat2)%>%
  summarize(n=length(cat2))%>%
  mutate(prop=n/3548)

trt_contrast<-ntrt_sum%>%
  bind_rows(drttrt_sum, irrtrt_sum, multnuttrt_sum)

theme_set(theme_bw(16))

ggplot(data=trt_contrast, aes(x=trt_type, y=prop, fill=cat2))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Treatment")+
  ylab("Proportion of Sp. Responses")+
  scale_fill_manual(name="Response", values=c("red", "gray",'blue'))

ggplot(data=filter(cat_trts, trt_type=="N"|trt_type=="mult_nutrient"|trt_type=="irr"|trt_type=="drought"), aes(x=trt_type, y=avediff))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Treatment")+
  ylab("Average Control-Trt Diff")

#exploring species specfic responses
###example with A. gerardii
ange<-categories%>%
  filter(species_matched=="Andropogon gerardii")%>%
  left_join(trts2)

angecat_sum<-ange%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/130)

ange_trt<-ange%>%
  filter(trt_type=="N"|trt_type=="mult_nutrient"|trt_type=="irr"|trt_type=="drought")%>%
  group_by(trt_type)%>%
  summarize(ave_diff=mean(diff),
            n=length(diff),
            sd=sd(diff))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=ange_trt, aes(x=trt_type, y=ave_diff, label=n))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Treatment")+
  ylab("Average Control-Trt Diff")+
  geom_label(hjust=1)+
  geom_errorbar(aes(ymin=ave_diff-se, ymax=ave_diff+se), position = position_dodge(0.9), width=0.5)



#common species - get list of common species to explore
common<-categories%>%
  group_by(species_matched)%>%
  summarize(n=length(species_matched))

##erigeron canadensis
erca<-categories%>%
  filter(species_matched=="Erigeron canadensis")%>%
  left_join(trts2)

ercacat_sum<-erca%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/209)

erca_trt<-erca%>%
  filter(trt_type=="N"|trt_type=="mult_nutrient"|trt_type=="irr"|trt_type=="drought")%>%
  group_by(trt_type)%>%
  summarize(ave_diff=mean(diff),
            n=length(diff),
            sd=sd(diff))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=erca_trt, aes(x=trt_type, y=ave_diff, label=n))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Treatment")+
  ylab("Average Control-Trt Diff")+
  geom_label(hjust=1)+
  geom_errorbar(aes(ymin=ave_diff-se, ymax=ave_diff+se), position = position_dodge(0.9), width=0.5)


#for poa pratensis
popr<-categories%>%
  filter(species_matched=="Poa pratensis")%>%
  left_join(trts2)

poprcat_sum<-popr%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/233)

popr_trt<-popr%>%
  filter(trt_type=="N"|trt_type=="mult_nutrient"|trt_type=="irr"|trt_type=="drought")%>%
  group_by(trt_type)%>%
  summarize(ave_diff=mean(diff),
            n=length(diff),
            sd=sd(diff))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=popr_trt, aes(x=trt_type, y=ave_diff, label=n))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Treatment")+
  ylab("Average Control-Trt Diff")+
  geom_label(hjust=1)+
  geom_errorbar(aes(ymin=ave_diff-se, ymax=ave_diff+se), position = position_dodge(0.9), width=0.5)


