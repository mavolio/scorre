#investigate winners losers

library(tidyverse)
library(gridExtra)


my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"

dat<-read.csv(paste(my.wd, "CoRRE data/CoRRE_raw_abundance_Nov2019.csv",sep=""))

reldat<-read.csv(paste(my.wd, "CoRRE data/CoRRE_relative_abundance_Nov2019.csv",sep=""))

pplots<-dat%>%
  filter(site_code=="KNZ"&project_name=="pplots")%>%
  select(-X)%>%
  filter(treatment=="N1P0"|treatment=="N2P3")

ave<-pplots%>%
  group_by(treatment_year, treatment, genus_species)%>%
  summarise(mcov=mean(abundance))

c<-ave%>%
  filter(treatment=="N1P0")%>%
  rename(ccov=mcov)%>%
  ungroup()%>%
  select(-treatment)

ctdiff<-ave%>%
  filter(treatment!="N1P0")%>%
  left_join(c)

ctdiff$ccov[is.na(ctdiff$ccov)] <- 0

ctdiff<-ctdiff%>%
  mutate(diff=ccov-mcov)

# splist<-unique(ctdiff$genus_species)
# 
# slopes<-data.frame()
# 
# for i in 1:length(splist){
#   
#   subset<-ctdiff%>%
#     filter(genus_species==splist[i])
#   
#   
#}

test<-ctdiff%>%
  filter(genus_species=="ambrosia psilostachya")
  
ggplot(data=test, aes(x=treatment_year, y=diff))+
  geom_point()

pplotsw<-dat%>%
  filter(site_code=="KNZ"&project_name=="pplots")%>%
  select(-X)%>%
  filter(treatment=="N1P0"|treatment=="N2P3")%>%
  filter(genus_species=="ambrosia psilostachya")%>%
  mutate(time=ifelse(treatment_year<6, 1, 2))

ggplot(data=pplotsw, aes(x=treatment, y=abundance))+
  geom_point()+
  facet_wrap(~time)


pplotsw<-dat%>%
  filter(site_code=="KNZ"&project_name=="pplots")%>%
  select(-X)%>%
  filter(treatment=="N1P0"|treatment=="N2P3")%>%
  filter(genus_species=="ambrosia psilostachya")%>%
  mutate(time=ifelse(treatment_year<6, 1, 2))

ave_ct<-pplotsw%>%
  group_by(treatment, time)%>%
  summarize(mean=mean(abundance))

c1<-ave_ct[1, 3]
t1<-ave_ct[3, 3]

logrr<-log(t1/c1)


pplots_test<-dat%>%
  filter(site_code=="KNZ"&project_name=="pplots")%>%
  select(-X)%>%
  filter(treatment=="N1P0"|treatment=="N2P3")%>%
  mutate(time=ifelse(treatment_year<6, 1, 2))%>%
  group_by(treatment, time, genus_species)%>%
  summarize(mean=mean(abundance))

testC<-pplots_test%>%
  filter(treatment=="N1P0")%>%
  ungroup()%>%
  select(-treatment)%>%
  rename(cmean=mean)


test2<-pplots_test%>%
  filter(treatment!="N1P0")%>%
  ungroup()%>%
  select(-treatment)%>%
  full_join(testC)

test2$cmean[is.na(test2$cmean)] <- 0
test2$mean[is.na(test2$mean)] <- 0

test3<-test2%>%
  mutate(logrr=log10((mean+1)/(cmean+1)))%>%
  filter(time==2)

#how to deal with species that are not in BOTH treatment and control? 
#option 1, make thier abundance very low 
#option 2, add 1 to all values
#option 3 ???

###want the x-axis to be one of dominace not just abundnace in controls and  to include frequency across the landscape
abundplot<-ggplot(data=test3, aes(x=log10(cmean+1), y=log10(mean+1)))+
  geom_point()+
  geom_abline(aes(intercept=0, slope=1))

log<-ggplot(data=test3, aes(x=log10(cmean+1), y=logrr))+
  geom_point()+
  geom_hline(yintercept=0)

grid.arrange(abundplot, log, ncol=2)



###for all the data
trts<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE_treatment_summary.csv")
sp<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE_TRY_species_list.csv")%>%
  select(genus_species, species_matched)%>%
  unique

sp$genus_species <- trimws(sp$genus_species, which="right")

alldat<-dat%>%
  select(-X)%>%
  #filter(genus_species=="andropogon gerardii")%>%
  left_join(trts)%>%
  left_join(sp)%>%
  na.omit()
  
  
ave<-alldat%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani, species_matched)%>%
  summarize(mean=mean(abundance))

testC<-ave%>%
  filter(plot_mani==0)%>%
  ungroup()%>%
  select(-treatment, -plot_mani)%>%
  rename(cmean=mean)
test2<-ave%>%
  filter(plot_mani!=0)%>%
  ungroup()%>%
  select(-treatment, -plot_mani)%>%
  full_join(testC)

test2$cmean[is.na(test2$cmean)] <- 0
test2$mean[is.na(test2$mean)] <- 0

test3<-test2%>%
  mutate(logrr=log10((mean+1)/(cmean+1)))%>%
  mutate(spcolor=ifelse(species_matched=="Andropogon gerardii", "ag", ifelse(species_matched=="Ambrosia psilostachya", "ap", ifelse(species_matched=="Triodanis perfoliata", "tp", "a"))))

#how to deal with species that are not in BOTH treatment and control? 
#option 1, make thier abundance very low 
#-->option 2, add 1 to all values - we did this b/c we add 1 to everything.
#option 3 ???
###want the x-axis to be one of dominace not just abundnace in controls and  to include frequency across the landscape
abundplot<-
  
  ggplot(data=test3, aes(x=log10(cmean+1), y=log10(mean+1), color=spcolor))+
  geom_point(data=subset(test3, spcolor=="a"))+
  geom_point(data=subset(test3, spcolor!="a"))+
  scale_color_manual(values=c("gray", "red", 'blue', "darkgreen"))+
  geom_abline(aes(intercept=0, slope=1))

log<-ggplot(data=test3, aes(x=log10(cmean+1), y=logrr, color=spcolor))+
  geom_point(data=subset(test3, spcolor=="a"))+
  geom_point(data=subset(test3, spcolor!="a"))+
  scale_color_manual(values=c("gray", "red", 'blue', "darkgreen"))+
  geom_hline(yintercept=0)

grid.arrange(abundplot, log, ncol=2)


##now doing this with dominance on the x-axis and y-axis
#facet wrap by treatment type.
allreldat<-reldat%>%
  select(-X)%>%
  left_join(trts)%>%
  left_join(sp)%>%
  na.omit()

relave<-allreldat%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani, species_matched)%>%
  summarize(mean=mean(relcov))

Crelave<-relave%>%
  filter(plot_mani==0)%>%
  ungroup()%>%
  select(-treatment, -plot_mani)%>%
  rename(cmean=mean)

Trelave<-relave%>%
  filter(plot_mani!=0)%>%
  ungroup()%>%
  select(-plot_mani)

controlplots<-allreldat%>%
  filter(plot_mani==0)%>%
  select(site_code, project_name, community_type, plot_id)%>%
  unique()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(ncplots=length(plot_id))
  
control_freq<-allreldat%>%
  filter(plot_mani==0)%>%
  select(site_code, project_name, community_type, species_matched, plot_id)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, species_matched)%>%
  summarize(nplots=length(plot_id))%>%
  left_join(controlplots)%>%
  mutate(freq=nplots/ncplots)

control_dom<-control_freq%>%
  left_join(Crelave)%>%
  mutate(DCi=(cmean+freq)/2)%>%
  select(site_code, project_name, community_type, species_matched, DCi)

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
  select(site_code, project_name, community_type, treatment, species_matched, treatDCi)

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

CT<-treat_dom%>%
  full_join(control_dom)

CT$DCi[is.na(CT$DCi)] <- 0
CT$treatDCi[is.na(CT$treatDCi)] <- 0

CT_plot<-CT%>%
  mutate(spcolor=ifelse(species_matched=="Andropogon gerardii", "ag", ifelse(species_matched=="Ambrosia psilostachya", "ap", ifelse(species_matched=="Triodanis perfoliata", "tp", "a"))))


abundplot_dom2<-
  
  ggplot(data=CT_plot, aes(x=DCi, y=treatDCi, color=spcolor))+
  geom_point(data=subset(CT_plot, spcolor=="a"))+
  geom_point(data=subset(CT_plot, spcolor!="a"))+
  scale_color_manual(values=c("gray", "red", 'blue', "darkgreen"))+
  geom_abline(aes(intercept=0, slope=1))

grid.arrange(abundplot, abundplot_dom2, ncol=2)

#getting residuals
residuals<-as.data.frame(with(toplot2_test2,(resid(lm(treatDCi-DCi~0)))))

#Getting residuals from a 1:1 line is the same as subtracting x-y.
#postive values, increased abundance in treatment
#negative values, decreased abudnance in treatment

CT_diff<-CT%>%
  mutate(diff=treatDCi-DCi,
         absdiff=abs(diff))%>%
  select(-nplots, -ncplots, -freq, -cmean)

categories<-CT_diff%>%
  mutate(cat=ifelse(DCi!=0&DCi<0.08&treatDCi>0.6, "Sup. Win.", 
             ifelse(DCi==0&treatDCi>0.6, "Sup. Dup. Win", 
             ifelse(DCi>0.6&treatDCi==0, "Sup. Dup. Lose",
             ifelse(DCi>0.6&treatDCi!=0&treatDCi<0.1, "Sup. Lose",
             ifelse(DCi!=0&DCi<0.08&treatDCi!=0&treatDCi<0.08, "Rare Tol.", 
             ifelse(DCi>0.6&treatDCi>0.6, "Dom. Tol.", 
             ifelse(DCi==0&treatDCi<0.8, "Appear Rare",
             ifelse(DCi==0&treatDCi>0.8, "Appear", 
             ifelse(DCi!=0&DCi<0.08&treatDCi==0,"Disapp. Rare",
             ifelse(DCi>0.08&treatDCi==0, "Disapp.", 
             ifelse(diff>0.1, "Winner", 
             ifelse(absdiff>0.1&diff<0, "Loser", 
             ifelse(absdiff<0.1, "Tol.", "none"))))))))))))),
         DR=ifelse(DCi>0.6, "Dom", 
            ifelse(DCi!=0&DCi<0.08, "Rare", "none")))

dr_sum<-categories%>%
  group_by(DR)%>%
  summarize(n=length(DR))

cat_sum<-categories%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/17984)
         