#investigate winners losers

library(tidyverse)
library(gridExtra)


my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"

dat<-read.csv(paste(my.wd, "CoRRE data/CoRRE_raw_abundance_Nov2019.csv",sep=""))

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
  summarize(mean=mean(relcov))

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



###for a. gerardii
trts<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE_treatment_summary.csv")
sp<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE_TRY_species_list.csv")%>%
  select(genus_species, species_matched)%>%
  unique

sp$genus_species <- trimws(sp$genus_species, which="right")

ag_test<-dat%>%
  select(-X)%>%
  #filter(genus_species=="andropogon gerardii")%>%
  left_join(trts)%>%
  left_join(sp)%>%
  na.omit()
  
  
ave<-ag_test%>%
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
  mutate(spcolor=ifelse(species_matched=="Andropogon gerardii", "ag", ifelse(species_matched=="Ambrosia psilostachya", "ap", ifelse(species_matched=="Triodanis perfoliata", "tp", "z"))))

#how to deal with species that are not in BOTH treatment and control? 
#option 1, make thier abundance very low 
#option 2, add 1 to all values
#option 3 ???
###want the x-axis to be one of dominace not just abundnace in controls and  to include frequency across the landscape
abundplot<-
  
  ggplot(data=test3, aes(x=log10(cmean+1), y=log10(mean+1), color=spcolor))+
  geom_point()+
  scale_color_manual(values=c("red", "blue", "green4", "gray"))+
  geom_abline(aes(intercept=0, slope=1))

log<-ggplot(data=test3, aes(x=log10(cmean+1), y=logrr, color=spcolor))+
  geom_point()+
  scale_color_manual(values=c("red", "blue", "green4", "gray"))+
  geom_hline(yintercept=0)

grid.arrange(abundplot, log, ncol=2)
