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
trts<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE_treatment_summary.csv")%>%
  select(-X)
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
  left_join(sp)%>%# this drops the unknowns
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
  full_join(control_dom)#this introduces NA into the treatment column when it is present in the controls but not the treated plots.

CT$DCi[is.na(CT$DCi)] <- 0
CT$treatDCi[is.na(CT$treatDCi)] <- 0

CT_plot<-CT%>%
  mutate(spcolor=ifelse(species_matched=="Andropogon gerardii", "ag", ifelse(species_matched=="Ambrosia psilostachya", "ap", ifelse(species_matched=="Linum sulcatum", "ls", "a"))))


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
         absdiff=abs(diff))

categories<-CT_diff%>%
  na.omit%>%
  mutate(cat=ifelse(DCi!=0&DCi<0.08&treatDCi>0.6, "Sup. Win.", 
             ifelse(DCi==0&treatDCi>0.6, "Sup. Dup. Win", 
             ifelse(DCi>0.6&treatDCi!=0&treatDCi<0.1, "Sup. Lose",
             ifelse(DCi!=0&DCi<0.08&treatDCi!=0&treatDCi<0.08, "Rare Tol.", 
             ifelse(DCi>0.6&treatDCi>0.6, "Dom. Tol.", 
             ifelse(DCi==0&treatDCi<0.08, "Appear Rare",
             ifelse(DCi==0&treatDCi>0.08, "Appear", 
             ifelse(DCi>0.6&treatDCi<0.6, "Loose Dom.",
             ifelse(DCi<0.08&treatDCi>0.08, "No Long. Rare",
             ifelse(diff>0.1, "Winner", 
             ifelse(absdiff>0.1&diff<0, "Loser", 
             ifelse(absdiff<0.1, "Tol.", "none")))))))))))),
         DR=ifelse(DCi>0.6, "Dom", 
            ifelse(DCi!=0&DCi<0.08, "Rare", "none")),
         AppDis=ifelse(DCi==0&treatDCi>0, "Appear",
                ifelse(DCi>0&treatDCi==0, "Disapp.", "none")),
         Pres=ifelse(DCi!=0&DCi>0.08&DCi<0.6, 1, 0))

cat_sum<-categories%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/17670)

dr_sum<-categories%>%
  group_by(DR)%>%
  summarize(n=length(DR))

ad_sum<-categories%>%
  group_by(AppDis)%>%
  summarize(n=length(AppDis))

sum(categories$Pres)

###example with A. gerardii
ange<-categories%>%
  filter(species_matched=="Andropogon gerardii")

angecat_sum<-ange%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/103)

#common species
common<-categories%>%
  group_by(species_matched)%>%
  summarize(n=length(species_matched))


erca<-categories%>%
  filter(species_matched=="Erigeron canadensis")

ercacat_sum<-erca%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/152)

popr<-categories%>%
  filter(species_matched=="Poa pratensis")
poprcat_sum<-popr%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/157)

###contrasting n-addition versus drought
trts2<-trts%>%
  select(site_code, project_name, community_type, treatment, trt_type)%>%
  unique


cat_trts<-categories%>% #this is a problem for when species are lost b/c control has no treatment.
  left_join(trts2)%>%
  mutate(cat2=ifelse(cat=="Appear"|cat=="Winner"|cat=="No Long. Rare"|cat=="Appear Rare"|cat=="Sup. Dup. Win", "Winner", ifelse(cat=="Loser"|cat=="Loose Dom.", "Loser", ifelse(cat=="Tol."|cat=="Dom. Tol."|cat=="Rare Tol.", "Tol.", "None"))))

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

ggplot(data=trt_contrast, aes(x=cat2, y=prop, fill=trt_type))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Species Response")+
  ylab("Proportion of Sp. Responses")+
  scale_fill_manual(name="Treatment", values=c("lightblue", "darkblue",'darkred','red'))

####now doing this for every year of an experiment.

relave_yr<-allreldat%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani, treatment_year, species_matched)%>%
  summarize(mean=mean(relcov))

Crelave_yr<-relave_yr%>%
  filter(plot_mani==0)%>%
  ungroup()%>%
  select(-treatment, -plot_mani)%>%
  rename(cmean=mean)

Trelave_yr<-relave_yr%>%
  filter(plot_mani!=0)%>%
  ungroup()%>%
  select(-plot_mani)

controlplots_yr<-allreldat%>%
  filter(plot_mani==0)%>%
  select(site_code, project_name, community_type, plot_id, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment_year)%>%
  summarize(ncplots=length(plot_id))

control_freq_yr<-allreldat%>%
  filter(plot_mani==0)%>%
  select(site_code, project_name, community_type, species_matched, plot_id, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, species_matched, treatment_year)%>%
  summarize(nplots=length(plot_id))%>%
  left_join(controlplots)%>%
  mutate(freq=nplots/ncplots)

control_dom_yr<-control_freq_yr%>%
  left_join(Crelave_yr)%>%
  mutate(DCi=(cmean+freq)/2)%>%
  select(site_code, project_name, community_type, species_matched, treatment_year, DCi)

treatplots_yr<-allreldat%>%
  filter(plot_mani!=0)%>%
  select(site_code, project_name, community_type, treatment, plot_id, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment, treatment_year)%>%
  summarize(ntplots=length(plot_id))

treat_freq_yr<-allreldat%>%
  filter(plot_mani!=0)%>%
  select(site_code, project_name, community_type, treatment, species_matched, plot_id, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment, species_matched, treatment_year)%>%
  summarize(nplots=length(plot_id))%>%
  left_join(treatplots)%>%
  mutate(freq=nplots/ntplots)

treat_dom_yr<-treat_freq_yr%>%
  left_join(Trelave_yr)%>%
  mutate(treatDCi=(mean+freq)/2)%>%
  select(site_code, project_name, community_type, treatment, species_matched, treatment_year, treatDCi)

CT_yr<-treat_dom_yr%>%
  full_join(control_dom_yr)#this introduces NA into the treatment column when it is present in the controls but not the treated plots.

CT_yr$DCi[is.na(CT_yr$DCi)] <- 0
CT_yr$treatDCi[is.na(CT_yr$treatDCi)] <- 0


CT_diff_yr<-CT_yr%>%
  mutate(diff=treatDCi-DCi,
         absdiff=abs(diff))

categories_yr<-CT_diff_yr%>%
  na.omit%>%
  mutate(cat=ifelse(DCi!=0&DCi<0.08&treatDCi>0.6, "Sup. Win.", 
             ifelse(DCi==0&treatDCi>0.6, "Sup. Dup. Win", 
             ifelse(DCi>0.6&treatDCi!=0&treatDCi<0.1, "Sup. Lose",
             ifelse(DCi!=0&DCi<0.08&treatDCi!=0&treatDCi<0.08, "Rare Tol.", 
             ifelse(DCi>0.6&treatDCi>0.6, "Dom. Tol.", 
             ifelse(DCi==0&treatDCi<0.08, "Appear Rare",
             ifelse(DCi==0&treatDCi>0.08, "Appear", 
             ifelse(DCi>0.6&treatDCi<0.6, "Loose Dom.",
             ifelse(DCi<0.08&treatDCi>0.08, "No Long. Rare",
             ifelse(diff>0.1, "Winner", 
             ifelse(absdiff>0.1&diff<0, "Loser", 
             ifelse(absdiff<0.1, "Tol.", "none")))))))))))),
         DR=ifelse(DCi>0.6, "Dom", 
                   ifelse(DCi!=0&DCi<0.08, "Rare", "none")),
         AppDis=ifelse(DCi==0&treatDCi>0, "Appear",
                       ifelse(DCi>0&treatDCi==0, "Disapp.", "none")),
         Pres=ifelse(DCi!=0&DCi>0.08&DCi<0.6, 1, 0),
         cat2=ifelse(cat=="Appear"|cat=="Winner"|cat=="No Long. Rare"|cat=="Appear Rare"|cat=="Sup. Dup. Win", "Winner", ifelse(cat=="Loser"|cat=="Loose Dom."|cat=="Sup. Lose", "Loser", ifelse(cat=="Tol."|cat=="Dom. Tol."|cat=="Rare Tol.", "Tol.", "None"))))

cat_sum_yr<-categories_yr%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/72221)

dr_sum_yr<-categories_yr%>%
  group_by(DR)%>%
  summarize(n=length(DR))

ad_sum_yr<-categories_yr%>%
  group_by(AppDis)%>%
  summarize(n=length(AppDis))

cat2_sum_yr<-categories_yr%>%
  group_by(cat2)%>%
  summarise(n=length(cat2))%>%
  mutate(prop=n/72221)
