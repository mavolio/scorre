library(tidyverse)
library(codyn)

setwd('C:\\Users\\mavolio2\\Dropbox\\CoRRE_database\\Data\\CompiledData')

theme_set(theme_bw(12))

dat<-read.csv('RawAbundance.csv')
trts<-read.csv('ExperimentInfo.csv')

alldat<-dat %>% 
  left_join(trts)

proddat<-read.csv('ANPP2021.csv')

#replicate communities - maybe switch foci a bit on diff community types
N_commtypes<-alldat%>%
  filter(community_type!=0&n>0, site_code!='ASGA'&site_code!='KUFS') %>% 
  group_by(site_code, project_name, community_type, treatment) %>% 
  summarize(max=max(treatment_year), n=length(treatment_year)) %>% 
  filter(max>5) %>% 
  select(site_code, project_name, community_type) %>%
  unique()

SubsetDat1<-alldat %>% 
  filter(site_code=='KNZ'&project_name=='change'|site_code=='KNZ'&project_name=='pplots'|site_code=='KNZ'&project_name=='BGP'|site_code=='KNZ'&project_name=='RPHs'|site_code=='ARC'|site_code=='DL'&project_name=='NSFC'|site_code=='DL'&project_name=='GCME'|site_code=='DL'&project_name=='GCME2'|site_code=='NWT'&project_name=='snow')

SubsetDat2<-alldat %>% 
  right_join(N_commtypes)

SubsetDat_trts<-SubsetDat1 %>% 
  bind_rows(SubsetDat2) %>% 
  filter(plot_mani==0|n>0) %>% 
  select(site_code, project_name, community_type, treatment, plot_mani, n, p, k) %>% 
  unique()

SubsetDat<-SubsetDat1 %>% 
  bind_rows(SubsetDat2) %>% 
  mutate(keep=
           ifelse(site_code=='KNZ'&project_name=='BGP'&treatment %in% c('u_u_c','u_u_n'), 1,
           ifelse(site_code=='KNZ'&project_name=='pplots'&treatment %in% c('N1P0','N2P0'), 1, 
           ifelse(site_code=='KNZ'&project_name=='RPHs'&treatment %in% c('control','N'), 1, 
           ifelse(site_code=='KNZ'&project_name=='change'&treatment %in% c('0','5', '10'), 1,
           ifelse(site_code %in% c('ARC', 'CAR'),  1, 
           ifelse(site_code=='CDR'&project_name=='e001', 1, 
           ifelse(site_code=='DL'&project_name=='GCME'&treatment %in% c('C', 'N', 'NP'), 1,
           ifelse(site_code=='DL'&project_name=='GCME2'&treatment %in% c('C', 'N', 'NP'), 1,
           ifelse(site_code=='DL'&project_name=='NSFC'&treatment %in% c('C', 'N'), 1, 
           ifelse(site_code=='LATNJA'&treatment %in% c('CONTROL', 'N'), 1,
           ifelse(site_code=='NWT'&project_name=='snow'&treatment %in% c('XXX', 'XNX'), 1,
           ifelse(site_code=='NWT'&project_name=='bowman', 1, 
           ifelse(site_code=='CAU'&treatment %in% c('Cont', 'N', 'NP'), 1, 0)))))))))))))) %>% 
  filter(keep==1)


check<-SubsetDat %>% 
  select(site_code, project_name, community_type, treatment, plot_mani, n, p, k) %>% 
  unique()

###make control always be control
newdat<-SubsetDat %>% 
  mutate(spc=paste(site_code, project_name, community_type, sep='::')) %>% 
  select(site_code, project_name, community_type, spc, calendar_year, treatment_year, treatment, plot_id, block, abundance, genus_species, plot_mani, n, p, k, CO2, mow_clip, successional, trt_type) %>% 
  mutate(treat2=ifelse(plot_mani==0, 'control', ifelse(site_code=='NWT'&n>0&p==0, 'n_25', ifelse(treatment %in% c('N', 'N2P0', 'u_u_n', 'XNX', 1:10), paste('n', n, sep="_"), treatment)))) %>% 
  mutate(plot_id2=ifelse(project_name=='NSFC', paste(treatment, plot_id, sep='_'), plot_id))

info<-newdat %>% 
  select(site_code, project_name, community_type, spc, treatment_year, calendar_year, treatment, treat2, plot_id, block, plot_mani, n, p, k, trt_type) %>% 
  unique()

list<-unique(newdat$spc)

table(newdat$treat2)

###seeing number of replicated treatments at each site
numexp<-newdat %>% 
  select(site_code, spc, treat2) %>% 
  unique() %>% 
  group_by(site_code, treat2) %>% 
  summarize(n=length(spc))

#there are 28 C-T comparisons across   


###Community outcomes

#how many have pre-treatment data? only 3, so not enough to use pretreatment data
pretrt<-newdat %>% 
  filter(treatment_year==0) %>% 
  select(site_code, spc) %>% 
  unique()

#doing control-treatment differences over time taking into account blocking when possible

blocked<-info %>% 
  select(spc, block) %>% 
  unique()

test<-newdat %>% 
  filter(site_code=='KNZ'&project_name=='RPHs') %>% 
  select(treatment, plot_id, block) %>%
  unique()

#experiments with blocking
blocked<-newdat %>%
  filter(project_name=='MAT2'|project_name=='change'|project_name=='RPHs')

scplistB<-unique(blocked$spc)

blocked_diff<-data.frame()

for (i in 1:length(scplistB)){

  site_proj=scplistB[i]

  subset<-blocked %>%
    filter(spc==scplistB[i])

  diff<-RAC_difference(subset, time.var='treatment_year', abundance.var = 'abundance', replicate.var = 'plot_id', species.var = 'genus_species', treatment.var = 'treat2', reference.treatment = 'control', block.var = 'block') %>%
    mutate(spc=site_proj)

  blocked_diff <-blocked_diff %>%
    bind_rows(diff)
}

CT_diff1<-blocked_diff %>%
  group_by(spc, treat22, treatment_year) %>%
  summarise_at(vars(richness_diff:species_diff), mean, na.rm=T)

#experiments that are not blocked
Notblocked<-newdat %>%
  filter(project_name!='MAT2'&project_name!='change'&project_name!='RPHs')
  
scplist<-unique(newdat$spc)

notblocked_diff<-data.frame()

for (i in 1:length(scplist)){
  
  site_proj=scplist[i]
  
  subset<-newdat %>% 
    filter(spc==scplist[i])
  
  diff<-RAC_difference(subset, time.var='treatment_year', abundance.var = 'abundance', replicate.var = 'plot_id2', species.var = 'genus_species', treatment.var = 'treat2', reference.treatment = 'control') %>% 
    mutate(spc=site_proj)
  
  notblocked_diff <-notblocked_diff %>% 
    bind_rows(diff)
}

CT_diff2<-notblocked_diff %>% 
  group_by(spc, plot_id2, treat22, treatment_year) %>% 
  summarise_at(vars(richness_diff:species_diff), mean, na.rm=T) %>% 
  group_by(spc, treat22, treatment_year) %>% 
  summarise_at(vars(richness_diff:species_diff), mean, na.rm=T)

CT_diff <- CT_diff1 %>% 
  bind_rows(CT_diff2) %>% 
  separate(spc, into=c('site_code', 'project_name', 'community_type'), sep="::", remove = F) %>% 
  mutate(facet=paste(site_code, treat22, sep="_"))

ggplot(data=CT_diff, aes(x=treatment_year, y=richness_diff, color=spc))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = 0)+
  facet_wrap(~facet)


test<-subset %>% 
  select(block, treat2, plot_id) %>% 
  unique()
group_by(spc, plot_id, treat22, treatment_year) %>% 
  summarise_at(vars(richness_diff:species_diff), mean, na.rm=T)

##productivity responses
subsetdatasets<-info %>% 
  select(site_code, project_name, community_type, spc, treatment, plot_id, treat2, trt_type) %>% 
  unique()

prodSubset<-proddat %>% 
  right_join(subsetdatasets) %>% 
  drop_na() %>% 
  filter(site_code!='NWT'&site_code!='DL')

control_prod<-prodSubset %>% 
  group_by(site_code, project_name, community_type, treatment, treat2, trt_type, treatment_year) %>% 
  summarise(contprod=mean(anpp)) %>% 
  filter(treat2=='control') %>% 
  ungroup() %>% 
  select(site_code, project_name, community_type, contprod, treatment_year)
  

proddiff<-prodSubset %>% 
  group_by(site_code, project_name, community_type, spc, treatment, treat2, trt_type, treatment_year) %>% 
  filter(treat2!='control') %>% 
  summarise(treatprod=mean(anpp)) %>% 
  left_join(control_prod) %>% 
  mutate(prod_diff=(treatprod-contprod)/contprod)%>% 
  mutate(facet=paste(site_code, treat2, sep="_"))

ggplot(data=proddiff, aes(x=treatment_year, y=prod_diff, color=spc))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = 0)+
  facet_wrap(~facet, scales='free')
