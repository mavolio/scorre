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
rep<-alldat%>%
  filter(community_type!=0) %>% 
  group_by(site_code, project_name, community_type, treatment) %>% 
  summarize(max=max(treatment_year), n=length(treatment_year))

###subuset out sites and experiments
arc<-alldat %>% 
  filter(site_code=='ARC')

#not enough data to use
asga<-alldat %>% 
  filter(site_code=='ASGA') %>% 
  filter(project_name=='clonal'&treatment=="mixed_CO"|project_name=='clonal'&treatment=='mixed_UN'|project_name=='Exp1'&treatment=='1_0_CO'|project_name=='Exp1'&treatment=='1_0_UN') %>% 
  mutate(trt2=ifelse(treatment=='mixed_CO'|treatment=='1_0_CO', 'control', 'Nit')) %>% 
  select(-treatment) %>% 
  rename(treatment=trt2)

#study too short
bt<-alldat %>% 
  filter(site_code=="Bt") %>% 
  filter(project_name=='NPKDNet'&treatment=="control"|project_name=='NPKDNet'&treatment=='NPK'|project_name=='NutNet'&treatment=='Control'|project_name=='NutNet'&treatment=='NPK') %>% 
  mutate(trt2=ifelse(treatment=='Control'|treatment=='control', 'control', 'Nit')) %>% 
  select(-treatment) %>% 
  rename(treatment=trt2)

test<-bt %>% 
  select(project_name, treatment, n, p, k, plot_mani) %>% 
  unique()

cdr<-alldat %>% 
  filter(site_code=="CDR") %>% 
  filter(project_name=='e001'|project_name=='e002') %>% 
  mutate(trt1=paste('N', n, sep="_")) %>% 
  mutate(trt2=ifelse(trt1=="N_0", 'control', trt1)) %>% 
  select(-treatment, -trt1) %>% 
  rename(treatment=trt2)

dl<-alldat %>% 
  filter(site_code=="DL") %>% 
  filter(treatment=='C'|treatment=='M'|treatment=='N'|treatment=='MN')

dl_nsfc<-dl %>% 
  filter(project_name=='NSFC') %>% 
  select(-block) %>% 
  mutate(block=plot_id) %>% 
  mutate(plot=paste(plot_id, treatment, sep="_")) %>% 
  select(-plot_id) %>% 
  rename(plot_id=plot)

dl2<-dl %>% 
  filter(project_name!='NSFC') %>% 
  bind_rows(dl_nsfc)

#studies too short
jrn<-alldat %>% 
  filter(site_code=='JRN') %>% 
  filter(treatment=='T'|treatment=='C'|treatment=='P3N0'|treatment=='P3N1') %>% 
  mutate(trt2=ifelse(treatment=='C'|treatment=='P3N0', 'control', 'Nit')) %>% 
  select(-treatment) %>% 
  rename(treatment=trt2)

knz<-alldat %>% 
  filter(site_code=='KNZ') %>% 
  filter(treatment=='u_u_c'|treatment==0|treatment=="N1P0"|treatment=='u_u_n'|treatment=='10'|treatment=='N2P0') %>% 
  mutate(trt2=ifelse(treatment=='u_u_c'|treatment=='N1P0'|treatment==0, 'control', 'Nit')) %>% 
  select(-treatment) %>% 
  rename(treatment=trt2)

nwt<-alldat %>% 
  filter(site_code=='NWT') %>% 
  filter(treatment=='XXX'|treatment=='Control'|treatment=="N"|treatment=='XNX') %>% 
  filter(n<11) %>% 
  mutate(trt2=ifelse(treatment=='XXX'|treatment=='Control', 'control', 'Nit')) %>% 
  select(-treatment) %>% 
  rename(treatment=trt2)

scl<-alldat %>% 
  filter(site_code=='SCL')%>% 
  filter(treatment=='N0'|treatment=='OO'|treatment=="N1"|treatment=='OF') %>% 
  mutate(trt2=ifelse(treatment=='N0'|treatment=='OO', 'control', 'Nit')) %>% 
  select(-treatment) %>% 
  rename(treatment=trt2)

#not good overlap in when data was collected - might be able to use for productivity
serc<-alldat %>% 
  filter(site_code=="SERC")%>% 
  filter(treatment=='A'|treatment=='E'|treatment=="t1"|treatment=='t3') %>% 
  mutate(trt2=ifelse(treatment=='A'|treatment=='t1', 'control', 'Elev')) %>% 
  select(-treatment) %>% 
  rename(treatment=trt2)


###make control always be control
newdat<-arc %>% 
  bind_rows(cdr, dl2, knz, nwt, scl) %>% 
  mutate(spc=paste(site_code, project_name, community_type, sep='::')) %>% 
  select(site_code, project_name, spc, calendar_year, treatment_year, treatment, plot_id, block, abundance, genus_species, plot_mani, n, p, k, CO2, mow_clip, successional, trt_type) %>% 
  mutate(treat2=ifelse(treatment %in% c('C', 'control', 'CT'), 'control', treatment))

info<-newdat %>% 
  select(site_code, project_name, spc, treatment_year, calendar_year, treatment, plot_id, block, plot_mani, n, p, k, CO2, mow_clip, successional, trt_type) %>% 
  unique()

list<-unique(newdat$spc)

table(newdat$treat2)

###seeing number of replicated treatments at each site
numexp<-newdat %>% 
  select(site_code, spc, treatment) %>% 
  unique() %>% 
  group_by(site_code, treatment) %>% 
  summarize(n=length(spc))

#there are 28 C-T comparisons across   


###Community outcomes

#how many have pre-treatment data? only 5, so not enough to use pretreatment data
pretrt<-newdat %>% 
  filter(treatment_year==0) %>% 
  select(site_code, spc) %>% 
  unique()

#doing control-treatment differences over time taking into account blocking when possible

#experiments with blocking
blocked<-newdat %>% 
  filter(site_code %in% c('ASGA')|site_code=='ARC'&project_name=='MAT2'|site_code=='KNZ'&project_name=='change'|site_code=="DL"&project_name=='NSFC')

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
  filter(site_code %in% c('Bt', 'CDR', 'DL', 'JRN', 'NWT', 'SERC')|site_code=='ARC'&project_name=='MNT'|site_code=='KNZ'&project_name=='pplots'|site_code=='KNZ'&project_name=='BGP')
  
scplistNB<-unique(Notblocked$spc)

notblocked_diff<-data.frame()

for (i in 1:length(scplistNB)){
  
  site_proj=scplistNB[i]
  
  subset<-Notblocked %>% 
    filter(spc==scplistNB[i])
  
  diff<-RAC_difference(subset, time.var='treatment_year', abundance.var = 'abundance', replicate.var = 'plot_id', species.var = 'genus_species', treatment.var = 'treat2', reference.treatment = 'control') %>% 
    mutate(spc=site_proj)
  
  notblocked_diff <-notblocked_diff %>% 
    bind_rows(diff)
}

CT_diff2<-notblocked_diff %>% 
  group_by(spc, plot_id, treat22, treatment_year) %>% 
  summarise_at(vars(richness_diff:species_diff), mean, na.rm=T) %>% 
  group_by(spc, treat22, treatment_year) %>% 
  summarise_at(vars(richness_diff:species_diff), mean, na.rm=T)

CT_diff <- CT_diff1 %>% 
  bind_rows(CT_diff2) %>% 
  separate(spc, into=c('site_code', 'project_name', 'community_type'), sep="::", remove = F) %>% 
  mutate(facet=paste(site_code, treat22, sep="_"))

ggplot(data=CT_diff, aes(x=treatment_year, y=species_diff, color=spc))+
  geom_point()+
  geom_line()+
  facet_wrap(~facet)


test<-subset %>% 
  select(block, treat2, plot_id) %>% 
  unique()
group_by(spc, plot_id, treat22, treatment_year) %>% 
  summarise_at(vars(richness_diff:species_diff), mean, na.rm=T)

