library(tidyverse)#

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data')
setwd('C:\\Users\\megha\\Dropbox\\sDiv_sCoRRE_shared\\CoRRE data')

correProjectInfo <- read.csv('CoRRE_project_summary.csv')
correTrtInfo <- read.csv('CoRRE_treatment_summary.csv')
correRelAbund <- read.csv('CoRRE_relative_abundance_Nov2019.csv')%>%
  select(-X)
correCleanNames<-read.csv("sCoRRE_tnrs_matching.csv")%>%
  select(-X, -AccSpeciesName)%>%
  rename(genus_species=genus.species)%>%
  unique()
  

#number of total plots, species, observations (plot*years)
all <- correRelAbund%>%
  left_join(correTrtInfo)%>%
  left_join(correProjectInfo)
plots <- all%>%
  select(site_code, project_name, community_type, plot_id)%>%
  unique()
years <- all%>%
  select(site_code, project_name, community_type, plot_id, calendar_year)%>%
  unique()


#histogram of experiment lengths
ggplot(data=correProjectInfo, aes(x=experiment_length)) +
  geom_histogram() + xlab('Experiment Length (yrs)') + ylab('Count')
summary(correProjectInfo$experiment_length)

#list of treatment types and number of experiments for which they are present
trtTypes <- correTrtInfo%>%
  select(site_code, project_name, community_type, trt_type)%>%
  unique()%>%
  filter(trt_type!="control")%>%
  group_by(trt_type)%>%
  summarize(n=length(trt_type))

write.csv(trtTypes,"trt_type.csv", row.names=F)

#list of species and how often they occur in the dataset
species<-correRelAbund%>%
  left_join(correCleanNames)%>%
  filter(!is.na(Name_matched))%>%
  mutate(site_proj_comm=paste(site_code, project_name, community_type))

cont_dat<-species%>%
  select(site_proj_comm, Name_matched, calendar_year)%>%
  unique()%>%
  group_by(site_proj_comm, Name_matched)%>%
  summarize(n=length(calendar_year))%>%
  filter(n!=1)

length(unique(cont_dat$Name_matched))

species_occurances_trt<-species%>%
  select(site_proj_comm, treatment, Name_matched)%>%
  unique()%>%
  group_by(site_proj_comm, Name_matched)%>%
  summarize(num_trts=length(treatment))

write.csv(species_occurances_trt, "basic dataset info/species_occurances_treatments.csv", row.names=F)
  
species_occurances_sites<-species%>%
  select(site_code, Name_matched)%>%
  unique()%>%
  group_by(Name_matched)%>%
  summarize(num_sites=length(site_code))

write.csv(species_occurances_sites, "basic dataset info/species_occurances_sites.csv", row.names=F)
