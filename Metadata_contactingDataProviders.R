library(readxl)
library(tidyverse)

#contact info
contact<-read_xlsx('C:\\Users\\mavolio2\\Dropbox\\CoRRE_database\\Data\\Data_Summaries\\Data_log.xlsx') %>% 
  rename(update='2020_update')

contact2<-contact%>% 
  filter(update %in% c('new_site', 'original_corre'))

#experiment details
details<-read_xlsx('C:\\Users\\mavolio2\\Dropbox\\CoRRE_database\\Data\\Data_Summaries\\Experiment_List.xlsx') %>% 
  rename(treats='what they did',
         length='experiment length in years (not including pre-treatment)',
         cite='citation (all papers can be found in Publications of Experiments folder)')

treatments<-read.csv('C:\\Users\\mavolio2\\Dropbox\\CoRRE_database\\Data\\CompiledData\\ExperimentInfo_March2024.csv') %>% 
  select(site_code, project_name, trt_type, n:plant_trt) %>% 
  unique() %>% 
  filter(trt_type!='control') 

#site location
loc<-read.csv('C:\\Users\\mavolio2\\Dropbox\\CoRRE_database\\Data\\CompiledData\\SiteLocationClimate.csv')

all<-contact2 %>% 
  left_join(details) %>% 
  left_join(loc) %>% 
  select(site_code, Location, project_name, treats, length, cite, Contact, 'co-authors')

write.csv(all, 'C:\\Users\\mavolio2\\Dropbox\\CoRRE_database\\Data\\CompiledData\\ForCoAuthors.csv', row.names=F)
