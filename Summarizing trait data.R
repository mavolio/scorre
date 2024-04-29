#Summarizing trait data

setwd('C://Users//mavolio2//Dropbox//CoRRE_database//Data//CleanedData//Traits//')

library(tidyverse)

#error family < 3 drops ~20000 observations
Data_error<-read.csv('CoRRE_allTraitData_Oct2023.csv') %>% 
  filter(error_risk_overall<3|is.na(error_risk_overall)) %>%
  mutate(trait_value2=as.numeric(trait_value))

mean<-Data_error %>% 
  group_by(species_matched, trait) %>% 
  summarize(mean=mean(trait_value2)) %>% 
  group_by(trait) %>% 
  summarise(n=length(mean))

categorical<-Data_error %>% 
  filter(trait_value!="uncertain") %>% 
  group_by(trait) %>% 
  summarise(n=length(trait_value))

setwd('C:\\Users\\mavolio2\\Dropbox\\CoRRE_database\\Data\\CompiledData')

info <- read.csv('ExperimentInfo.csv') %>% 
  filter(trt_type!="control") %>% 
  mutate(Trt_type2=ifelse(trt_type %in% c("CO2"), "CO2", 
                   ifelse(trt_type %in% c("drought"), "drought", 
                   ifelse(trt_type %in% c("irr"), "irg", 
                   ifelse(trt_type %in% c("temp"), "temp",
                   ifelse(trt_type=="N", "N", 
                   ifelse(trt_type=="P", "P", 
                  ifelse(trt_type %in% c("CO2*temp", "burn*graze","burn*mow_clip","drought*mow_clip","drought*temp*mow_clip","herb_removal*mow_clip","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","P*burn","P*mow_clip","temp*mow_clip","drought*temp","irr*temp",'C', 'C*stone', 'irr*plant_mani', 'irr*plant_mani*herb_removal',  'N*fungicide', 'N*plant_mani', 'N*plant_mani*disturbance','N*plant_mani*mow_clip', 'N*stone','N*temp*fungicide', 'P*plant_mani', 'precip_vari*temp'), "mult treats",
                  ifelse(trt_type %in% c("mult_nutrient","N*P","drought*CO2*temp","irr*CO2","irr*CO2*temp","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N*CO2","N*drought","N*P*burn","N*P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip",'N*P*burn*mow_clip', 'mult_nutrient*drought', 'N*P*burn*mow_clip', 'N*P*plant_mani','mult_nutrient*fungicide', 'mult_nutrient*plant_mani', 'mult_nutrient*plant_mani*herb_removal'), "mult_resources",
                  ifelse(trt_type %in% c("herb_removal", "mow_clip","burn","burn*graze","disturbance","burn*mow_clip", 'fungicide',  'lime', 'plant_mani', 'plant_mani*disturbance', 'plant_mani*herb_removal','plant_mani*mow_clip','stone', 'temp*fungicide'), "non-resource", 
                   ifelse(trt_type %in% c('precip_vari','K', 'light'), "other resource", "missed")))))))))))

anpp<-read.csv('ANPP2021.csv') %>% 
  left_join(info) %>% 
  filter(plot_mani!=0) %>%
  select(site_code, project_name, community_type, Trt_type2) %>% 
  unique() %>% 
  group_by(Trt_type2) %>% 
  summarize(n=length(site_code))
  
  
abundance<-read.csv("RawAbundance.csv") %>% 
  left_join(info) %>% 
  filter(plot_mani!=0) %>%
  select(site_code, project_name, community_type, Trt_type2) %>% 
  unique() %>% 
  group_by(Trt_type2) %>% 
  summarize(n=length(site_code))
