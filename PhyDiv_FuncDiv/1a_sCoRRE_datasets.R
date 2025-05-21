################################################################################
##  1a_sCoRRE_datasets.R: Selecting the appropriate datasets from CoRRE Database for analysis.
##
##  Author: Kimberly Komatsu
##  Date created: May 21, 2025
################################################################################

library(tidyverse)

#### read data ####

#spp names
names <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data\\trait data\\corre2trykey_2021.csv') %>% 
  select(genus_species, species_matched) %>% 
  unique()

#trait data
contTraitSpp <- read.csv('https://pasta.lternet.edu/package/data/eml/edi/1533/3/169fc12d10ac20b0e504f8d5ca0b8ee8') %>% 
  select(family, species) %>% 
  unique()

#trt data - subset to relevant treatments
trt_analysis <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\RawAbundanceMarch2024.csv') %>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id) %>%
  unique() %>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\ExperimentInfo_March2024.csv')) %>%
  group_by(site_code, project_name, community_type) %>%
  mutate(experiment_length=max(treatment_year)) %>%
  ungroup() %>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, 
         trt_type, experiment_length, plot_mani, n, p, CO2, precip, temp) %>%
  mutate(alltrts=ifelse(trt_type %in% c("control", "CO2","CO2*temp", "mow_clip","burn","burn*graze","disturbance",
                                        "burn*mow_clip","drought","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip",
                                        "herb_removal","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip",
                                        "irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip",
                                        "N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N","mult_nutrient","N*P","P",
                                        "N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze",
                                        "P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp",
                                        "N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip",
                                        "N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal",
                                        "mult_nutrient*herb_removal*mow_clip","temp","temp*mow_clip","drought*temp","irr*temp","irr"),1,0)) %>%
  filter(alltrts==1) %>%
  mutate(dist=ifelse(trt_type %in% c("mow_clip","burn","burn*graze","disturbance","burn*mow_clip"), 1, 0), #unify codes across datasets
         # tCO2=ifelse(trt_type %in% c("CO2"), 1, 0),
         drought=ifelse(trt_type %in% c("drought"), 1, 0),
         # therb_removal=ifelse(trt_type %in% c("herb_removal"), 1, 0),
         irg=ifelse(trt_type %in% c("irr"), 1, 0),
         # ttemp=ifelse(trt_type %in% c("temp"), 1, 0),
         # tn=ifelse(trt_type %in% c("N"), 1, 0),
         # tp=ifelse(trt_type %in% c("P"), 1, 0),
         multtrts=ifelse(trt_type %in% c("CO2*temp", "burn*graze","burn*mow_clip","drought*CO2*temp","drought*mow_clip",
                                         "drought*temp*mow_clip","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip",
                                         "irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip",
                                         "N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N*CO2","N*mow_clip","N*burn",
                                         "N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought",
                                         "N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp",
                                         "N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn",
                                         "P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip",
                                         "temp*mow_clip","drought*temp","irr*temp","mult_nutrient","N*P"),1,0)) %>%
  mutate(trt_type2=ifelse(dist==1, 'disturbance', ifelse(multtrts==1, 'multiple trts', trt_type))) %>%
  select(site_code, project_name, community_type, treatment, alltrts, dist, drought, irg, multtrts, trt_type2, plot_mani) %>%
  unique()

#community data
comm2 <- read.table("C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\RelativeCoverMarch2024.csv", header=T, sep=",", fill = TRUE) %>% 
  right_join(trt_analysis) %>%
  mutate(trt_binary=ifelse(plot_mani>0, 1, 0)) %>% 
  mutate(drop=ifelse(site_code=="CDR" & treatment %in% c(2, 3, 4, 5, 7),
                     1, 0)) %>% #drop some of the CDR e001 and e002 treatments to prevent over-representation
  filter(drop==0) %>%
  filter(!(project_name %in% c('gap'))) #remove pulse light expt

#total cover
totCover <- comm2 %>% 
  left_join(names) %>% 
  filter(!is.na(species_matched)) %>% #remove species identified only to genus level
  select(-genus_species) %>% 
  rename(species=species_matched) %>%
  left_join(contTraitSpp) %>% 
  filter(!is.na(family)) %>% #remove species without continuous trait data
  group_by(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block, plot_id) %>% 
  summarize(totcov=sum(relcov), .groups='drop') #total cover of remaining species (7569 of plots are <0.8 totcov, will be dropped; 23.4% of plots)

hist(totCover$totcov)

comm <- comm2 %>% 
  left_join(names) %>% 
  filter(!is.na(species_matched)) %>% #remove species identified only to genus level (drops 22,966 data points, 6.4% of data)
  select(-genus_species) %>% 
  rename(species=species_matched) %>% 
  left_join(contTraitSpp) %>% 
  filter(!is.na(family)) %>% #remove species without continuous trait data (drops 38,288 data points, 11.4% of data)
  left_join(totCover) %>% 
  mutate(plot_id=ifelse(project_name=='NSFC', paste(plot_id, treatment, sep='__'), plot_id)) %>% #rename plot id for NSFC expt
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block, plot_id, species, relcov, totcov) %>% 
  mutate(plot_id2 = paste(site_code, project_name, community_type, calendar_year, plot_id, sep = "::")) %>%  #create new plot identifier
  group_by(plot_id2) %>% 
  mutate(richness=length(species)) %>% 
  ungroup() %>% 
  filter(richness>1, #remove plots with only one species (drops 700 data points, 0.1% of data)
         totcov>0.8) #remove plots with less than 80% cover of species with known trait values (drops 59,7749 data points, 20.2% of data points)

#lists
spp <- comm %>% select(species) %>% unique() #1559 spp
plots <- comm %>% select(site_code, project_name, community_type, plot_id) %>% unique() #3408 plots
expt <- comm %>% select(site_code, project_name, community_type) %>% unique() #125 experiments
sites <- unique(comm$site_code) #66 sites


#### save data ####

# saveRDS(comm, file = "PhyDiv_FuncDiv/PD_FD_commData.rds") # saving subset of CoRRE database for analysis
# saveRDS(spp, file = "PhyDiv_FuncDiv/sppList.rds") # saving species list for analysis