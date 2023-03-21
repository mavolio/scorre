################################################################################
##  sCoRRE_taxonomicDiversityMetrics.R: Examining differences in taxonomic diveristy within the CoRRE database.
##
##  Author: Kimberly Komatsu
##  Date created: December 16, 2022
################################################################################

library(data.table)
library(codyn)
library(hillR)
library(tidyverse)


setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared')  #kim's computer

##### functions and themes #####
###standard error function
se <- function(x, na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}

##### data #####
#treatment data
trt <- read.csv('CoRRE data\\CoRRE data\\community composition\\CoRRE_RawAbundance_Jan2023.csv') %>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id) %>%
  unique() %>%
  left_join(read.csv('CoRRE data\\CoRRE data\\basic dataset info\\ExperimentInfo.csv')) %>%
  group_by(site_code, project_name, community_type) %>%
  mutate(experiment_length=max(treatment_year)) %>%
  ungroup() %>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, trt_type, experiment_length, plot_mani, n, p, CO2, precip, temp)

#species relative cover data
relCover <- read.csv('CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeCover_Jan2023.csv') %>%
  mutate(replicate=paste(site_code, project_name, community_type, treatment, plot_id, sep='::')) #creating identifying column of each plot


##### calculate diversity metrics ##### 
#getting community diversity metrics for each plot
richness <- community_structure(relCover, time.var="treatment_year", abundance.var="relcov", replicate.var="replicate") %>%
  separate(replicate, into=c("site_code", "project_name", "community_type", "treatment", "plot_id"), sep='::')

#hill numbers
sppMatrix <- relCover %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::')) %>% 
  select(site_proj_comm, treatment, treatment_year, plot_id, genus_species, relcov) %>% 
  group_by(site_proj_comm, treatment, treatment_year, plot_id, genus_species) %>% 
  summarise(relcov=mean(relcov)) %>% #average for CHY EDGE, which has multiple values per plot
  ungroup() %>% 
  filter(genus_species!='')

label <- sppMatrix %>%
  select(site_proj_comm) %>%
  unique()

hillNumbers <- data.frame(row.names=1) 

#calculate richness for each site
for(i in 1:length(label$site_proj_comm)) {
  subset <- sppMatrix[sppMatrix$site_proj_comm==as.character(label$site_proj_comm[i]),] %>% 
            pivot_wider(names_from=genus_species, values_from=relcov, values_fill=0)
  
  hill <- hill_taxa(subset[,-1:-4], q = 1, MARGIN = 1, base = exp(1))
  
  info <- subset[,1:4] %>% 
    bind_cols(hill) %>% 
    rename(hill=...5)

  hillNumbers <- rbind(hillNumbers, info)  
}

hillNumbers2 <- hillNumbers %>% 
  separate(site_proj_comm, into=c('site_code', 'project_name', 'community_type'), sep='::')


##### combine taxonomic diversity #####
rDiv <- richness %>% 
  left_join(hillNumbers2)

# write.csv(rDiv, 'CoRRE_taxonomicDiversity_2023-03-21.csv', row.names=F)











