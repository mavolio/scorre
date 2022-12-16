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


setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\paper 2_PD and FD responses\\data\\')  #kim's computer

##### functions and themes #####
###standard error function
se <- function(x, na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}

##### data #####
#treatment data
trt <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\RawAbundance.csv') %>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id) %>%
  unique() %>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\ExperimentInfo.csv')) %>%
  group_by(site_code, project_name, community_type) %>%
  mutate(experiment_length=max(treatment_year)) %>%
  ungroup() %>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, trt_type, experiment_length, plot_mani, n, p, CO2, precip, temp)

#species relative cover data
relCover <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeCover_Dec2021.csv') %>%
  mutate(replicate=paste(site_code, project_name, community_type, plot_id, sep='::')) #creating identifying column of each plot

#getting community diversity metrics for each plot
richness <- community_structure(relCover, time.var="treatment_year", abundance.var="relcov", replicate.var="replicate") %>%
  separate(replicate, into=c("site_code", "project_name", "community_type", "plot_id"), sep='::')

#hill numbers
sppMatrix <- relCover %>% 
  mutate(site_proj_comm_year_plot=paste(site_code, project_name, community_type, treatment_year, plot_id, sep='::')) %>% 
  select(site_proj_comm_year_plot, genus_species, relcov) %>% 
  pivot_wider(names_from=genus_species, values_from=relcov, values_fill=0)
#START HERE fix fill with null problem

hill <- hill_taxa(sppMatrix, q = 1, MARGIN = 1, base = exp(1))


#combine taxonomic diversity
rDiv <- richness # %>% left_join(hill)

# write.csv(rDiv, 'CoRRE_taxonomicDiversity_2022-12-15.csv', row.names=F)











