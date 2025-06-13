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


##### functions and themes #####
###standard error function
se <- function(x, na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}

##### data #####

#spp names
names <- readRDS('PhyDiv_FuncDiv/sppList.rds')

#community data
comm <- readRDS('PhyDiv_FuncDiv/PD_FD_comm.rds') %>% 
  mutate(replicate=paste(site_code, project_name, community_type, treatment, plot_id, sep='::'))

#treatment data
trt <- readRDS('PhyDiv_FuncDiv/trt_info.rds')


##### calculate diversity metrics ##### 
#getting community diversity metrics for each plot
richness <- community_structure(comm, time.var="calendar_year", abundance.var="relcov", replicate.var="replicate") %>%
  separate(replicate, into=c("site_code", "project_name", "community_type", "treatment", "plot_id"), sep='::')

#hill numbers
sppMatrix <- comm %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::')) %>% 
  select(site_proj_comm, treatment, calendar_year, plot_id, species, relcov) %>% 
  filter(species!='')

label <- sppMatrix %>%
  select(site_proj_comm) %>%
  unique()


#calculate richness for each site
hillNumbers <- data.frame(row.names=1) 

for(i in 1:length(label$site_proj_comm)) {
  subset <- sppMatrix[sppMatrix$site_proj_comm==as.character(label$site_proj_comm[i]),] %>% 
            pivot_wider(names_from=species, values_from=relcov, values_fill=0)
  
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


#save output:
saveRDS(rDiv, file = "PhyDiv_FuncDiv/rDiv.rds")










