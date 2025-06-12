################################################################################
##  PhyloDiv_metrics.R: Calculating phylogenetic diversity metrics for the CoRRE database.
##
##  Author: Kimberly Komatsu, Josep Padulles Cubino
##  Date created: December 10, 2019
################################################################################

#load packages:
library(EDIutils)
library(rlist)
library(matrixStats)
library(picante)
library(tidyverse)


#### read data ####

#spp names
names <- readRDS('PhyDiv_FuncDiv/sppList.rds')

#community data
comm <- readRDS('PhyDiv_FuncDiv/PD_FD_commData.rds')

#treatment data
trt <- readRDS('PhyDiv_FuncDiv/trt_info.rds')


#### calculate phylogenetic diversity metrics using one single tree (scenario 3) ####

#load tree
scorre.tree <- read.tree("PhyDiv_FuncDiv/corre.phylo.tree.S3_20250521.tre")

sites <- unique(comm$site_code) #66 sites

phylogeneticDiversityMetrics <- NULL

for (i in 1:length(sites)){ #loop to calculate metrics for each site independently

  print(i*100/length(sites))
  
  comm2 <- comm %>% 
    filter(site_code == sites[i]) %>%  #subset plots within each site
    filter(treatment_year>0) %>% #only keep treatment data
    select(plot_id2, species, relcov) %>% 
    # mutate(relcov=ifelse(relcov>0, 1, 0)) %>% 
    unique() %>% 
    group_by(plot_id2, species) %>% 
    summarize(relcov=mean(relcov)) %>% 
    ungroup() %>% 
    pivot_wider(names_from=species, values_from=relcov, values_fill=0) #make species matrix
    
  colnames(comm2) <- gsub(" ", "_", colnames(comm2)) #add underscore in column names
  
  #Prune tree with only species in our site:
  tree <- keep.tip(scorre.tree, colnames(comm2[,-1]))
  
  #distance matrix
  distance <- as.matrix(cophenetic(tree))
  
  #species matrix
  spp <- comm2 %>% 
    column_to_rownames('plot_id2') %>% 
    as.matrix(.)
  
  #calculate abundance-weighted phylogenetic diversity metrics:
  mpd.raw <- as.data.frame(mpd(samp=spp, dis=distance,  abundance.weighted = T))
  
  pd.raw <- cbind(comm2[,1], mpd.raw) %>% 
    rename(mpd="mpd(samp = spp, dis = distance, abundance.weighted = T)") %>% 
    mutate(permutation=0)
  
  pd.ses <- {}
  sesVector <- c(1:999)
  for(n in 1:length(sesVector)){
    distanceShuffle <- as.data.frame(distance) %>% 
      rownames_to_column(var = "identifier") %>% 
      mutate(spp_shuffling = sample(identifier, size = n(), replace = FALSE)) %>% 
      select(-identifier) %>% 
      column_to_rownames(var="spp_shuffling") %>% 
      as.matrix(.)
    
    mpd.ses <- as.data.frame(mpd(samp=spp, dis=distanceShuffle,  abundance.weighted = T))
    
    ses <- cbind(comm2[,1], mpd.ses) %>% 
      rename(mpd="mpd(samp = spp, dis = distanceShuffle, abundance.weighted = T)") %>% 
      mutate(permutation=sesVector[n])
    
    pd.ses <- rbind(pd.ses, ses)
    
    rm(list=ls()[grep("SES", ls())])
  }
  
  pd.all <- rbind(pd.raw, pd.ses)
  
  #mean MPD and MNTD for null distribution to create SES MPD
  sesSubsetMean <- pd.all %>% 
    separate(plot_id2, into=c('site_code', 'project_name', 'community_type', 'calendar_year', 'plot_id'), sep='::') %>% 
    mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::')) %>% 
    mutate(plot_id2=paste(site_proj_comm, calendar_year, plot_id, sep='::')) %>% 
    mutate(calendar_year=as.integer(calendar_year)) %>% 
    group_by(plot_id2, site_proj_comm, calendar_year, plot_id) %>% 
    summarise(across(c('mpd'), list(mean=mean, sd=sd))) %>% 
    ungroup()
  
  #calculate SES values for MPD 
  mpdSubsetSES <- pd.all %>% 
    full_join(sesSubsetMean) %>% 
    filter(permutation==0) %>% 
    mutate(MPD_phylo_ses=(mpd-mpd_mean)/mpd_sd) %>%
    select(plot_id2, site_proj_comm, calendar_year, plot_id, mpd, MPD_phylo_ses) %>% 
    rename(MPD_phylo=mpd)
  
  #bind values into RR dataframe
  phylogeneticDiversityMetrics <- rbind(phylogeneticDiversityMetrics, mpdSubsetSES)

}

#save output:
saveRDS(phylogeneticDiversityMetrics, file = "PhyDiv_FuncDiv/phylogeneticDiversityMetrics.rds")
