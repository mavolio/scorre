################################################################################
##  PhyloDiv_metrics.R: Calculating phylogenetic diversity metrics for the CoRRE database.
##
##  Author: Kimberly Komatsu, Josep Padulles Cubino
##  Date created: December 10, 2019
################################################################################

#load packages:
library(rlist)
library(matrixStats)
library(picante)
library(tidyverse)

#### set directory ####
setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data\\')  #kim's computer
my.wd <- "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/" #padu


#### read data ####

#spp names
names <- read.csv('trait data\\corre2trykey_2021.csv') %>% 
  select(genus_species, species_matched) %>% 
  unique()

#community data
comm <- read.table("CoRRE data\\community composition\\CoRRE_RelativeCover_Jan2023.csv", header=T, sep=",", fill = TRUE) %>% 
  left_join(names) %>% 
  left_join(read.csv('trait data\\sCoRRE categorical trait data_12142022.csv')) %>% 
  filter(growth_form!='moss' & leaf_type!='microphyll') %>% 
  mutate(plot_id=ifelse(project_name=='NSFC', paste(plot_id, treatment, sep='__'), plot_id)) %>% 
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block, plot_id, genus_species, species_matched, relcov) %>% 
  filter(!grepl("(sp.)$", species_matched)) %>%  #remove species identified only at the genus level
  mutate(plot_id2 = paste(site_code, project_name, community_type, calendar_year, plot_id, sep = "::")) %>%  #create new plot identifier
  group_by(plot_id2) %>% 
  mutate(richness=length(genus_species)) %>% 
  ungroup() %>% 
  filter(richness>1) #remove plots with only one species


#species list
spp <- comm %>% 
  select(genus_species, species_matched) %>% 
  unique()


#create list of sites
sites <- unique(comm$site_code)

#treatment data
trt <- read.csv('CoRRE data\\community composition\\CoRRE_RawAbundance_Jan2023.csv') %>%
  select(site_code, project_name, community_type, calendar_year, calendar_year, treatment, plot_id) %>%
  unique() %>%
  left_join(read.csv('CoRRE data\\basic dataset info\\ExperimentInfo.csv')) %>%
  group_by(site_code, project_name, community_type) %>%
  mutate(experiment_length=max(calendar_year)) %>%
  ungroup() %>%
  select(site_code, project_name, community_type, calendar_year, calendar_year, treatment, plot_id, trt_type, experiment_length, plot_mani, n, p, CO2, precip, temp) %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::')) %>% 
  mutate(plot_id2=paste(site_proj_comm, calendar_year, plot_id, sep='::')) %>% 
  mutate(plot_id=ifelse(project_name=='NSFC', paste(plot_id, treatment, sep='__'), plot_id)) %>% 
  unique()


#### calculate phylogenetic diversity metrics using one single tree (scenario 3) (non-weighted by abundances) ####

#load tree
scorre.tree <- read.tree("Phylogenies\\scorre.phylo.tree.S3_20230427.tre")

phylogeneticDiversityMetrics <- NULL

for (i in 1:length(sites)){ #loop to calculate metrics for each site independently

  print(i*100/length(sites))
  
  comm2 <- comm %>% 
    filter(site_code == sites[i]) %>%  #subset plots within each site
    filter(treatment_year>0) %>% #only keep treatment data
    select(plot_id2, species_matched, relcov) %>% 
    mutate(relcov=ifelse(relcov>0, 1, 0)) %>% 
    unique() %>% 
    pivot_wider(names_from=species_matched, values_from=relcov, values_fill=0) #make species matrix
    
  colnames(comm2) <- gsub(" ", "_", colnames(comm2)) #add underscore in column names
  
  #Prune tree with only species in our site:
  tree <- keep.tip(scorre.tree, colnames(comm2[,-1]))
  
  #distance matrix
  distance <- as.data.frame(cophenetic(tree))
  
  #species matrix
  spp <- comm2 %>% 
    column_to_rownames('plot_id2')
  
  #calculate phylogenetic diversity metrics:
  mpd.raw <- as.data.frame(mpd(samp=spp, dis=distance,  abundance.weighted = F))
  mntd.raw <- as.data.frame(mntd(samp=spp, dis=distance,  abundance.weighted = F))
  
  pd.raw <- cbind(comm2[,1], mpd.raw, mntd.raw) %>% 
    rename(mpd="mpd(samp = spp, dis = distance, abundance.weighted = F)",
           mntd="mntd(samp = spp, dis = distance, abundance.weighted = F)") %>% 
    mutate(permutation=0)
  
  pd.ses <- {}
  sesVector <- c(1:999)
  for(n in 1:length(sesVector)){
    distanceShuffle <- distance %>% 
      rownames_to_column(var = "identifier") %>% 
      mutate(spp_shuffling = sample(identifier, size = n(), replace = FALSE)) %>% 
      select(-identifier) %>% 
      column_to_rownames(var="spp_shuffling")
    
    mpd.ses <- as.data.frame(mpd(samp=spp, dis=distanceShuffle,  abundance.weighted = F))
    mntd.ses <- as.data.frame(mntd(samp=spp, dis=distanceShuffle,  abundance.weighted = F))
    
    ses <- cbind(comm2[,1], mpd.ses, mntd.ses) %>% 
      rename(mpd="mpd(samp = spp, dis = distanceShuffle, abundance.weighted = F)",
             mntd="mntd(samp = spp, dis = distanceShuffle, abundance.weighted = F)") %>% 
      mutate(permutation=sesVector[n])
    
    pd.ses <- rbind(pd.ses, ses)
    
    rm(list=ls()[grep("SES", ls())])
  }
  
  pd.all <- rbind(pd.raw, pd.ses)
  
  #mean MPD and MNTD for null distribution to create SES MPD and MNTD
  sesSubsetMean <- pd.all %>% 
    separate(plot_id2, into=c('site_code', 'project_name', 'community_type', 'calendar_year', 'plot_id'), sep='::') %>% 
    mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::')) %>% 
    mutate(plot_id2=paste(site_proj_comm, calendar_year, plot_id, sep='::')) %>% 
    mutate(calendar_year=as.integer(calendar_year)) %>% 
    group_by(plot_id2, site_proj_comm, calendar_year, plot_id) %>% 
    summarise(across(c('mpd', 'mntd'), list(mean=mean, sd=sd))) %>% 
    ungroup()
  
  #calculate SES values for MPD and MNTD
  mpdMNTDSubsetSES <- pd.all %>% 
    full_join(sesSubsetMean) %>% 
    mutate(MNTD_phylo_ses=(mntd-mntd_mean)/mntd_sd,
           MPD_phylo_ses=(mpd-mpd_mean)/mpd_sd) %>%
    select(plot_id2, site_proj_comm, calendar_year, plot_id, mpd, mntd, MPD_phylo_ses, MNTD_phylo_ses) %>% 
    rename(MPD_phylo=mpd,
           MNTD_phylo=mntd)
  
  #subset control plots
  mpdMNTDSubsetSESctl <- mpdMNTDSubsetSES %>% 
    left_join(trt) %>% 
    filter(plot_mani==0) %>% 
    group_by(site_proj_comm, calendar_year) %>% 
    summarise(across(c('MNTD_phylo', 'MPD_phylo', 'MNTD_phylo_ses', 'MPD_phylo_ses'), list(mean=mean))) %>% 
    ungroup()
  
  #calculate lnRR for ses values of MPD and MNTD
  mpdMNTDSubsetRRses <- mpdMNTDSubsetSES %>% 
    left_join(trt) %>% 
    full_join(mpdMNTDSubsetSESctl) %>% 
    mutate(RR_MNTD_phylo=ifelse(plot_mani>0, log(MNTD_phylo/MNTD_phylo_mean), NA),
           RR_MPD_phylo=ifelse(plot_mani>0, log(MPD_phylo/MPD_phylo_mean), NA),
           RR_MNTD_phylo_ses=ifelse(plot_mani>0, ((MNTD_phylo_ses-MNTD_phylo_ses_mean)/MNTD_phylo_ses_mean), NA), #percent difference for ses due to neg values
           RR_MPD_phylo_ses=ifelse(plot_mani>0, ((MPD_phylo_ses-MPD_phylo_ses_mean)/MPD_phylo_ses_mean), NA)) %>% #percent difference for ses due to neg values
    select(plot_id2, site_proj_comm, calendar_year, plot_id, MNTD_phylo, MNTD_phylo_ses, MPD_phylo, 
           MPD_phylo_ses, RR_MNTD_phylo, RR_MPD_phylo, RR_MNTD_phylo_ses, RR_MPD_phylo_ses)
  
  #lnRR MPD and MNTD for null distribution to create SES lnRR MPD and SES lnRR MNTD
  mpdMNTDSubsetPermCtl <- pd.all %>% 
    left_join(trt) %>% 
    filter(plot_mani==0) %>% 
    group_by(site_proj_comm, calendar_year, permutation) %>% 
    summarize_at(vars(mpd, 
                      mntd), 
                 list(mean=mean), na.rm=T) %>% #average across plots
    ungroup() %>% 
    rename(MNTD_phylo_ctl=mntd_mean,
           MPD_phylo_ctl=mpd_mean)
  
  RRmpdMNTDSubset <- pd.all %>% 
    left_join(trt) %>% 
    filter(plot_mani!=0) %>% 
    left_join(mpdMNTDSubsetPermCtl) %>% 
    mutate(RR_MNTD_phylo=log(mntd/MNTD_phylo_ctl),
           RR_MPD_phylo=log(mpd/MPD_phylo_ctl))
  
  RRmpdMNTDSubsetRaw <- RRmpdMNTDSubset %>% 
    filter(permutation==0) %>% 
    select(site_proj_comm, calendar_year, plot_id, RR_MNTD_phylo, RR_MPD_phylo)
  
  mpdMNTDSubsetSESrr <- RRmpdMNTDSubset %>% 
    filter(permutation>0) %>% 
    group_by(site_proj_comm, calendar_year, plot_id) %>% 
    summarize_at(vars(RR_MNTD_phylo, 
                      RR_MPD_phylo), 
                 list(mean=mean, sd=sd), na.rm=T) %>% #average across permutation
    ungroup() %>% 
    left_join(RRmpdMNTDSubsetRaw) %>% 
    mutate(SES_RR_MNTD_phylo=(RR_MNTD_phylo-RR_MNTD_phylo_mean)/RR_MNTD_phylo_sd,
           SES_RR_MPD_phylo=(RR_MPD_phylo-RR_MPD_phylo_mean)/RR_MPD_phylo_sd) %>% 
    select(site_proj_comm, calendar_year, plot_id, SES_RR_MNTD_phylo, SES_RR_MPD_phylo)
  
  allSubset <- mpdMNTDSubsetRRses %>% 
    full_join(mpdMNTDSubsetSESrr) %>% 
    left_join(trt)
  
  #bind values into RR dataframe
  phylogeneticDiversityMetrics <- rbind(phylogeneticDiversityMetrics, allSubset)

}

#save output:
write.table(phylogeneticDiversityMetrics, paste(my.wd, "CoRRE_PD_metrics_non_weighted_April2023.csv", sep=","))