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

#### set directory ####
setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data\\')  #kim's computer
my.wd <- "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/" #padu


#### read data ####

#spp names
names <- read.csv('trait data\\corre2trykey_2021.csv') %>% 
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
comm2 <- read.table("CoRRE data\\community composition\\CoRRE_RelativeCover_Jan2023.csv", header=T, sep=",", fill = TRUE) %>% 
  right_join(trt_analysis) %>%
  mutate(trt_binary=ifelse(plot_mani>0, 1, 0)) %>% 
  mutate(drop=ifelse(site_code=="CDR" & treatment %in% c(2, 3, 4, 5, 7),
                     1, 0)) %>% #drop some of the CDR e001 and e002 treatments to prevent over-representation
  filter(drop==0) %>%
  filter(!(project_name %in% c('gap'))) #remove pulse light expt

#drop data points
dropExp <- comm2 %>% 
  left_join(names) %>% 
  filter(!is.na(species_matched)) %>% #remove species identified only to genus level
  select(-genus_species) %>% 
  rename(species=species_matched) %>%
  left_join(contTraitSpp) %>% 
  filter(!is.na(family)) %>% #remove species without continuous trait data
  group_by(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block, plot_id) %>% 
  summarize(totcov=sum(relcov), .groups='drop') #total cover of remaining species (7394 of plots are <0.8 totcov, will be dropped; 24.1% of plots)

hist(dropExp$totcov)

comm <- comm2 %>% 
  left_join(names) %>% 
  filter(!is.na(species_matched)) %>% #remove species identified only to genus level (drops 21,612 data points, 6.4% of data)
  select(-genus_species) %>% 
  rename(species=species_matched) %>% 
  left_join(contTraitSpp) %>% 
  filter(!is.na(family)) %>% #remove species without continuous trait data (drops 37,553 data points, 11.8% of data)
  left_join(dropExp) %>% 
  mutate(plot_id=ifelse(project_name=='NSFC', paste(plot_id, treatment, sep='__'), plot_id)) %>% #rename plot id for NSFC expt
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block, plot_id, species, relcov, totcov) %>% 
  mutate(plot_id2 = paste(site_code, project_name, community_type, calendar_year, plot_id, sep = "::")) %>%  #create new plot identifier
  group_by(plot_id2) %>% 
  mutate(richness=length(species)) %>% 
  ungroup() %>% 
  filter(richness>1, #remove plots with only one species (drops 688 data points, 0.2% of data)
         totcov>0.8) #remove plots with less than 80% cover of species with known trait values (drops 58,849 data points, 21.1% of data points)

#species list
spp <- comm %>% select(species) %>% unique() #1518 spp

#plot list
plots <- comm %>% select(site_code, project_name, community_type, plot_id) %>% unique() #3284 plots

#experiment list
expt <- comm %>% select(site_code, project_name, community_type) %>% unique() #121 experiments

#site list
sites <- unique(comm$site_code) #65 sites

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


#### calculate phylogenetic diversity metrics using one single tree (scenario 3) ####

#load tree
scorre.tree <- read.tree("Phylogenies\\scorre.phylo.tree.S3_20230427.tre")

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
  
  #calculate phylogenetic diversity metrics:
  mpd.raw <- as.data.frame(mpd(samp=spp, dis=distance,  abundance.weighted = T))
  mntd.raw <- as.data.frame(mntd(samp=spp, dis=distance,  abundance.weighted = T))
  
  pd.raw <- cbind(comm2[,1], mpd.raw, mntd.raw) %>% 
    rename(mpd="mpd(samp = spp, dis = distance, abundance.weighted = T)",
           mntd="mntd(samp = spp, dis = distance, abundance.weighted = T)") %>% 
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
    mntd.ses <- as.data.frame(mntd(samp=spp, dis=distanceShuffle,  abundance.weighted = T))
    
    ses <- cbind(comm2[,1], mpd.ses, mntd.ses) %>% 
      rename(mpd="mpd(samp = spp, dis = distanceShuffle, abundance.weighted = T)",
             mntd="mntd(samp = spp, dis = distanceShuffle, abundance.weighted = T)") %>% 
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
    filter(permutation==0) %>% 
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
write.csv(phylogeneticDiversityMetrics, "CoRRE_PD_metrics_weighted_April2023.csv", row.names=F)
