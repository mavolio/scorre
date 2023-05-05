################################################################################
##  sCoRRE_functionalDiversityMetrics.R: Calculating functional diversity metrics.
##
##  Authors: Kimberly Komatsu, Magda Garbowski, Kevin Wilcox, Josep Padulles Cubino
##  Date created: April 7, 2021
################################################################################

library(FD)
library(car)
library(gawdis)
library(tidyverse)

##### set working directory #####
setwd("~/Dropbox/sDiv_sCoRRE_shared/")
setwd("C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\") # Kim's laptop/desktop
setwd("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/") #Padu's wd
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Working groups\\sDiv\\") # Kevin's laptop wd


##### defining functions #####
## Standard Error Function:
se <- function(x, na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}


##### data import and cleaning #####

# trait data
traits <- read.csv('CoRRE data\\trait data\\AllTraits\\CoRRE_allTraitData_April2023.csv') %>% 
  select(species_matched, leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation) %>% 
  filter(growth_form!="moss", species_matched!="") %>% #keeping lycophytes
  mutate(mycorrhizal=ifelse(mycorrhizal_type=="none", 'no', ifelse(mycorrhizal_type=="uncertain", "unk", "yes"))) %>% 
  select(-mycorrhizal_type) %>% 
  rename(mycorrhizal_type=mycorrhizal) %>% 
  mutate(photo_path=ifelse(photosynthetic_pathway=="possible C4"|photosynthetic_pathway=="possible C4/CAM", "C4", ifelse(photosynthetic_pathway=="possible CAM", "CAM",photosynthetic_pathway))) %>% 
  select(-photosynthetic_pathway) %>%
  rename(photosynthetic_pathway=photo_path) %>%
  drop_na() #only keep trait data that is complete for all traits (drops 600 species with only categorical trait data)


##### testing normality #####
hist(traits$leaf_C.N)
qqPlot(traits$leaf_C.N)
shapiro.test(traits$leaf_C.N)

hist(traits$LDMC)
qqPlot(traits$LDMC)
shapiro.test(traits$LDMC)

hist(traits$SLA)
qqPlot(traits$SLA)
shapiro.test(traits$SLA)

hist(traits$plant_height_vegetative)
qqPlot(traits$plant_height_vegetative)
shapiro.test(traits$plant_height_vegetative)

hist(traits$rooting_depth)
qqPlot(traits$rooting_depth)
shapiro.test(traits$rooting_depth)

hist(traits$seed_dry_mass)
qqPlot(traits$seed_dry_mass)
shapiro.test(traits$seed_dry_mass)


##### log transform and scale continuous traits #####
traitsScaled <- traits %>%
  mutate_at(vars(leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass), log) %>% 
  mutate_at(vars(leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass), scale) #scale continuous traits

colnames(traitsScaled) <- c('species_matched', 'leaf_C.N', 'LDMC', 'SLA', 'plant_height_vegetative', 'rooting_depth', 'seed_dry_mass', 'growth_form', 'lifespan', 'clonal', 'n_fixation', 'mycorrhizal_type', 'photosynthetic_pathway')

#making categorical traits factors 
traitsScaled[,c(8:13)] <- lapply(traitsScaled[,c(8:13)], as.factor)

#testing normality
hist(traitsScaled$leaf_C.N)
qqPlot(traitsScaled$leaf_C.N)
shapiro.test(traitsScaled$leaf_C.N)
#log W = 0.9714, p-value < 2.2e-16
#sqrt W = 0.89825, p-value < 2.2e-16

hist(traitsScaled$LDMC)
qqPlot(traitsScaled$LDMC)
shapiro.test(traitsScaled$LDMC)
#log W = 0.92104, p-value < 2.2e-16
#sqrt W = 0.97428, p-value < 2.2e-16

hist(traitsScaled$SLA)
qqPlot(traitsScaled$SLA)
shapiro.test(traitsScaled$SLA)
#log W = 0.96681, p-value < 2.2e-16
#sqrt W = 0.92173, p-value < 2.2e-16

hist(traitsScaled$plant_height_vegetative)
qqPlot(traitsScaled$plant_height_vegetative)
shapiro.test(traitsScaled$plant_height_vegetative)
#log W = 0.99327, p-value = 3.754e-07
#sqrt W = 0.82253, p-value < 2.2e-16

hist(traitsScaled$rooting_depth)
qqPlot(traitsScaled$rooting_depth)
shapiro.test(traitsScaled$rooting_depth)
#log W = 0.99503, p-value = 1.445e-05
#sqrt W = 0.8344, p-value < 2.2e-16

hist(traitsScaled$seed_dry_mass)
qqPlot(traitsScaled$seed_dry_mass)
shapiro.test(traitsScaled$seed_dry_mass)
#log W = 0.99679, p-value = 0.001055
#sqrt W = 0.7173, p-value < 2.2e-16


##### relative cover datasets #####

# species relative cover data
relCoverRaw <- read.csv("CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeCover_Jan2023.csv") %>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>%
  mutate(plot_id=ifelse(project_name=='NSFC', paste(plot_id, treatment, sep='__'), plot_id)) %>% 
  select(site_code:community_type, site_proj_comm, calendar_year:relcov)

# corre to try species names key
corre_to_try <- read.csv("CoRRE data\\trait data\\corre2trykey_2021.csv") %>%
  select(genus_species, species_matched) %>%
  unique()

# merge species names and remove all mosses -- moss key to remove mosses from species comp data
moss_sp_vec <- read.csv("CoRRE data\\trait data\\sCoRRE categorical trait data_12142022.csv") %>%
  select(species_matched, leaf_type) %>%
  mutate(moss = ifelse(leaf_type=="moss", "moss","non-moss")) %>%
  filter(moss=="moss") %>%
  pull(species_matched)

relCovClean <- relCoverRaw %>%
  left_join(corre_to_try, by="genus_species") %>%
  filter(!species_matched  %in% moss_sp_vec) %>%
  mutate(plot_id=ifelse(site_proj_comm=='DL_NSFC_0', paste(plot_id, treatment, sep='__'), plot_id))

rm(moss_sp_vec)

##### treatment data #####
trt <- read.csv('CoRRE data\\CoRRE data\\community composition\\CoRRE_RawAbundance_Jan2023.csv') %>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id) %>%
  unique() %>%
  left_join(read.csv('CoRRE data\\CoRRE data\\basic dataset info\\ExperimentInfo.csv')) %>%
  group_by(site_code, project_name, community_type) %>%
  mutate(experiment_length=max(treatment_year)) %>%
  ungroup() %>%
  mutate(plot_id=ifelse(project_name=='NSFC', paste(plot_id, treatment, sep='__'), plot_id)) %>% 
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, trt_type, experiment_length, plot_mani, n, p, CO2, precip, temp) %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='_'))


##### calculate functional dispersion - loop through sites #####
functionalDiversityMetrics <- {}
site_vector <- unique(relCovClean$site_code) # do this for site_code only to reshuffle accurately

for(s in 1:length(site_vector)){
  
  #relative cover data from each site
  relCoverSubset <- relCovClean %>%
    filter(site_code==site_vector[s]) #%>% 
    # mutate(relcov2=ifelse(relcov>0, 1, 0)) %>% 
    # select(-relcov) %>% 
    # rename(relcov=relcov2)
  
  #species vector for pulling traits from relative cover
  sppSubset <- data.frame(genus_species = unique(relCoverSubset$genus_species), dummy=1) %>%
    left_join(corre_to_try, by="genus_species") %>%
    unique() 
  
  sppSubsetVector <- sppSubset %>%
    na.omit() %>% 
    pull(species_matched) %>% 
    unique()
  
  #subset trait data to just include species in the relative cover data
  traitsSubset <- traitsScaled %>%
    filter(species_matched %in% sppSubsetVector)
  
  #dataframe with species present in both the trait database and the relative cover data base
  speciesSubsetKeep <- data.frame(species_matched = unique(traitsSubset$species_matched),
                                  dummy_traits=2) %>%
    arrange(species_matched)
  
  #vector of species not in trait database (but in relative abundance data) to remove from species abundance data
  speciesSubsetRemove <- sppSubset %>%
    full_join(speciesSubsetKeep, by="species_matched") %>%
    filter(is.na(dummy_traits)) %>%
    pull(genus_species)
  
  #abundance dataset with species removed that do not have trait information
  relCoverSubsetKeep <- relCoverSubset %>%
    filter(!genus_species %in% speciesSubsetRemove) #removing species without trait information
  
  #abundance data into wide format
  relCoverWideSubset <- relCoverSubsetKeep %>%
    select(-genus_species) %>%
    group_by(site_code, project_name, community_type, site_proj_comm, calendar_year, treatment_year, treatment,
             block, plot_id, data_type, version, species_matched) %>%
    summarize(relcov=sum(relcov, na.rm=T)) %>%
    ungroup() %>%
    spread(key=species_matched, value=relcov) %>%
    replace(is.na(.), 0)
  
  #plot information
  plotInfoSubset <- relCoverWideSubset %>%
    select(site_code:version)
  
  #cover data
  relCoverWideSubset2 <- relCoverWideSubset %>%
    select(-site_code:-version) %>% 
    mutate(identifier=paste(plotInfoSubset$site_proj_comm, plotInfoSubset$calendar_year, plotInfoSubset$plot_id, sep="::")) %>% 
    column_to_rownames("identifier") 
  
  #dbFD function requires species names in trait data frame be arranged A-Z and identical order to the abundance data 
  traitsSubsetArranged <- traitsSubset %>%
    arrange(species_matched) %>%
    column_to_rownames("species_matched") %>%
    # mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
    # mutate_at(.vars=c("growth_form", "photosynthetic_pathway","lifespan", "clonal", "mycorrhizal_type", "n_fixation"), funs(as.numeric(as.factor(.)))) %>%
    mutate_at(.vars=c("leaf_C.N", "LDMC", "SLA", "plant_height_vegetative", "rooting_depth", "seed_dry_mass"), funs(as.numeric(.))) %>% 
    select(growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation, 
           leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass)
  
  ### Calculate functional diversity metrics ###
  relCoverMatrixSubset <- as.matrix(relCoverWideSubset2)
  traitMatrixSubset <- as.matrix(gawdis(traitsSubsetArranged, w.type = "optimized", opti.maxiter = 200, groups.weight=T, groups = c(1,2,3,4,5,6,7,7,7,8,9,10))) # because of NAs, use optimized

  #FDis and RaoQ
  FDsubset <- dbFD(x=as.matrix(traitsSubsetArranged), # matrix of traits
                  a=as.matrix(relCoverWideSubset2), # matrix of species
                  w.abun=T, # weight by abundance
                  cor="cailliez", # use Cailliez correlations because Euclidean distances could be calculated
                  calc.FRic=F, calc.FDiv=F, calc.CWM=F)

  FDsubset2 <- do.call(cbind.data.frame, FDsubset) %>%
    rownames_to_column(var = "identifier") %>% 
    separate(identifier, into=c("site_proj_comm", "calendar_year","plot_id"), sep="::") %>%
    mutate(calendar_year = as.numeric(calendar_year))

  #MPD and MNTD
  mpdMNTDSubset <- data.frame(
    plotInfoSubset[,c("site_proj_comm", "calendar_year", "plot_id")],
    MNTD_traits = picante::mntd(relCoverMatrixSubset, traitMatrixSubset),
    MPD_traits = picante::mpd(relCoverMatrixSubset, traitMatrixSubset)) %>% 
    full_join(plotInfoSubset)
  
  # distanceSubset <- FD %>%
  #   full_join(mpdMNTD)
  
  #null distributions for MPD and MNTD
  ses <- {}
  sesVector <- c(1:2)
  for(n in 1:length(sesVector)){
    traitsSubsetSES <- traitsSubsetArranged %>%
      rownames_to_column(var = "identifier") %>% 
      mutate(spp_shuffling = sample(identifier, size = n(), replace = FALSE)) %>% 
      select(-identifier) %>% 
      column_to_rownames(var="spp_shuffling")
    
    traitMatrixSubsetSES <- as.matrix(gawdis(traitsSubsetSES, w.type = "optimized", opti.maxiter = 200, groups.weight=T, groups = c(1,2,3,4,5,6,7,7,7,8,9,10)))
    
    mpdMNTDSubsetSES <- data.frame(
      plotInfoSubset[,c("site_proj_comm", "calendar_year", "plot_id")],
      MNTD_traits = picante::mntd(relCoverMatrixSubset, traitMatrixSubsetSES),
      MPD_traits = picante::mpd(relCoverMatrixSubset, traitMatrixSubsetSES)) %>% 
      mutate(permutation=sesVector[n])
    
    ses <- rbind(ses, mpdMNTDSubsetSES) 
    
    rm(list=ls()[grep("SES", ls())])
  }
  
  #mean MPD and MNTD for null distribution to create SES MPD and MNTD
  sesSubsetMean <- ses %>% 
    group_by(site_proj_comm, calendar_year, plot_id) %>% 
    summarise(across(c('MNTD_traits', 'MPD_traits'), list(mean=mean, sd=sd))) %>% 
    ungroup()
  
  #calculate SES values for MPD and MNTD
  mpdMNTDSubsetSES <- mpdMNTDSubset %>% 
    full_join(sesSubsetMean) %>% 
    mutate(MNTD_traits_ses=(MNTD_traits-MNTD_traits_mean)/MNTD_traits_sd,
           MPD_traits_ses=(MPD_traits-MPD_traits_mean)/MPD_traits_sd) %>% 
    select(site_proj_comm, site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, MNTD_traits, MNTD_traits_ses, MPD_traits, MPD_traits_ses)
  
  #subset control plots
  mpdMNTDSubsetSESctl <- mpdMNTDSubsetSES %>% 
    left_join(trt) %>% 
    filter(plot_mani==0) %>% 
    group_by(site_proj_comm, calendar_year) %>% 
    summarise(across(c('MNTD_traits', 'MPD_traits', 'MNTD_traits_ses', 'MPD_traits_ses'), list(mean=mean))) %>% 
    ungroup()
  
  #calculate lnRR for ses values of MPD and MNTD
  mpdMNTDSubsetRRses <- mpdMNTDSubsetSES %>% 
    left_join(trt) %>% 
    full_join(mpdMNTDSubsetSESctl) %>% 
    mutate(RR_MNTD_traits=ifelse(plot_mani>0, log(MNTD_traits/MNTD_traits_mean), NA),
           RR_MPD_traits=ifelse(plot_mani>0, log(MPD_traits/MPD_traits_mean), NA),
           RR_MNTD_traits_ses=ifelse(plot_mani>0, ((MNTD_traits_ses-MNTD_traits_ses_mean)/MNTD_traits_ses_mean), NA), #percent difference for ses due to neg values
           RR_MPD_traits_ses=ifelse(plot_mani>0, ((MPD_traits_ses-MPD_traits_ses_mean)/MPD_traits_ses_mean), NA)) %>% #percent difference for ses due to neg values
    select(site_proj_comm, calendar_year, plot_id, MNTD_traits, MNTD_traits_ses, MPD_traits, 
           MPD_traits_ses, RR_MNTD_traits, RR_MPD_traits, RR_MNTD_traits_ses, RR_MPD_traits_ses)
  
  #lnRR MPD and MNTD for null distribution to create SES lnRR MPD and SES lnRR MNTD
  ses2 <- ses %>% 
    full_join(plotInfoSubset)
  
  #bind onto raw data
  mpdMNTDSubsetPerm <- mpdMNTDSubset %>% 
    mutate(permutation='0') %>% 
    rbind(ses2) 
  
  #average MPD and MNTD in control plots by permutation
  mpdMNTDSubsetPermCtl <- mpdMNTDSubsetPerm %>% 
    left_join(trt) %>% 
    filter(plot_mani==0) %>% 
    group_by(site_proj_comm, calendar_year, permutation) %>% 
    summarize_at(vars(MNTD_traits, 
                      MPD_traits), 
                 list(mean=mean), na.rm=T) %>% #average across plots
    ungroup() %>% 
    rename(MNTD_traits_ctl=MNTD_traits_mean,
           MPD_traits_ctl=MPD_traits_mean)
  
  RRmpdMNTDSubset <- mpdMNTDSubsetPerm %>% 
    left_join(trt) %>% 
    filter(plot_mani!=0) %>% 
    left_join(mpdMNTDSubsetPermCtl) %>% 
    mutate(RR_MNTD_traits=log(MNTD_traits/MNTD_traits_ctl),
           RR_MPD_traits=log(MPD_traits/MPD_traits_ctl))
  
  RRmpdMNTDSubsetRaw <- RRmpdMNTDSubset %>% 
    filter(permutation==0) %>% 
    select(site_proj_comm, calendar_year, plot_id, RR_MNTD_traits, RR_MPD_traits)
  
  mpdMNTDSubsetSESrr <- RRmpdMNTDSubset %>% 
    filter(permutation>0) %>% 
    group_by(site_proj_comm, calendar_year, plot_id) %>% 
    summarize_at(vars(RR_MNTD_traits, 
                      RR_MPD_traits), 
                 list(mean=mean, sd=sd), na.rm=T) %>% #average across permutation
    ungroup() %>% 
    left_join(RRmpdMNTDSubsetRaw) %>% 
    mutate(SES_RR_MNTD_traits=(RR_MNTD_traits-RR_MNTD_traits_mean)/RR_MNTD_traits_sd,
           SES_RR_MPD_traits=(RR_MPD_traits-RR_MPD_traits_mean)/RR_MPD_traits_sd) %>% 
    select(site_proj_comm, calendar_year, plot_id, SES_RR_MNTD_traits, SES_RR_MPD_traits)
  
  allSubset <- FDsubset2 %>% 
    full_join(mpdMNTDSubsetRRses) %>% 
    full_join(mpdMNTDSubsetSESrr) %>% 
    left_join(trt)
  
  #bind values into RR dataframe
  functionalDiversityMetrics <- rbind(functionalDiversityMetrics, allSubset)
  
  rm(list=ls()[grep("Subset", ls())])
  rm(list=ls()[grep("subset", ls())])
  rm(list=ls()[grep("ses", ls())])
}


# write.csv(functionalDiversityMetrics, 'paper 2_PD and FD responses\\data\\CoRRE_functionalDiversity_2023-04-25.csv',row.names=F)