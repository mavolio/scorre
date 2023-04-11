################################################################################
##  sCoRRE_functionalDiversityMetrics.R: Calculating functional diversity metrics.
##
##  Authors: Kimberly Komatsu, Magda Garbowski, Kevin Wilcox, Josep Padulles Cubino
##  Date created: April 7, 2021
################################################################################

library(FD)
library(car)
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
traits <- read.csv('CoRRE data\\trait data\\AllTraits\\CoRRE_allTraitData_March2023.csv') %>% 
  select(family, species_matched, leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation) %>% 
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

colnames(traitsScaled) <- c('family', 'species_matched', 'leaf_C.N', 'LDMC', 'SLA', 'plant_height_vegetative', 'rooting_depth', 'seed_dry_mass', 'growth_form', 'lifespan', 'clonal', 'n_fixation', 'mycorrhizal_type', 'photosynthetic_pathway')

#making categorical traits factors 
traitsScaled[,c(9:14)] <- lapply(traitsScaled[,c(9:14)], as.factor)

#testing normality
hist(traitsScaled$leaf_C.N)
qqPlot(traitsScaled$leaf_C.N)
shapiro.test(traitsScaled$leaf_C.N)
#log W = 0.96953, p-value < 2.2e-16
#sqrt W = 0.89825, p-value < 2.2e-16

hist(traitsScaled$LDMC)
qqPlot(traitsScaled$LDMC)
shapiro.test(traitsScaled$LDMC)
#log W = 0.94491, p-value < 2.2e-16
#sqrt W = 0.97428, p-value < 2.2e-16

hist(traitsScaled$SLA)
qqPlot(traitsScaled$SLA)
shapiro.test(traitsScaled$SLA)
#log W = 0.98562, p-value = 1.732e-12
#sqrt W = 0.92173, p-value < 2.2e-16

hist(traitsScaled$plant_height_vegetative)
qqPlot(traitsScaled$plant_height_vegetative)
shapiro.test(traitsScaled$plant_height_vegetative)
#log W = 0.98562, p-value = 1.732e-12
#sqrt W = 0.82253, p-value < 2.2e-16

hist(traitsScaled$rooting_depth)
qqPlot(traitsScaled$rooting_depth)
shapiro.test(traitsScaled$rooting_depth)
#log W = 0.98672, p-value = 7.119e-12
#sqrt W = 0.8344, p-value < 2.2e-16

hist(traitsScaled$seed_dry_mass)
qqPlot(traitsScaled$seed_dry_mass)
shapiro.test(traitsScaled$seed_dry_mass)
#log W = 0.99474, p-value = 5.24e-06
#sqrt W = 0.7173, p-value < 2.2e-16


##### relative cover datasets #####

# species relative cover data
relCoverRaw <- read.csv("CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeCover_Jan2023.csv") %>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="_")) %>%
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
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, trt_type, experiment_length, plot_mani, n, p, CO2, precip, temp) %>% 
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='_'))


##### calculate functional dispersion - loop through sites #####
distance <- {}
site_vector <- unique(relCovClean$site_code) # do this for site_code only to reshuffle accurately

for(s in 1:length(site_vector)){
  
  #relative cover data from each site
  relCoverSubset <- relCovClean %>%
    filter(site_code==site_vector[s]) %>% 
    mutate(relcov2=ifelse(relcov>0, 1, 0)) %>% 
    select(-relcov) %>% 
    rename(relcov=relcov2)
  
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
    select(-family) %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
    mutate_at(.vars=c("growth_form", "photosynthetic_pathway","lifespan", "clonal", "mycorrhizal_type", "n_fixation"), funs(as.numeric(as.factor(.)))) %>%
    mutate_at(.vars=c("leaf_C.N", "LDMC", "SLA", "plant_height_vegetative", "rooting_depth", "seed_dry_mass"), funs(as.numeric(.))) %>% 
    select(growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation, 
           leaf_C.N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass)
  
  ### Calculate functional diversity metrics ###
  relCoverMatrixSubset <- as.matrix(relCoverWideSubset2)
  traitMatrixSubset <- as.matrix(gowdis(traitsSubsetArranged))

  #FDis and RaoQ
  FDsubset <- dbFD(x=traitsSubsetArranged, # matrix of traits
                  a=relCoverWideSubset2, # matrix of species
                  w.abun=F, # don't weight by abundance
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
    
    traitMatrixSubsetSES <- as.matrix(gowdis(traitsSubsetSES))
    
    mpdMNTDSubsetSES <- data.frame(
      plotInfoSubset[,c("site_proj_comm", "calendar_year", "plot_id")],
      MNTD_traits = picante::mntd(relCoverMatrixSubset, traitMatrixSubsetSES),
      MPD_traits = picante::mpd(relCoverMatrixSubset, traitMatrixSubsetSES)) %>% 
      mutate(permutation=sesVector[n])
    
    ses <- rbind(ses, mpdMNTDSubsetSES) 
    
    rm(list=ls()[grep("SES", ls())])
  }
  
  #mean MPD and MNTD for null distribution to create SES MPD and MNTD
  ses2 <- ses %>% 
    full_join(plotInfoSubset)
  
  mpdMNTDSubsetSES <- mpdMNTDSubset %>% 
    mutate(permutation='0') %>% 
    rbind(ses2) 
  
#####START HERE: split into a dataframe that does the raw and SES of the values for each plot and then another that calculates lnRR and SES lnRR  
  
  #average MPD and MNTD in control plots by permutation
  sesCtl <- ses %>% 
    left_join(trt) %>% 
    filter(plot_mani==0) %>% 
    group_by(site_proj_comm, calendar_year, permutation) %>% 
    summarize_at(vars(MNTD_traits_permuted, 
                      MPD_traits_permuted), 
                 list(mean=mean), na.rm=T) %>% #average across plots and years
    ungroup() 
  
  sesRR <- ses %>% 
    left_join(trt) %>% 
    filter(plot_mani!=0) %>% 
    left_join(sesCtl) %>% 
    mutate(MNTD_traits_RR_ses=log(MNTD_traits_permuted/MNTD_traits_permuted_mean),
           MPD_traits_RR_ses=log(MPD_traits_permuted/MPD_traits_permuted_mean))
  
  #calculate SES values for MPD and MNTD
  mpdMNTDSubsetSES2 <- mpdMNTDSubset %>% 
    full_join(ses2) %>% 
    mutate(MNTD_traits_ses=(MNTD_traits_raw-MNTD_traits_permuted_mean)/MNTD_traits_permuted_sd,
           MPD_traits_ses=(MPD_traits_raw-MPD_traits_permuted_mean)/MPD_traits_permuted_sd) %>% 
    select(site_proj_comm, site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, MNTD_traits_raw, MNTD_traits_ses, MPD_traits_raw, MPD_traits_ses) %>% 
    full_join(FDsubset2)
  
  #bind values into original dataframes
  distance <- rbind(distance, mpdMNTDSubsetSES)
  
  rm(list=ls()[grep("Subset", ls())])
}


# write.csv(distance, 'paper 2_PD and FD responses\\data\\CoRRE_functionalDiversity_2023-03-30.csv',row.names=F)