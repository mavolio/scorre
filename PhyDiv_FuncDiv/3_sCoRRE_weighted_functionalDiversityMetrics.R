################################################################################
##  sCoRRE_functionalDiversityMetrics.R: Calculating functional diversity metrics.
##
##  Authors: Kimberly Komatsu, Magda Garbowski, Kevin Wilcox, Josep Padulles Cubino
##  Date created: April 7, 2021
################################################################################

library(cluster)
library(fundiversity)
library(car)
library(tidyverse)


##### defining functions #####
## Standard Error Function:
se <- function(x, na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}


##### data import and cleaning #####

#spp names
names <- readRDS('PhyDiv_FuncDiv/sppList.rds')

#community data
comm <- readRDS('PhyDiv_FuncDiv/PD_FD_commData.rds')

#treatment data
trt <- readRDS('PhyDiv_FuncDiv/trt_info.rds')

# trait data
correGExTraitsContinuous <- read.csv('https://pasta.lternet.edu/package/data/eml/edi/1533/3/169fc12d10ac20b0e504f8d5ca0b8ee8')  %>% 
  filter(error_risk_overall<2|is.na(error_risk_overall)) %>% #drops 326 trait values
  select(family, species, trait, trait_value) 

correGExTraitsCategorical <- read.csv('https://pasta.lternet.edu/package/data/eml/edi/1533/3/5ebbc389897a6a65dd0865094a8d0ffd') %>% 
  select(-family, -source, -error_risk_overall)

traits <- rbind(correGExTraitsCategorical, correGExTraitsContinuous) %>% 
  pivot_wider(names_from=trait, values_from=trait_value) %>% 
  select(-leaf_area, -leaf_dry_mass, -leaf_type, -leaf_compoundness, -stem_support) %>% 
  mutate(across(c(growth_form, photosynthetic_pathway, lifespan, clonal, 
                  mycorrhizal_type, n_fixation_type), 
                as.factor),
         across(c(LDMC, SLA, SRL, leaf_N, plant_height_vegetative, seed_dry_mass), 
                as.numeric)) %>% 
  na.omit() #only keep trait data that is complete for all traits (drops 1157 species with only categorical trait data, 28.4% of species)


##### testing normality #####
hist(traits$leaf_N)
qqPlot(traits$leaf_N)
shapiro.test(traits$leaf_N)

hist(traits$LDMC)
qqPlot(traits$LDMC)
shapiro.test(traits$LDMC)

hist(traits$SLA)
qqPlot(traits$SLA)
shapiro.test(traits$SLA)

hist(traits$plant_height_vegetative)
qqPlot(traits$plant_height_vegetative)
shapiro.test(traits$plant_height_vegetative)

hist(traits$SRL)
qqPlot(traits$SRL)
shapiro.test(traits$SRL)

hist(traits$seed_dry_mass)
qqPlot(traits$seed_dry_mass)
shapiro.test(traits$seed_dry_mass)


##### log transform and scale continuous traits #####
traitsScaled <- traits %>%
  mutate_at(vars(LDMC, SLA, SRL, leaf_N, plant_height_vegetative, seed_dry_mass), log) %>% 
  mutate(across(c(LDMC, SLA, SRL, leaf_N, plant_height_vegetative, seed_dry_mass),
                ~ (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE)))) #scale continuous traits to 0-1


#testing normality
hist(traitsScaled$leaf_N)
qqPlot(traitsScaled$leaf_N)
shapiro.test(traitsScaled$leaf_N)
#log W = 0.99666, p-value = 4.562e-06

hist(traitsScaled$LDMC)
qqPlot(traitsScaled$LDMC)
shapiro.test(traitsScaled$LDMC)
#log W = 0.9404, p-value < 2.2e-16

hist(traitsScaled$SLA)
qqPlot(traitsScaled$SLA)
shapiro.test(traitsScaled$SLA)
#log W = 0.97353, p-value < 2.2e-16

hist(traitsScaled$plant_height_vegetative)
qqPlot(traitsScaled$plant_height_vegetative)
shapiro.test(traitsScaled$plant_height_vegetative)
#log W = 0.99595, p-value = 4.063e-07

hist(traitsScaled$SRL)
qqPlot(traitsScaled$SRL)
shapiro.test(traitsScaled$SRL)
#log W = 0.88705, p-value < 2.2e-16

hist(traitsScaled$seed_dry_mass)
qqPlot(traitsScaled$seed_dry_mass)
shapiro.test(traitsScaled$seed_dry_mass)
#log W = 0.99368, p-value = 5.985e-10


##### calculate functional dispersion - loop through sites #####
functionalDiversityMetrics <- {}
site_vector <- unique(comm$site_code) # do this for site_code only to reshuffle accurately

for(s in 1:length(site_vector)){
  
  print(s*100/length(site_vector))
  
  #relative cover data from each site
  relCoverSubset <- comm %>%
    filter(site_code==site_vector[s])
  
  #species vector for pulling traits from relative cover
  sppSubset <- data.frame(species = unique(relCoverSubset$species), dummy=1) 
  
  sppSubsetVector <- sppSubset %>%
    na.omit() %>% 
    pull(species) %>% 
    unique()
  
  #subset trait data to just include species in the relative cover data
  traitsSubset <- traitsScaled %>%
    filter(species %in% sppSubsetVector)
  
  #dataframe with species present in both the trait database and the relative cover data base
  speciesSubsetKeep <- data.frame(species = unique(traitsSubset$species),
                                  dummy_traits=2) %>%
    arrange(species)
  
  #vector of species not in trait database (but in relative abundance data) to remove from species abundance data
  speciesSubsetRemove <- sppSubset %>%
    full_join(speciesSubsetKeep, by="species") %>%
    filter(is.na(dummy_traits)) %>%
    pull(species)
  
  #abundance dataset with species removed that do not have trait information
  relCoverSubsetKeep <- relCoverSubset %>%
    filter(!species %in% speciesSubsetRemove) #removing species without trait information
  
  #abundance data into wide format
  relCoverWideSubset <- relCoverSubsetKeep %>%
    group_by(site_code, project_name, community_type, 
             calendar_year, treatment_year, treatment,
             block, plot_id, plot_id2, species) %>%
    summarize(relcov=sum(relcov, na.rm=T)) %>%
    ungroup() %>%
    spread(key=species, value=relcov) %>%
    replace(is.na(.), 0)
  
  #plot information
  plotInfoSubset <- relCoverWideSubset %>%
    select(site_code:plot_id2)
  
  #cover data
  relCoverWideSubset2 <- relCoverWideSubset %>%
    mutate(identifier=paste(plotInfoSubset$site_code, plotInfoSubset$project_name, plotInfoSubset$community_type, plotInfoSubset$calendar_year,                            plotInfoSubset$plot_id, sep="::")) %>% 
    select(-site_code:-plot_id2) %>% 
    column_to_rownames("identifier") 
  
  #dbFD function requires species names in trait data frame be arranged A-Z and identical order to the abundance data 
  traitsSubsetArranged <- traitsSubset %>%
    arrange(species) %>%
    column_to_rownames("species")
  
  ### Calculate functional diversity metrics ###
  relCoverMatrixSubset <- as.matrix(relCoverWideSubset2)
  traitMatrixSubset <- as.matrix(daisy(traitsSubsetArranged, metric = "gower"))
  
  # RaoQ
  rao_q <- fd_raoq(sp_com=relCoverMatrixSubset, dist_matrix=traitMatrixSubset) %>%
    rename(identifier=site) %>% 
    separate(identifier, into=c("site_code", "project_name", "community_type", "calendar_year","plot_id"), sep="::") %>%
    mutate(calendar_year = as.numeric(calendar_year),
           permutation=0)
  
  # # FDis
  # FDis <- fd_fdis() # can't use because only takes numeric traits
  # 
  # #below is using FD package, so slow! will take days to run 
  # FDsubset <- dbFD(x=traitsSubsetArranged, # matrix of traits
  #                  a=relCoverWideSubset2, # matrix of species
  #                  w.abun=T, # don't weight by abundance
  #                  cor="cailliez", # use Cailliez correlations because Euclidean distances could be calculated
  #                  calc.FRic=F, calc.FGR=F, calc.CWM=F, calc.FDiv=F)
  # 
  # FDsubset2 <- do.call(cbind.data.frame, FDsubset) %>%
  #   rownames_to_column(var = "identifier") %>% 
  #   separate(identifier, into=c("site_proj_comm", "calendar_year","plot_id"), sep="::") %>%
  #   mutate(calendar_year = as.numeric(calendar_year)) %>% 
  #   select(-FEve) %>% 
  #   mutate(permutation=0)
  
  #null distributions for RaoQ
  ses <- {}
  sesVector <- c(1:999)
  for(n in 1:length(sesVector)){
    
    traitsSubsetSES <- traitsSubsetArranged %>%
      rownames_to_column(var = "identifier") %>% 
      mutate(spp_shuffling = sample(identifier, size = n(), replace = FALSE)) %>% 
      select(-identifier) %>% 
      column_to_rownames(var="spp_shuffling")
    
    traitsSubsetArrangedSES <- traitsSubsetSES %>%
      rownames_to_column("species") %>% 
      arrange(species) %>%
      column_to_rownames("species")
    
    traitMatrixSubsetSES <- as.matrix(daisy(traitsSubsetArrangedSES, metric = "gower"))
    
    ses2 <- fd_raoq(sp_com=relCoverMatrixSubset, 
                         dist_matrix=traitMatrixSubsetSES) %>%
      rename(identifier=site) %>% 
      separate(identifier, into=c("site_code", "project_name", "community_type", "calendar_year","plot_id"), sep="::") %>%
      mutate(calendar_year = as.numeric(calendar_year),
             permutation=sesVector[n])
    
    # FDis - so slow, takes days
    # FDsubsetSES <- dbFD(x=traitMatrixSubsetSES, # matrix of traits
    #                     a=relCoverWideSubset2, # matrix of species
    #                     w.abun=T, # weight by abundance
    #                     cor="cailliez", 
    #                  calc.FRic=F, calc.FDiv=F, calc.CWM=F)
    # 
    # FDsubsetSES2 <- do.call(cbind.data.frame, FDsubsetSES) %>%
    #   rownames_to_column(var = "identifier") %>% 
    #   separate(identifier, into=c("site_proj_comm", "calendar_year","plot_id"), sep="::") %>%
    #   mutate(calendar_year = as.numeric(calendar_year)) %>% 
    #   select(-FEve) %>% 
    #     mutate(permutation=sesVector[n])
    
    ses <- rbind(ses, ses2) 
    
    rm(list=ls()[grep("SES", ls())])
  }
  
  #mean RaoQ for null distribution to create SES RaoQ
  sesSubsetMean <- ses %>% 
    group_by(site_code, project_name, community_type, calendar_year, plot_id) %>% 
    summarise(across(c('Q'), list(mean=mean, sd=sd))) %>% 
    ungroup()
  
  #calculate SES values for RaoQ 
  raoqSubsetSES <- rao_q %>% 
    full_join(sesSubsetMean) %>% 
    mutate(RaoQ_ses=(Q-Q_mean)/Q_sd) %>% 
    select(site_code, project_name, community_type, calendar_year, plot_id, RaoQ_ses)
  
  allSubset <- rao_q %>% 
    select(-permutation) %>% 
    full_join(raoqSubsetSES) %>% 
    left_join(trt)
  
  #bind values into dataframe
  functionalDiversityMetrics <- rbind(functionalDiversityMetrics, allSubset)
  
  rm(list=ls()[grep("Subset", ls())])
  rm(list=ls()[grep("subset", ls())])
  rm(list=ls()[grep("ses", ls())])
}

#save output:
saveRDS(functionalDiversityMetrics, file = "PhyDiv_FuncDiv/functionalDiversityMetrics.rds")
