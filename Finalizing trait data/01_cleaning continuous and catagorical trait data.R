### Cleaning continuous trait data
### 
### Authors: Kimberly Komatsu (komatsuk@si.edu), Meghan Avolio (meghan.avolio@jhu.edu), Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Created: December 14, 2021; last updated: June 20, 2022
### Created with R version 4.1.1
###

###################################
### Set up working space
###################################

# rm(list=ls()) clean up workspace
#library(FD)
library(PerformanceAnalytics)
library(tidyverse)

# setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\CoRRE data\\") # Kevin's laptop wd
# setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data\\trait data\\Final TRY Traits\\') #Kim's 


############################
##### Continous traits #####
############################

### Continuous traits to include: LDMC, SLA, Vegetative_height, seed dry mass, seed number, rooting density, rooting depth 


##### Read in trait data #####
imputedRaw <- read.csv("Imputed Continuous_Traits\\Backtrans_GapFilled_sCorre.csv") %>%
  dplyr::select(-SRL) #need to remove SRL because input data was flawed


##### Look at histograms #####
# hist(imputedRaw$LDMC)
# hist(imputedRaw$SLA)
# hist(imputedRaw$plant_height_vegetative)
# hist(imputedRaw$seed_dry_mass)
# hist(imputedRaw$seed_number)
# hist(imputedRaw$rooting_depth)
# filter(imputedRaw, SLA == max(SLA))
# filter(imputedRaw, SLA > 200)


##### Removing Mosses #####
mossKey <- read.csv("sCoRRE categorical trait data_final_20211209.csv") %>%
  dplyr::select(species_matched, leaf_type) %>%
  mutate(moss = ifelse(leaf_type=="moss", "moss","non-moss")) %>%
  dplyr::select(-leaf_type)

imputedSubset <- imputedRaw %>%
  left_join(mossKey, by="species_matched") %>%
  mutate(moss=ifelse(moss=="moss","moss","non-moss")) %>%
  filter(moss!="moss") %>%
  dplyr::select(-moss) #%>%
# group_by(genus, family, species_matched) %>%
# summarize_at(vars(seed_dry_mass:seed_number), list(mean=mean, sd=sd), na.rm=T) %>%
# ungroup() #total of 1786 species remain

# #total species count = 2403 (so 617 species dropped [likely because they had zero trait data and therefore can't have any data imputed])
# sppList <- mossKey %>%
#   filter(moss=='non-moss')


# write.csv(imputedSubset, 'C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data\\trait data\\Final TRY Traits\\Imputed Continuous_Traits\\data to play with\\imputed_continuous_20220620.csv')


##### outlier check #####
hierarchy <- read.csv('Imputed Continuous_Traits\\hierarchy_info.txt')
# hierarchy2 <- filter(hierarchy, X0>0)

imputedMean <- read.delim('Imputed Continuous_Traits\\mean_gap_filled2.txt', header=T)
imputedSD <- read.delim('Imputed Continuous_Traits\\std_gap_filled2.txt', header=T)

#trying to link files, but won't bind becasue of column mismatch, checking with Franzi
# bhpmf_means <- imputedMean 
# Final_Mean_Traits <- cbind(hierarchy,bhpmf_means)


apply(imputedSD[], 2, quantile, probs=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

imputedSDLong <- imputedSD%>%
  pivot_longer(names_to='trait', values_to='sd', seed_dry_mass:seed_number)


ggplot(data=imputedSDLong, aes(x=trait, y=sd)) + geom_boxplot() + geom_point()


##### TO DO #####
# how were outliers removed before gap filling -- should maybe do it species by species, because there are some outliers for andro and achillea that might not be outliers in the entire dataset
# figure out how to merge the sd with the imputed backtransformed data
# determine how to do the sd cutoffs (we cannot just do quantiles because then the same amount of data is dropped from each trait, which is unfair to the good traits and generous to the bad traits)
# drop the outliers based on sd cutoffs and the extreme values
# correlations for each trait with imputed data (see code below)
# look up some values for species that we know and make sure they are right











##### Correlate imputed data with input data #####
originalTRY <- read.csv('Imputed Continuous_Traits\\Final_Input_sCorr.csv')%>%
  select(-X)%>%
  pivot_longer(names_to='trait', values_to='TRY_value', seed_dry_mass:seed_number)%>%
  na.omit()

traitsContCheck <- imputedSubset%>%
  select(-X)%>%
  pivot_longer(names_to='trait', values_to='imputed_value', seed_dry_mass:seed_number)%>%
  left_join(originalTRY)%>%
  na.omit()

ggplot(data=traitsContCheck, aes(x=TRY_value, y=imputed_value)) +
  geom_point() +
  geom_abline(slope=1) +
  facet_wrap(~trait, scales='free')
### traits that are very well correlated with input data
#seed_number

#getting figures for a few test species
ggplot(data=subset(traitsContCheck, species_matched=='Achillea millefolium'), aes(x=TRY_value, y=imputed_value)) +
  geom_point() +
  geom_abline(slope=1) +
  facet_wrap(~trait, scales='free')




# hist(imputedSubset$LDMC_mean)
# hist(imputedSubset$SLA_mean)
# hist(imputedSubset$plant_height_vegetative_mean)
# hist(imputedSubset$seed_dry_mass_mean)
# hist(imputedSubset$seed_number_mean)
# hist(imputedSubset$rooting_depth_mean)

rm(imputedRaw)

#########################################################################
### Categorical trait data 
#########################################################################

#### NOTES:
###
### Categorical traits to include: growth_form, life_span, mycorrhizal_type, n_fixation, clonal, photosynthetic_pathway 
###
#########


#all_traits_raw <- openxlsx::read.xlsx("sCoRRE categorical trait data_final_20211209.xlsx", sheet=1) %>%
traits_catag_clean <- read.csv("CoRRE data\\trait data\\Final Cleaned Traits\\sCoRRE categorical trait data_final_20211209.csv") %>%
  dplyr::select(species_matched, growth_form, photosynthetic_pathway, lifespan,  clonal, mycorrhizal_type, n_fixation) %>%
  mutate(photosynthetic_pathway = replace(photosynthetic_pathway, grep("possible", photosynthetic_pathway), NA)) %>%
  mutate(clonal = replace(clonal, clonal=="uncertain", NA)) %>%
  mutate(mycorrhizal_type = replace(mycorrhizal_type, mycorrhizal_type=="uncertain", NA)) %>%
  mutate(lifespan = replace(lifespan, lifespan=="uncertain", NA)) %>%
  filter(lifespan != "moss")


###########################################################
### Combine continuous and catagorical traits
###########################################################

traits_all <- dplyr::select(imputedSubset, -LDMC_sd:-rooting_depth_sd) %>%
  full_join(traits_catag_clean, by="species_matched")

# Find species that have continuous data but no categorical data -- 3 of these are mosses, so I will remove them from continuous data in that cleaning script
# TO DO: We will want to populate categorical data for three species: "Galium mollugo" "Heracleum sphondylium" "Trachypogon spicatus"