## Trait Coverage

# 1. Coverage by site
# 2. Coverage by plot
# 3. Coverage by manipulation
# 4. Coverage by experiment

## Libraries
library(here)
library(plyr)
library(tidyr)
library(ggplot2)

## Load Data
TRY_cat <- read.csv(here::here("Data", "TRY_trait_data_categorical.csv")) # categorical data
TRY_con <- read.csv(here::here("Data", "TRY_trait_data_continuous.csv")) # continuous data
Sp_lst <- read.csv(here::here("Data", "CoRRE_TRY_species_list.csv")) # need to use to combine trait data to CoRRE data
Rel_abun <- read.csv(here::here("Data", "CoRRE_relative_abundance_Nov2019.csv"), row.names = 1) # Relative abundance data for coverage


# Make sp names consistent between TRY and CoRRE
Rel_abun <- join(Rel_abun, Sp_lst[,c(1,2)], by = "genus_species", match = "first")
# Get rid of NA species
Rel_abun <- Rel_abun[-which(is.na(Rel_abun$species_matched)),] #~40,000 entries are NA
Rel_abun <- Rel_abun[,-1]
test <- aggregate(Rel_abun$relcov, by = list(site_code = Rel_abun$site_code, project_name = Rel_abun$project_name, community_type = Rel_abun$community_type,
                                             calendar_year = Rel_abun$calendar_year, treatment_year = Rel_abun$treatment_year, treatment = Rel_abun$treatment, 
                                             block = Rel_abun$block, plot_id = Rel_abun$plot_id, species_matched = Rel_abun$species_matched), 
                  FUN = sum)
names(test)[10] <- "relcov"
## Make dataframe with species as columns filled with relative abundances
abundcov <- spread(test, species_matched, relcov)
abundcov[is.na(abundcov)] <- 0

## Make dataframe with species as rows and traits as columns
TRY_cat$present <- 1
TRY_cat <- TRY_cat[,c(1,3,7)]
TRY_con$present <- 1
TRY_con <- TRY_con[,c(1,3,7)]
traits <- rbind(unique(TRY_cat), TRY_con) # dispersal mode is in TRY_cat x2
traits <- spread(traits, CleanTraitName, present)
traits <- traits[-which(is.na(traits$species_matched)),]
traits <- traits[order(traits$species_matched),]
rownames(traits) <- traits[,1]
traits <- traits[,-1]
traits[is.na(traits)] <- 0

###################
#### Coverage #####
###################
# Need to get rid of species with 0 trait values
abund_mat <- as.matrix(abundcov[,which(colnames(abundcov) %in% rownames(traits))])
trait_mat <- as.matrix(traits[which(colnames(abund_mat) %in% rownames(traits)),])
# Multiply abundance by trait matrix
mult <- abund_mat %*% trait_mat
# Add in site, project, ect
mult <- cbind(abundcov[,c(1:8)], mult)

long.mult <- gather(mult, trait, coverage, allelopathic:stomata_conductance)

ggplot(aes(x = site_code, y = coverage), data = long.mult) + 
  geom_boxplot() + 
  facet_grid(~trait)
