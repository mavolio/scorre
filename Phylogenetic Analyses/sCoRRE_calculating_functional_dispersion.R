### Calculating functional diversity and dispersion using categorical and continuous traits
### Last updated: Dec 13, 2021
### R version: 4.1.1

library(FD)
library(tidyverse)
# install.packages("mFD")
# library("mFD")

setwd("~/Dropbox/sDiv_sCoRRE_shared/")
my.wd<-"/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/" #Padu's address
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\CoRRE data\\") # Kevin's laptop wd

### Start with PPlots to get things going

### Species relative cover data
pplots_raw <- read.csv("CoRRE data\\community composition\\CoRRE_RelativeCover_Dec2021.csv") %>%
  filter(project_name=="pplots") %>%
  mutate(genus_species = replace(genus_species, genus_species == "baptisia brachteata", "baptisia bracteata"))

corre_to_try <- read.csv("CoRRE data\\trait data\\corre2trykey_2021.csv") %>%
  dplyr::select(genus_species, species_matched) %>%
  unique(.)

pplots_long <- pplots_raw %>%
  dplyr::left_join(corre_to_try, by="genus_species")

###
### Get species vector for pulling traits from pplots relative cover
###
pplots_sp_df <- data.frame(genus_species = unique(pplots_raw$genus_species), dummy=1) %>%
  left_join(corre_to_try, by="genus_species") %>%
  unique(.) 

pplots_sp_vec <- pplots_sp_df %>%
  pull(species_matched)

###
### Read in and clean trait data
###
all_traits_raw <- read.csv("CoRRE data\\trait data\\Final Cleaned Traits\\sCoRRE categorical trait data_final_20211209.csv") %>%
  dplyr::select(species_matched, growth_form, photosynthetic_pathway, lifespan,  clonal, mycorrhizal_type, n_fixation) %>%
  mutate(photosynthetic_pathway = replace(photosynthetic_pathway, grep("possible", photosynthetic_pathway), NA)) %>%
  mutate(clonal = replace(clonal, clonal=="uncertain", NA)) %>%
  mutate(mycorrhizal_type = replace(mycorrhizal_type, mycorrhizal_type=="uncertain", NA)) %>%
  mutate(lifespan = replace(lifespan, lifespan=="uncertain", NA)) %>%
  filter(lifespan != "moss")
  
  
#### NOTES
###
### Continuous traits to include: LDMC, SLA, Vegetative_height, seed dry mass, seed number, rooting density, rooting depth 
###
###
### Categorical traits to include: growth_form, life_span, mycorrhizal_type, n_fixation, clonal, photosynthetic_pathway 
###
#########

### Checking trait categories -- THEY LOOK GOOD NOW!
with(all_traits_raw, unique(growth_form))  
with(all_traits_raw, unique(photosynthetic_pathway))  
with(all_traits_raw, unique(lifespan))  
with(all_traits_raw, unique(clonal))  
with(all_traits_raw, unique(mycorrhizal_type))   
with(all_traits_raw, unique(n_fixation)) 

### Subset trait data to just include the species present in the pplots relative cover data
pplots_traits <- all_traits_raw %>%
  filter(species_matched %in% pplots_sp_vec)

### Get data frame with species in trait data base and in pplots abundance data base
species_in_trait_data <- data.frame(species_matched = unique(pplots_traits$species_matched),
                                    dummy_traits=2) ## there are less species in the unique trait dataset than in the species comp data because they're things like "unknown forb"

### Get vector of species not in trait database (but in pplots relative abundance data) to remove from species abundance data
pplots_sp_to_remove <- pplots_sp_df %>%
  full_join(species_in_trait_data, by="species_matched") %>%
  filter(is.na(species_matched)) %>%
  pull(genus_species)

### Abundance data set with species removed that do not have trait information
pplots_relcov_unkn_sp_rm <- pplots_long %>%
  filter(!genus_species %in% pplots_sp_to_remove) # removing species without trait information

### get abundance data wide format
pplots_wide <- pplots_relcov_unkn_sp_rm %>%
  dplyr::select(-genus_species) %>%
  spread(key=species_matched, value=relcov) %>%
  replace(is.na(.), 0)

pplots_plot_info <- pplots_wide %>%
  dplyr::select(site_code:version)

pplots_relcov <- pplots_wide %>%
  dplyr::select(-site_code:-version) 

row.names(pplots_relcov) <- paste(pplots_plot_info$calendar_year, pplots_plot_info$plot_id, sep="_")

### dbFD function requires species names in trait data frame be arranged A-Z and identical order to the abundance data 
pplots_traits <- pplots_traits %>%
  arrange(species_matched) %>%
  column_to_rownames("species_matched")

### Changing all traits to factors
pplots_traits$growth_form <- as.factor(pplots_traits$growth_form)
pplots_traits$photosynthetic_pathway <- as.factor(pplots_traits$photosynthetic_pathway)
pplots_traits$lifespan <- as.factor(pplots_traits$lifespan)
pplots_traits$clonal <- as.factor(pplots_traits$clonal)
pplots_traits$mycorrhizal_type <- as.factor(pplots_traits$mycorrhizal_type)
pplots_traits$n_fixation <- as.factor(pplots_traits$n_fixation)

### create distance matrix for incorporation into dbFD function
pplots_gowdis <- gowdis(pplots_traits)

### NOT WORKING, ARG
FD_full <- dbFD(x=pplots_gowdis, a=pplots_relcov)

FD_full <- dbFD(x=pplots_traits, a=pplots_relcov)

### Trouble getting species by species distance matrix in Euclidean form... here's the error message:
# Error in dbFD(x = pplots_traits, a = pplots_relcov) : 
#   Species x species distance matrix was still is not Euclidean after 'sqrt' correction. Use another correction method.
# In addition: Warning messages:
#   1: In is.euclid(x.dist) : Zero distance(s)
# 2: In is.euclid(x.dist) : Zero distance(s)
# 3: In is.euclid(x.dist2) : Zero distance(s)


################## End of current script ##########################





### check abundance data
with(pplots_abun, table(calendar_year, plot_id)) ### one plot id per year, looking good!


####################
####################
####################

#### Running example
test_traits <- data.frame(trt1=runif(4, 0, 10),
                          trt2=runif(4, 0, 1),
                          trt3=runif(4, 0, 1),
                          trt4=runif(4, 0, 50),
                          trt5=as.factor(c("a","a","b","a")))
rownames(test_traits) <- c("sp1","sp2","sp3","sp4")

test_abun <- data.frame(sp1=runif(5, 0, 1),
                        sp2=runif(5, 0, 1),
                        sp3=runif(5, 0, 1),
                        sp4=runif(5, 0, 1))
rownames(test_abun) <- c("plot1","plot2","plot3","plot4","plot5")

ex1 <- dbFD(test_traits, test_abun)
cbind(ex1$FDiv, ex1$FEve)

# mixed trait types, NA's
data(dummy)
print(dummy)
ex1 <- dbFD(dummy$trait, dummy$abun)
ex1
# add variable weights
# 'cailliez' correction is used because 'sqrt' does not work
w<-c(1, 5, 3, 2, 5, 2, 6, 1)
ex2 <- dbFD(dummy$trait, dummy$abun, w, corr="cailliez")



# Working with just pplot from KNZ to splot.window()# Working with just pplot from KNZ to start
pplots_knz<-read.csv("paper 2_PD and FD responses/data/pplots species abundance wide.csv")
environmental_knz<-read.csv("paper 2_PD and FD responses/data/pplots environmental data.csv")

# Code from above 
contTraits <- read.csv('Trait Data/TRY Data/Gap_Filled/TRY_new.csv')%>%
  rename(species_matched=Species)%>%
  select(-X.1, -X, -Family, -Genus, -ObservationID)%>%
  group_by(species_matched)%>%
  summarise_all(funs(mean))%>%
  ungroup()

contTraitsSubset <- contTraits%>%
  rename(ssd=X4, rooting_depth=X6, SLA=X11, leaf_C_mass=X13, leaf_N_mass=X14, leaf_P_mass=X15, stem_diameter=X21, seed_mass=X26, seed_length=X27, leaf_thickness=X46, LDMC=X47, leaf_dry_mass=X55, germination_rate=X95, leaf_length=X144, leaf_width=X145, leaf_CN=X146, stem_conduit_density=X169, stem_conduit_diameter=X281, seed_number=X138, SRL=X1080)%>%
  select(-X18, -X50, -X78, -X163, -X223, -X224, -X237, -X282, -X289, -X3112, -X3113, -X3114, -X3120)

traits <- read.csv('CoRRE data/CoRRE data/trait data/sCoRRE categorical trait data - traits_complete_pre spot check_03102021.csv')%>%
  full_join(contTraitsSubset) %>%
  drop_na()%>%
  filter(leaf_P_mass<20, stem_diameter<0.5, seed_mass<50, seed_number<10000, leaf_width<40, stem_conduit_density<1000, stem_conduit_diameter<200)

traitsOutliersRemoved <- traits %>%
  filter(!leaf_type %in% c("microphyll","frond")) %>%
  filter(!species_matched %in% c("Centrolepis aristata", "Centrolepis strigosa", "Acorus calamus"))

trait_subset<-traitsOutliersRemoved [, c("species_matched","seed_mass", "seed_number","lifespan", "clonal",
                                         "rooting_depth", "LDMC", "SLA", "photosynthetic_pathway",
                                         "mycorrhizal_type","n_fixation")] # Missing veg height and root tissue density 

# Read in relative abundance data
sp_name_key <- read.csv("CoRRE data/CoRRE data/trait data/corre2trykey.csv")
rel_abun_df <- read.csv("CoRRE data/CoRRE data/community composition/CoRRE_RelativeAbundanceMar2021.csv") %>%
  left_join(dplyr::select(sp_name_key, genus_species, species_matched), by="genus_species") %>%
  drop_na(species_matched)

abund_species_vector <- unique(rel_abun_df$species_matched)
abund_species_vector<-as.data.frame(abund_species_vector)
abund_trait_merge<-merge(trait_subset, as.data.frame(abund_species_vector), by.x = "species_matched", by.y = "abund_species_vector")
cols<-c("lifespan", "clonal","photosynthetic_pathway", "mycorrhizal_type", "n_fixation")
abund_trait_merge[cols] <- lapply(abund_trait_merge[cols], as.factor)
abund_trait_merge$mycorrhizal_type[abund_trait_merge$mycorrhizal_type==""] <- NA
abund_trait_merge_cc<-abund_trait_merge[complete.cases(abund_trait_merge), ]

# Something funky is going on here - need to keep changing this back to a dataframe
sps_knz<-as.data.frame(colnames(pplots_knz))
sps_knz<-sps_knz[-1,]
sps_knz<-as.data.frame(sps_knz)
sps_knz<-gsub("\\.", " ", sps_knz$sps_knz)
sps_knz<-as.data.frame(sps_knz)

knz_trait_sps<-merge(sps_knz, abund_trait_merge_cc, by.x = "sps_knz", by.y = "species_matched", all.x = TRUE)
knz_trait_missing <- knz_trait_sps[rowSums(is.na(knz_trait_sps)) > 0,]
knz_trait_missing_sps_names<-knz_trait_missing$sps_knz
knz_trait_sps_complete<-knz_trait_sps[complete.cases(knz_trait_sps),]

rownames(knz_trait_sps_complete)<-knz_trait_sps_complete$sps_knz
knz_trait_sps_complete<-knz_trait_sps_complete[,-1]

# change factors to numeric - is this sketchy for the dbFD function? 
cols<-c("lifespan", "clonal","photosynthetic_pathway", "mycorrhizal_type", "n_fixation")
knz_trait_sps_complete[cols] <- lapply(knz_trait_sps_complete[cols], as.numeric)

row.names(pplots_knz)<-pplots_knz$X
pplots_knz<-pplots_knz[,-1]
colnames(pplots_knz)<-gsub("\\.", " ", colnames(pplots_knz))

pplots_knz<-pplots_knz[ , !names(pplots_knz) %in% knz_trait_missing_sps_names ]
pplots_knz[is.na(pplots_knz)] <- 0

test <- dbFD(knz_trait_sps_complete, pplots_knz) # I think this worked?

test_FDis<-as.data.frame(test$FDis)

test_FDis_w_trt<-cbind(environmental_knz, test_FDis) # make sure these line up - plot ID's not available in datasets 

# drop some treatments 
test_FDis_w_trt_1<-as.data.frame(test_FDis_w_trt[!(test_FDis_w_trt$treatment %in% c("NCP0", "NCP1", "NCP2", "NCP3", "NNP0", "NNP1", "NNP2", "NNP3")),])
names(test_FDis_w_trt_1)[12] <- "FDis"

ggplot(test_FDis_w_trt_1, aes(x = treatment, y = FDis)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_wrap(~calendar_year)

# To do next: 
# (1) Check how best to deal with categorical traits 
# (2) Do this for the whole dataset 
