
#####################################################
### Working to get functional dispersion estimate ###
#####################################################

library(FD)
library(tidyverse)

setwd("~/Dropbox/sDiv_sCoRRE_shared/")

# Working with just pplot from KNZ to start

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

# Next to do: 
# (1) Check how best to deal with categorical traits 
# (2) Merge back with plot treatments
# (3) Do this for the whole dataset 
# (4) Plot to see what is going on

test_FDis<-as.data.frame(test$FDis)

test_FDis_w_trt<-cbind(environmental_knz, test_FDis) # make sure these line up - plot ID's not available in datasets 

#drop some treatments 
test_FDis_w_trt_1<-as.data.frame(test_FDis_w_trt[!(test_FDis_w_trt$treatment %in% c("NCP0", "NCP1", "NCP2", "NCP3", "NNP0", "NNP1", "NNP2", "NNP3")),])
names(test_FDis_w_trt_1)[12] <- "FDis"

ggplot(test_FDis_w_trt_1, aes(x = treatment, y = FDis)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_wrap(~calendar_year)
