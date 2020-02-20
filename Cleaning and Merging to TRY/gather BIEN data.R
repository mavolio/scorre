library(BIEN)
library(tidyverse)
library(plyr)


species_list<-BIEN_list_all()

corre_species <- read_csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/sCoRRE_tnrs_matching.csv")#species names are standardized to tnrs which is what BIEN uses
corre_species <- corre_species[,c("genus species","Name_matched")] #helps with merging later to prevent duplicate entries due to multiple subspecies
corre_species <- unique(corre_species)


sp.vector <- unique(corre_species$Name_matched)
bien_data <- BIEN_trait_species(species=sp.vector)
      #try grabbing all traits for desired species, THEN subset to desired traits
continuous_for_corre <- subset(bien_data, 
                         trait_name == "seed mass"|
                         trait_name ==  "maximum whole plant height"|
                         trait_name ==  "leaf carbon content per leaf nitrogen content"|
                         trait_name ==   "leaf photosynthetic rate per leaf dry mass"|
                         trait_name ==   "leaf phosphorus content per leaf dry mass"|
                         trait_name ==  "leaf area"|
                         trait_name ==   "leaf carbon content per leaf dry mass"|
                         trait_name ==  "leaf life span"|
                         trait_name == "leaf stomatal conductance for H2O per leaf area"|
                         trait_name ==  "leaf nitrogen content per leaf area"|
                         trait_name ==  "leaf photosynthetic rate per leaf area"|
                         trait_name ==  "leaf nitrogen content per leaf dry mass"|
                         trait_name ==  "leaf phosphorus content per leaf area"|
                         trait_name ==  "leaf area per leaf dry mass"|
                         trait_name ==  "whole plant height"|
                         trait_name ==  "leaf dry mass per leaf fresh mass"|
                         trait_name ==  "leaf dry mass"|
                         trait_name ==   "leaf carbon content per leaf dry mass"|
                         trait_name ==   "leaf thickness"|
                         trait_name ==  "root dry mass"|
                         trait_name ==  "seed length")
continuous_for_corre$trait_value <- as.numeric(continuous_for_corre$trait_value)


categorical_for_corre <- subset(bien_data, 
                                 trait_name == "whole plant growth form"|
                                 trait_name =="whole plant vegetative phenology"|
                                 trait_name ==  "flower pollination syndrome"|
                                 trait_name ==  "whole plant sexual system"|
                                  trait_name ==   "leaf compoundness"|
                                 trait_name ==  "whole plant dispersal syndrome")


#standardize units to fit TRY
    #convert LDMC (BIEN mg/g    TRY g/g)
     #leaf N per area (BIEN kg/m2  g/m2)
    #leaf C per area (BIEN kg/m2  g/m2)
    #leaf P per area (BIEN kg/m2  g/m2)
continuous_for_corre$cleaned_trait_value <- ifelse(
  continuous_for_corre$trait_name == "leaf dry mass per leaf fresh mass", continuous_for_corre$trait_value*1000, 
  ifelse(
    continuous_for_corre$trait_name == "leaf nitrogen content per area", continuous_for_corre$trait_value*0.001, 
    ifelse(continuous_for_corre$trait_name == "leaf carbon content per area", continuous_for_corre$trait_value*0.001,
           ifelse(continuous_for_corre$trait_name == "leaf phosphorous content per area", continuous_for_corre$trait_value*0.001,
      continuous_for_corre$trait_value))))
    

#Change BIEN trait names to fit TRY trait names
  #categorical

categorical_for_corre$TRY_trait <- revalue(categorical_for_corre$trait_name, c(
        "whole plant growth form" = "Plant life form (Raunkiaer life form)",
        "whole plant vegetative phenology" = "Plant vegetative phenology (leaf phenology)",
        "flower pollination syndrome" = "Pollination syndrome",
        "whole plant sexual system" = "Flower secual syndrome (dichogamy, cleistogamy, dioecious, monoecious)",
        "whole plant dispersal syndrome" = "Dispersal syndrome",
        "leaf compoundness" = "Leaf compoundness"
        ))


  #continuous

continuous_for_corre$TRY_trait <- revalue(continuous_for_corre$trait_name, c(
  "seed mass" = "Seed dry mass",
  "maximum whole plant height" = "Plant height generative",
  "leaf carbon content per leaf nitrogen content" = "Leaf carbon/nitrogen (C/N) ratio",
  "leaf photosynthetic rate per leaf dry mass" = "Leaf photosynthesis rate per leaf dry mass",
  "leaf phosphorus content per leaf dry mass" = "Leaf phosphorus content per leaf dry mass",
  "leaf area" = "Leaf area (in case of compound leaves undefined if leaf or leaflet, undefined if pe", ##I know that the name is cut off but this is how it appears in the TRY_downloaded_traits_full_list.csv
  "leaf carbon content per leaf dry mass" = "Leaf carbon (C) content per leaf dry mass",
  "leaf life span" = "Leaf lifespan (longevity)",
  "leaf stomatal conductance for H2O per leaf area" = "Stomata conductance per leaf area",
  "leaf nitrogen content per leaf area" = "Leaf nitrogen (N) content per leaf area",
  "leaf photosynthetic rate per leaf area" = "Leaf photosynthesis rate per leaf area",
  "leaf nitrogen content per leaf dry mass" = "Leaf nitrogen (N) content per leaf dry mass",
  "leaf phosphorus content per leaf area" = "Leaf phosphorus (P) content per leaf area",
  "leaf area per leaf dry mass" = "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiol",
  "whole plant height" = "Plant height vegetative",
  "leaf dry mass per leaf fresh mass" = "Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)",
  "leaf dry mass" = "Leaf dry mass (single leaf)",
  "leaf thickness" = "Leaf thickness",
  # "root dry mass" -> root_dry_mass ##no downloaded TRY trait that is only dry mass and not proportional to something like volume
  "seed length" = "Seed length"
))


#re-attach to CoRRE species names
#merge each dataframe with corre_species
final_categorical <- merge(categorical_for_corre, corre_species, by.x="scrubbed_species_binomial",by.y="Name_matched",all.x = TRUE)
#trim unnecessary columns to simplify the spreadsheets before uploading to Dropbox
final_categorical <- final_categorical[,c("trait_value","method",#"id",
        "TRY_trait","genus species")]
final_categorical <- unique(final_categorical)




final_continuous <- merge(continuous_for_corre, corre_species, by.x="scrubbed_species_binomial",by.y="Name_matched",all.x = TRUE)
#trim unnecessary columns to simplify the spreadsheets before uploading to Dropbox
final_continuous <- final_continuous[,c("cleaned_trait_value","method",#"project_pi",#"id","unit",
                                        "TRY_trait","genus species")]
final_continuous <- unique(final_continuous)

#write.csv(final_categorical,"BIEN_categorical_traits_2-20-20.csv")
#write.csv(final_continuous,"BIEN_continuous_traits_2-20-20.csv")



x <- data.frame(final_continuous$id[duplicated(final_continuous$id)])

y <- subset(bien_data, id == 13699987 | id ==13875237 | id==13717490)




BIEN_continuous_traits_1_18_20 <- read_csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/TRY/BIEN_continuous_traits_1-18-20.csv")

trial <- BIEN_continuous_traits_1_18_20[,c("trait_value","method","TRY_trait","genus species")]
trial.1 <- unique(trial)












