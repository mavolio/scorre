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
  #continuous

continuous_for_corre$TRY_trait <- revalue(continuous_for_corre$trait_name, c(
  "seed mass" = "seed_dry_mass",
  "maximum whole plant height" = "plant_height_generative",
  "leaf carbon content per leaf nitrogen content" = "leaf_C:N",
  "leaf photosynthetic rate per leaf dry mass" = "40",
  "leaf phosphorus content per leaf dry mass" = "leaf_P",
  "leaf area" = "leaf_area",
  "leaf carbon content per leaf dry mass" = "leaf_C",
  "leaf life span" = "leaf_longevity",
  "leaf stomatal conductance for H2O per leaf area" = "stomatal_conductance",
  "leaf nitrogen content per leaf area" = "50",
  "leaf photosynthetic rate per leaf area" = "photosynthesis rate",
  "leaf nitrogen content per leaf dry mass" = "leaf_N",
  "leaf phosphorus content per leaf area" = "51",
  "leaf area per leaf dry mass" = "SLA",
  "whole plant height" = "plant_height_vegetative",
  "leaf dry mass per leaf fresh mass" = "LDMC",
  "leaf dry mass" = "leaf_dry_mass",
  "leaf thickness" = "leaf_thickness"
  #,
  #"seed length" = "seed_length"
))


#re-attach to CoRRE species names
#merge each dataframe with corre_species



final_continuous <- merge(continuous_for_corre, corre_species, by.x="scrubbed_species_binomial",by.y="Name_matched",all.x = TRUE)
#trim unnecessary columns to simplify the spreadsheets before uploading to Dropbox
#final_continuous <- final_continuous[,c("cleaned_trait_value","method",#"project_pi",#"id","unit",
                                 #       "TRY_trait","genus species")]
#final_continuous <- unique(final_continuous)


#write.csv(final_continuous,"BIEN_continuous_traits_2-20-20.csv")


#Merge continuous traits with tree dataframe and remove trees
species_families_trees_compelete <- read.csv("C:/Users/ohler/Downloads/species_families_trees_compelete.csv")


with_tree <- merge(final_continuous, species_families_trees_compelete, by.x="scrubbed_species_binomial",by.y = "species_matched", all.x=TRUE)

without_tree <- subset(with_tree, tree.non.tree != "tree")



#make a genus column
without_tree <- separate(without_tree, col = "scrubbed_species_binomial", into = c("genus","species"), " ",remove=FALSE)
without_tree <- without_tree[,c("scrubbed_species_binomial","genus","TRY_trait","cleaned_trait_value","id","family")]



#make wide form
without_tree <- unique(without_tree)
wide <- spread(without_tree, key="TRY_trait",value="cleaned_trait_value")
wide <- subset(wide, select = -c(id) )



#remove duplicate rows
wide <- unique(wide)




#write.csv(wide, "BIEN_for_scorre.csv")



######################################
#Categorical
categorical_for_corre <- subset(bien_data, 
                                trait_name == "whole plant growth form"|
                                  trait_name =="whole plant vegetative phenology"|
                                  trait_name ==  "flower pollination syndrome"|
                                  trait_name ==  "whole plant sexual system"|
                                  trait_name ==   "leaf compoundness"|
                                  trait_name ==  "whole plant dispersal syndrome")

#categorical

categorical_for_corre$TRY_trait <- revalue(categorical_for_corre$trait_name, c(
  "whole plant growth form" = "Plant life form (Raunkiaer life form)",
  "whole plant vegetative phenology" = "Plant vegetative phenology (leaf phenology)",
  "flower pollination syndrome" = "Pollination syndrome",
  "whole plant sexual system" = "Flower secual syndrome (dichogamy, cleistogamy, dioecious, monoecious)",
  "whole plant dispersal syndrome" = "Dispersal syndrome",
  "leaf compoundness" = "Leaf compoundness"
))


final_categorical <- merge(categorical_for_corre, corre_species, by.x="scrubbed_species_binomial",by.y="Name_matched",all.x = TRUE)
#trim unnecessary columns to simplify the spreadsheets before uploading to Dropbox
final_categorical <- final_categorical[,c("trait_value","method",#"id",
                                          "TRY_trait","genus species")]
final_categorical <- unique(final_categorical)


#write.csv(final_categorical,"BIEN_categorical_traits_2-20-20.csv")









