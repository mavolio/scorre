library(BIEN)
library(tidyverse)



#x <- BIEN_trait_traitbyspecies(trait = "whole plant height", species = "Carex capitata")



#x <- BIEN_trait_genus("Pectocarya heterocarpa", all.taxonomy = FALSE,
 #                     political.boundaries = FALSE, source.citation = F)



species_list<-BIEN_list_all()





#sev_traits <- read.csv("~/Misc Grad school/Traits/sev_all_traits_Nov2018_with_metadata_TEF.csv")
#sev_traits <- separate(sev_traits, taxon.new, c("genus","species"), sep = "_", remove = TRUE,convert = FALSE)

#sev_traits <- unite(sev_traits, taxon.forbien, c("genus","species"), sep=" ")

corre_species <- read_csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/sCoRRE_tnrs_matching.csv")#species names are standardized to tnrs which is what BIEN uses


sp.vector <- unique(corre_species$Name_matched)
#sp.vector <- species_list$species
trait.vector <- c(
  "whole plant growth form",
  "seed mass",
  "maximum whole plant height",
  "leaf carbon content per leaf nitrogen content",
  "whole plant vegetative phenology",
  "flower pollination syndrome",
  "whole plant sexual system",
  "whole plant dispersal syndrome",
  "leaf photosythetic rate per leaf dry mass",
  "leaf phosphorus content per leaf dry mass",
  "leaf area",
  "leaf carbon content per leaf dry mass",
  "leaf life span",
  "leaf stomatal conductance for H2O per leaf area",
  "leaf nitrogen content per leaf area",
  "leaf photosythetic rate per leaf area",
  "leaf nitrogen content per leaf dry mass",
  "leaf phosphorus content per leaf area",
  "leaf area per leaf dry mass",
  "whole plant height",
  "leaf dry mass per leaf fresh mass",
  "leaf dry mass",
  "leaf carbon content per leaf dry mass",
  "leaf thickness",
  "leaf compoundness",
  "root dry mass",
  "seed length")

bien_data <- BIEN_trait_species(species=sp.vector)
      #try grabbing all traits for desired species, THEN subset to desired traits
continuous_for_corre <- subset(bien_data, 
                         trait_name == "seed mass"|
                         trait_name ==  "maximum whole plant height"|
                         trait_name ==  "leaf carbon content per leaf nitrogen content"|
                         trait_name ==   "leaf photosythetic rate per leaf dry mass"|
                         trait_name ==   "leaf phosphorus content per leaf dry mass"|
                         trait_name ==  "leaf area"|
                         trait_name ==   "leaf carbon content per leaf dry mass"|
                         trait_name ==  "leaf life span"|
                         trait_name == "leaf stomatal conductance for H2O per leaf area"|
                         trait_name ==  "leaf nitrogen content per leaf area"|
                         trait_name ==  "leaf photosythetic rate per leaf area"|
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


#standardize methods to fit TRY
    #convert LDMC (BIEN mg/g    TRY g/g)
     #leaf N per area (BIEN kg/m2  g/m2)
    #leaf C per area (BIEN kg/m2  g/m2)
    #leaf P per area (BIEN kg/m2  g/m2)
    #figure out stomatal conductance
continuous_for_corre$cleaned_trait_value <- ifelse(
  continuous_for_corre$trait_name == "leaf dry mass per leaf fresh mass", continuous_for_corre$trait_value*1000, 
  ifelse(
    continuous_for_corre$trait_name == "leaf nitrogen content per area", continuous_for_corre$trait_value*0.001, 
    ifelse(continuous_for_corre$trait_name == "leaf carbon content per area", continuous_for_corre$trait_value*0.001,
           ifelse(continuous_for_corre$trait_name == "leaf phosphorous content per area", continuous_for_corre$trait_value*0.001,
      continuous_for_corre$trait_value))))
    
    
  

#Change BIEN trait names to fit TRY trait names
  #categorical
    #"whole plant growth form" -> lifeform (? don't know if this is the correct try name)
    #"whole plant vegetative phenology" -> (must find try name)
    #"flower pollination syndrome" -> Pollination syndrome (must find try name)
    #"whole plant sexual system" -> (mst find try name)
    #"whole plant dispersal syndrome" -> dispersal_mode
    #  "leaf compoundness" -> Leaf compoundness (must find try name)





  #continuous
    #"seed mass"     -> seed_dry_mass
 #"maximum whole plant height" -> plant_height_regenerative
  # "leaf carbon content per leaf nitrogen content" -> leaf_C:N
  #  "leaf photosythetic rate per leaf dry mass" -> Leaf photosynthesis rate per leaf dry mass (must find try name)
  #  "leaf phosphorus content per leaf dry mass" -> (must find try name)
  # "leaf area" -> leaf_area
  #  "leaf carbon content per leaf dry mass" -> Leaf carbon (C) content per leaf dry mass (must find try name)
  # "leaf life span" -> leaf_longevity
  #"leaf stomatal conductance for H2O per leaf area" -> stomata_conductance
  # "leaf nitrogen content per leaf area" -> Leaf nitrogen (N) content per leaf area (must find try name)
  # "leaf photosythetic rate per leaf area" -> Leaf photosynthesis rate per leaf area (must find try name)
  # "leaf nitrogen content per leaf dry mass" -> Leaf nitrogen (N) content per leaf dry mass (must find try name)
  # "leaf phosphorus content per leaf area" -> Leaf phosphorus (P) content per leaf area (must find try name)
  # "leaf area per leaf dry mass" -> SLA
  # "whole plant height" -> plant_height_regenerative
   # "leaf dry mass per leaf fresh mass" -> LDMC
   # "leaf dry mass" -> leaf_dry_mass
  #  "leaf carbon content per leaf dry mass" -> (must find try name)
  #  "leaf thickness" -> leaf_thickness
  # "root dry mass" -> root_dry_mass
  # "seed length" -> Seed length (must find try name)



#re-attach to CoRRE species names
#format so it can be easily combined with TRY data for CoRRE

