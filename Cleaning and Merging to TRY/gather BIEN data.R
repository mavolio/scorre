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
final_categorical <- merge(categorical_for_corre, corre_species, by.x="scrubbed_species_binomial",by.y="Name_matched",all.x = TRUE)
#trim unnecessary columns to simplify the spreadsheets before uploading to Dropbox
final_categorical <- final_categorical[,c("trait_value","method",#"id",
        "TRY_trait","genus species")]
final_categorical <- unique(final_categorical)




final_continuous <- merge(continuous_for_corre, corre_species, by.x="scrubbed_species_binomial",by.y="Name_matched",all.x = TRUE)
#trim unnecessary columns to simplify the spreadsheets before uploading to Dropbox
#final_continuous <- final_continuous[,c("cleaned_trait_value","method",#"project_pi",#"id","unit",
                                 #       "TRY_trait","genus species")]
#final_continuous <- unique(final_continuous)

#write.csv(final_categorical,"BIEN_categorical_traits_2-20-20.csv")
#write.csv(final_continuous,"BIEN_continuous_traits_2-20-20.csv")



x <- data.frame(final_continuous$id[duplicated(final_continuous$id)])

y <- subset(bien_data, id == 13699987 | id ==13875237 | id==13717490)




BIEN_continuous_traits_1_18_20 <- read_csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/TRY/BIEN_continuous_traits_1-18-20.csv")

trial <- BIEN_continuous_traits_1_18_20[,c("trait_value","method","TRY_trait","genus species")]
trial.1 <- unique(trial)






#####################################


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







#Use Meghan's script for cleaning TRY data to clean BIEN data
library(tidyverse)
library(data.table)

theme_set(theme_bw(12))

#meghan's
setwd("C://Users/mavolio2/Dropbox/converge_diverge/datasets/Traits/Try Data Nov 2019")
setwd("C://Users/megha/Dropbox/converge_diverge/datasets/Traits/Try Data Nov 2019")
#kim's desktop
setwd('C:\\Users\\komatsuk\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\Traits\\Try Data Nov 2019')

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\Traits\\Try Data Nov 2019')



#read in try data
dat<-fread("7764.txt",sep = "\t",data.table = FALSE,stringsAsFactors = FALSE,strip.white = TRUE)



#generate list of units for ALL TRY traits
units <- dat%>%
  select(OriglName, OrigUnitStr, TraitName, UnitName)%>%
  unique()
# write.csv(units, 'TRY_all traits_units.csv')


#removing trait outliers
dat2<-dat%>%
  select(DatasetID,DataID, ObsDataID, ObservationID, AccSpeciesID, AccSpeciesName, TraitID, OriglName, TraitName, OrigValueStr, OrigUnitStr, StdValue, UnitName, ErrorRisk)%>%
  mutate(ErrorRisk2=ifelse(is.na(ErrorRisk), 0, ErrorRisk))%>%
  filter(ErrorRisk2<8)%>%
  filter(!is.na(TraitID))

#mering corre with try

#read in species links to merge try with core
key<-read.csv("corre2trykey.csv")%>%
  select(species_matched, AccSpeciesID, AccSpeciesName)%>%
  unique()

#merge list of species in Corre with TRY
dat3<-dat2%>%
  right_join(key)%>%
  select(-ErrorRisk, -ErrorRisk2)


####selecting desired continuous traits
#figuring out how clean traits are
tests<-dat3%>%
  filter(TraitID %in% c(3109))%>%
  #select(TraitID, OrigUnitStr, UnitName)%>%
  #unique()
  select(TraitID, UnitName, OrigUnitStr, OriglName, TraitName)%>%
  #filter(UnitName!="g/m2/d")%>%
  group_by(TraitID, UnitName, OrigUnitStr, OriglName, TraitName)%>%
  summarize(n=length(TraitID))




#testing consistent units for each trait and ranking traits by proirity
priority<-read.csv("trait_priority.csv")%>%
  rename(TraitID=TRY.trait.ID)

Traits_Units<-cont_traits%>%
  select(TraitID, TraitName, CleanTraitName, UnitName)%>%
  unique%>%
  mutate(Units=ifelse(CleanTraitName==3121, "g(W)/g(DM)", ifelse(CleanTraitName==3122, "g(W)/g(DM)", UnitName)))%>%
  select(-UnitName)%>%
  left_join(priority)

write.csv(Traits_Units, "For Franzi/ContTraitUnits.csv", row.names = F)

#subset out dead plants
#get list of dead plants
health<-dat%>%
  select(DatasetID,DataID, ObsDataID, AccSpeciesID, AccSpeciesName, TraitID, OriglName, TraitName, OrigValueStr, OrigUnitStr, StdValue, UnitName, ErrorRisk)%>%
  mutate(ErrorRisk2=ifelse(is.na(ErrorRisk), 0, ErrorRisk))%>%
  filter(ErrorRisk2<8)%>%
  filter(DataID==1961)%>%
  mutate(drop=ifelse(OrigValueStr=="Dead", 1, 0))%>%
  select(ObsDataID, drop)%>%
  unique()

table(health$OrigValueStr)

#merge with our list, there is no overalp
healthy<-cont_traits2%>%
  full_join(health)

##drop out trees that not seedlings
#read in which species are trees
treesp<-read.csv("species_families_trees_compelete.csv")%>%
  mutate(AccSpeciesName=species_matched)

#get list of tree observations that were made on seedlings
develop<-dat%>%
  select(DatasetID,DataID, ObsDataID, AccSpeciesID, AccSpeciesName, TraitID, OriglName, TraitName, OrigValueStr, OrigUnitStr, StdValue, UnitName, ErrorRisk)%>%
  mutate(ErrorRisk2=ifelse(is.na(ErrorRisk), 0, ErrorRisk))%>%
  filter(ErrorRisk2<8)%>%
  filter(DataID==413)%>%
  right_join(treesp)%>%
  filter(tree.non.tree=="tree")%>%
  mutate(drop=ifelse(OrigValueStr=="seedlings"|OrigUnitStr==0|OrigValueStr=="seedling"|OrigValueStr=="Seedling (0 - 1 y)",  0, 1))%>%
  select(ObsDataID, drop)%>%
  unique()

table(unique(develop$drop))

#merge to drop tree obseravtions that are not seedlings - none were.
developed<-cont_traits2%>%
  full_join(develop)

##drop out plant that were not measured in natural conditions - there is no overlap

#get list of observations that were not in natural settings
setting<-dat%>%
  select(DatasetID,DataID, ObsDataID, AccSpeciesID, AccSpeciesName, TraitID, OriglName, TraitName, OrigValueStr, OrigUnitStr, StdValue, UnitName, ErrorRisk)%>%
  mutate(ErrorRisk2=ifelse(is.na(ErrorRisk), 0, ErrorRisk))%>%
  filter(ErrorRisk2<8)%>%
  filter(DataID==327)%>%
  mutate(drop=ifelse(OrigValueStr=="Botanical garden"|OrigValueStr=="botanical garden (Bergius Botanical Garden, Stockholm, Sweden)"|OrigValueStr=="Botanical gardens, greenhouses and other atypical habitats"|OrigValueStr=="Chamber"|OrigValueStr=="climate chamber"|OrigValueStr=="Climate chamber"|OrigValueStr=="Climate Chamber"|OrigValueStr=="Climate chamber, non-limiting conditions, (cf. dataset reference)"|OrigValueStr=="climate chambers"|OrigValueStr=="Common Garden"|OrigValueStr=="Controlled climate chamber"|OrigValueStr=="controlled environment room"|OrigValueStr=="drought treatment"|OrigValueStr=="FACE"|OrigValueStr=="FE"|OrigValueStr=="C"|OrigValueStr=="Field Experiment"|OrigValueStr=="FW"|OrigValueStr=="G"|OrigValueStr=="GH"|OrigValueStr=="Glasshouse"|OrigValueStr=="Greehouse"|OrigValueStr=="Green house"|OrigValueStr=="greenhouse"|OrigValueStr=="Greenhouse"|OrigValueStr=="Greenhouse, grrowth container"|OrigValueStr=="groth chamber"|OrigValueStr=="growth-chamber"|OrigValueStr=="growth chamber"|OrigValueStr=="Growth chamber"|OrigValueStr=="Growth Chamber"|OrigValueStr=="growth chambers"|OrigValueStr=="Growth chambers"|OrigValueStr=="Growth exp"|OrigValueStr=="hydroponic"|OrigValueStr=="Irrigation"|OrigValueStr=="Irrigation and N fertilisation (100 kg/ha)"|OrigValueStr=="LAU_Ploughed/mown"|OrigValueStr=="LAU_Ploughed/mown and fertilized"|OrigValueStr=="mesocosm"|OrigValueStr=="mini-ecosystem"|OrigValueStr=="N"|OrigValueStr=="natural environment, high warming +4C, preccipitation ambient"|OrigValueStr=="natural environment, high warming +4C, preccipitation ambient -50%"|OrigValueStr=="natural environment, high warming +4C, preccipitation ambient +50%"|OrigValueStr=="natural environment, low warming +1.5C, preccipitation ambient"|OrigValueStr=="natural environment, low warming +1.5C, preccipitation ambient -50%"|OrigValueStr=="natural environment, low warming +1.5C, preccipitation ambient +50%"|OrigValueStr=="natural environment, medium warming +2.5C, preccipitation ambient"|OrigValueStr=="natural environment, medium warming +2.5C, preccipitation ambient -50%"|OrigValueStr=="natural environment, medium warming +2.5C, preccipitation ambient +50%"|OrigValueStr=="natural environment, no warming, preccipitation ambient -50%"|OrigValueStr=="natural environment, no warming, preccipitation ambient +50%"|OrigValueStr=="natural grassland, experimental nutrient NP addition"|OrigValueStr=="nutrient addition experiment"|OrigValueStr=="Open Top"|OrigValueStr=="open-top chamber"|OrigValueStr=="Open top chambers"|OrigValueStr=="OTC"|OrigValueStr=="plantation"|OrigValueStr=="PM"|OrigValueStr=="pot"|OrigValueStr=="Pot-grown"|OrigValueStr=="Pots outside"|OrigValueStr=="pots, outside in natural environment"|OrigValueStr=="shade houses"|OrigValueStr=="university campus"|OrigValueStr=="Uzbekistan: Irrigated desert land"|OrigValueStr=="VER_permanent extensively mown meadow"|OrigValueStr=="VER_permanent meadow mown and fertilized"|OrigValueStr=="VER_permanent meadows mown and fertilized"|OrigValueStr=="water stress experiment"|OrigValueStr=="water treatment", 1, 0))%>%
  select(ObsDataID, drop)%>%
  unique()

table(setting$drop)


#merge with out traits - there is no overalp
settings<-cont_traits2%>%
  full_join(setting)


#add info on genus family
fam<-read.csv("species_families.csv")

splist<-key%>%
  select(species_matched)%>%
  unique%>%
  left_join(fam)%>%
  separate(remove=F, species_matched, into=c("genus", "species"), sep=" ")%>%
  select(family, genus, species_matched)%>%
  left_join(treesp)

cont_traits2<-cont_traits%>%
  left_join(splist)%>%
  filter(tree.non.tree=="non-tree")%>%
  select(-tree.non.tree)


cont_traits3<-cont_traits2%>%
  select(DatasetID, ObservationID, family, genus, species_matched, CleanTraitName, StdValue)
###getting dataset to give to Frazni

#investigating problem traits
d453<-cont_traits3%>%#this dataset has 3 obs per plant, but no way to link leaves so we are averaging
  filter(DatasetID==453)%>%
  group_by(DatasetID, ObservationID, species_matched, CleanTraitName, family, genus)%>%
  summarise(StdValue=mean(StdValue))

d428<-cont_traits3%>% #this dataset has two height values per plant, we are taking the largest
  filter(DatasetID==428&CleanTraitName=="plant_height_vegetative"|DatasetID==428&CleanTraitName=="root_P")%>%
  group_by(DatasetID, ObservationID, species_matched, CleanTraitName, family, genus)%>%
  summarise(StdValue=max(StdValue))


cont_traits4<-cont_traits3%>%
  filter(DatasetID!=453)%>%
  mutate(remove=ifelse(DatasetID==428&CleanTraitName=="plant_height_vegetative", 1, 
                       ifelse(DatasetID==428&CleanTraitName=="root_P", 1, 0)))%>%
  filter(remove==0)%>%
  select(-remove)%>%
  bind_rows(d453)%>%
  bind_rows(d428)

cont_traits5<-cont_traits4%>%
  mutate(present=1)%>%
  group_by(DatasetID, ObservationID, species_matched, CleanTraitName)%>%
  summarize(n=sum(present))

probtraits<-subset(cont_traits5, n>1)%>%
  ungroup()%>%
  select(CleanTraitName, n)%>%
  unique()
#spread(CleanTraitName, StdValue)

ttraits<-cont_traits4%>%
  ungroup()%>%
  group_by(DatasetID, ObservationID, family, genus, species_matched)%>%
  spread(CleanTraitName, StdValue, fill=NA)


write.csv(ttraits, "For Franzi/TRY_trait_data_continuous.csv", row.names = F)
write.csv(cont_traits4, "For Franzi/TRY_trait_data_continuous_long.csv", row.names = F)




