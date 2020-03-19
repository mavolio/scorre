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
  select(DatasetID,DataID, ObsDataID, AccSpeciesID, AccSpeciesName, TraitID, OriglName, TraitName, OrigValueStr, OrigUnitStr, StdValue, UnitName, ErrorRisk)%>%
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
  filter(TraitID %in% c(3120))%>%
  select(TraitID, OrigUnitStr, UnitName)%>%
  unique()
  select(TraitID, UnitName, TraitName)%>%
  #filter(UnitName!="g/m2/d")%>%
  group_by(TraitID, UnitName, TraitName)%>%
   summarize(n=length(TraitID))

#subsetting out traits and naming them
cont_traits<-dat3%>%
  filter(TraitID %in% c(4, 6, 9, 12, 13, 14, 15, 26, 27, 40, 41, 44, 45, 46, 47, 48, 50, 51, 52, 53, 55, 56, 57, 58, 66, 77, 80, 82, 83, 84, 106, 111, 138, 145, 146, 185, 186, 269, 270, 363, 475, 570, 614, 683, 1080, 1104, 1781, 2809, 3106, 3107, 3108, 3109, 3110, 3111, 3112, 3113, 3114, 3115, 3116, 3117,3120, 3121, 3122))%>%
  mutate(remove=ifelse(TraitID==48&UnitName=='', 1, ifelse(TraitID==3107&UnitName=='cm', 1, ifelse(TraitID==53&UnitName=='g/m2/d',1, ifelse(TraitID==4&UnitName=='', 1, ifelse(TraitID==3116&UnitName=='', 1, ifelse(TraitID==3122&UnitName=='', 1, ifelse(TraitID==3121&UnitName=='', 1, 0))))))))%>%#remove problem data
  filter(remove==0)%>%
  select(-remove)%>%
  mutate(CleanTraitName=ifelse(TraitID==4, 'stem_spec_density', 
                        ifelse(TraitID==6, 'rooting_depth', 
                        ifelse(TraitID==9, 'root:shoot', 
                        ifelse(TraitID==12, 'leaf_longevity', 
                        ifelse(TraitID==13, 'leaf_C',
                        ifelse(TraitID==14, 'leaf_N',
                        ifelse(TraitID==15, 'leaf_P',
                        ifelse(TraitID==26, 'seed_dry_mass', 
                        ifelse(TraitID==27, 'seed_length', 
                        ifelse(TraitID==41, 'dark_resp_rate',
                        ifelse(TraitID==44, 'leaf_K',
                        ifelse(TraitID==45, 'stomatal_conductance',
                        ifelse(TraitID==46, 'leaf_thickness', 
                        ifelse(TraitID==47, 'LDMC', 
                        ifelse(TraitID==48, 'leaf_density', 
                        ifelse(TraitID==53, 'photosynthesis_rate',
                        ifelse(TraitID==55, 'leaf_dry_mass', 
                        ifelse(TraitID==56, 'leaf_N:P', 
                        ifelse(TraitID==66, 'seed_terminal_velocity',
                        ifelse(TraitID==77, 'RGR', 
                        ifelse(TraitID==80, 'root_N', 
                        ifelse(TraitID==82, 'root_density', 
                        ifelse(TraitID==83, 'root_diameter', 
                        ifelse(TraitID==84, 'root_C', 
                        ifelse(TraitID==111, 'leaf_transp_rate',
                        ifelse(TraitID==138, 'seed_number',
                        ifelse(TraitID==145, 'leaf_width', 
                        ifelse(TraitID==146, 'leaf_C:N', 
                        ifelse(TraitID==186, 'Vc_max',
                        ifelse(TraitID==269, 'J_max',
                        ifelse(TraitID==363, 'root_dry_mass',
                        ifelse(TraitID==683, 'root_P', 
                        ifelse(TraitID==1080, 'SRL',
                        ifelse(TraitID==2809, 'seedbank_duration',
                        ifelse(TraitID==3106, 'plant_height_vegetative', 
                        ifelse(TraitID==3107, 'plant_height_generative', 
                        ifelse(TraitID==3116, 'SLA', 
                        ifelse(TraitID==3110, 'leaf_area',
                        ifelse(TraitID==3120, 'water_content', 
                        TraitID))))))))))))))))))))))))))))))))))))))))

#testing consistent units for each trait
test2<-cont_traits%>%
  select(TraitID, TraitName, CleanTraitName, UnitName)%>%
  unique

#ranking traits by proirity
priority<-read.csv("trait_priority.csv")%>%
  rename(TraitID=TRY.trait.ID)

cont_traits2<-cont_traits%>%
  left_join(priority)

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
  right_join(health)

##drop out trees that not seedlings
#read in which species are trees
treesp<-read.csv("species_families_trees_compelete.csv")%>%
  rename(AccSpeciesName=species_matched)

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
  right_join(develop)

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
  right_join(setting)


#add info on genus family
fam<-read.csv("species_families.csv")

splist<-key%>%
  select(species_matched)%>%
  unique%>%
  left_join(fam)%>%
  separate(remove=F, species_matched, into=c("genus", "species"), sep=" ")%>%
  select(family, genus, species_matched)

cont_traits3<-cont_traits2%>%
  left_join(splist)

###getting dataset to give to Frazni
cont_traits4<-cont_traits3%>%
  select(ObservationID, AccSpeciesName, family, genus, species_matched, CleanTraitName, StdValue)%>%
  spread(CleanTraitName, StdValue)


write.csv(traits_cont, "C://Users/mavolio2/Dropbox/SDiv_sCoRRE_shared/CoRRE - community and anpp data/TRY_trait_data_continuous.csv", row.names = F)
