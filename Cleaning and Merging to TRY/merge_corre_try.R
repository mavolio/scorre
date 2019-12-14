library(tidyverse)
library(data.table)
library(Hmisc)
library(utf8)

#from Habacuc

#meghan's
setwd("C:/Users/mavolio2/Dropbox/converge_diverge/datasets/LongForm/fixing species names")
setwd("C:/Users/megha/Dropbox/converge_diverge/datasets/LongForm/fixing species names")

#kim's
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm\\fixing species names')

checkcorre<-read.csv("../SpeciesRelativeAbundance_Nov2019.csv")

#load clean taxonomy for try
taxdat <- read_csv("taxon_updates.csv")

#get rid of species matched by taxize
ii <- is.na(taxdat$species_Taxonstand) & !is.na(taxdat$species_taxize)

taxdat$species_matched[ii] <- taxdat$species[ii]

#select species and matched species columns
taxdat %>%
  select(species, species_matched) -> taxdat

#load corre
corre <- read_csv("CoRRE_SpList_Sept2019.csv")
#preprocessing to utf8
corre %>%
  mutate(species= as_utf8(genus_species))->corre
#capitalize all records
corre$species <- Hmisc::capitalize(tolower(corre$genus_species))

#load taxonomy for corre
taxcorre <- read_csv("Species_to_check_cleaned_2.csv")%>%
  filter(remove==0)#this filters out only to sp., unknowns, and non-vasular plants except ferns.
#get rid of taxize matches - we manually did all this
# ii <- is.na(taxcorre$species_Taxonstand) & !is.na(taxcorre$species_taxize)
# 
# taxcorre$species_matched[ii] <- taxcorre$species[ii]
#select submitted name and matched name
taxcorre %>%
  select(species, species_matched) -> taxcorre
#join corre to updated taxonomy
corre2 <- left_join(corre, taxcorre)%>%
  select(-species)%>%
  na.omit()

#read try 
try <- fread("TryAccSpecies.txt",sep = "\t",data.table = FALSE,stringsAsFactors = FALSE,strip.white = TRUE)
#preprocess try
try %>%
  mutate(species= as_utf8(AccSpeciesName))->try
#capitalize records
try$species <- Hmisc::capitalize(tolower(try$AccSpeciesName))#why are we doing this?
try$match<-ifelse(try$species==try$AccSpeciesName,1,0)

#join try to updated taxonomu
try <- left_join(try,taxdat, by = c("AccSpeciesName"="species"))%>%
  select(-species)

#join corre to try
corre2try <- left_join(corre2,try, by="species_matched")%>%
  unique()

# write_csv(corre2try, path = "corre2trykey.csv")

#make comma separted row to submit to try 

try_list <- corre2try[["AccSpeciesID"]][!is.na(corre2try$AccSpeciesID)]

# write_delim(x = as.data.frame(t(try_list)), "list_for_try.csv",delim = ",",col_names = FALSE)




###generating list for phylogeney
#want to include all columns, and indicate where a moss/lichen, plus include anything that is identified to genera
taxcorreAll <- read.csv("Species_to_check_cleaned_2.csv")%>%
  filter(remove!=3)%>% #this filters out unknowns, but keeps mosses/lichens and anything that was IDed to genus.
  select(species, species_matched, remove)%>%
  mutate(type=ifelse(remove==2, 'moss/lichen', ifelse(remove==1, 'identified genus', 'identified species')))%>%
  select(-remove)
#join corre to updated taxonomy
correTaxonomyAll <- left_join(corre, taxcorreAll)%>%
  select(-species)%>%
  na.omit()%>%
  unique()
# write.csv(correTaxonomyAll, 'CoRRE_TRY_species_list.csv')
