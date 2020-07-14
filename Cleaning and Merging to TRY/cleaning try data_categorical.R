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
key<-read.csv("corre2trykey.csv")%>%
  select(species_matched, AccSpeciesID, AccSpeciesName)%>%
  unique()

length(unique(key$species_matched))

splist<-key%>%
  select(species_matched)%>%
  unique

dat3<-dat2%>%
  right_join(key)%>%
  select(-ErrorRisk, -ErrorRisk2)


##how many traits for sp
sdivtrt<-read.csv("TRY_traits_type_11252019.csv")

traitnum<-dat3%>%
  select(species_matched, TraitID)%>%
  unique()%>%
  group_by(TraitID)%>%
  summarise(nusp=length(species_matched))%>%
  right_join(sdivtrt)

# write.csv(traitnum, "try_traits_export_nov2019.csv", row.names=F)


###lifespan
trait59<-dat3%>%
  filter(TraitID==59&OrigValueStr!="")%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=="1"|OrigValueStr=="always annual"|OrigValueStr=="ann"|OrigValueStr=="annual"|OrigValueStr=="Annual"|OrigValueStr=="annual-winter annual"|OrigValueStr=="annuals"|OrigValueStr=="winter annual"|OrigValueStr=="winter annuals"|OriglName=="Plant phenology: Annual"&OrigValueStr=="yes"|OriglName=="Plant phenology: Perennial"&OrigValueStr=="no"|OrigValueStr=="summer annuals", "Annual", 
        ifelse(OrigValueStr=="2"|OrigValueStr=="1, 2"|OrigValueStr=="1,2"|OrigValueStr=="1-2"|OriglName=="Plant phenology: Biennial"&OrigValueStr=="yes"|OrigValueStr=="always annual, always biennial"|OrigValueStr=="always biennial"|OrigValueStr=="annual-winter annual, biennial"|OrigValueStr=="annual, Biennial"|OrigValueStr=="Annual, Biennial"|OrigValueStr=="annual/bieenial"|OrigValueStr=="annual/biennial"|OrigValueStr=='annual/bisannual'|OrigValueStr=="biannual"|OrigValueStr=="biasannual"|OrigValueStr=="biennial"|OrigValueStr=="sometimes annual, always biennial"|OrigValueStr=="winter annual-biennial"|OrigValueStr=="always annual, always biennial, always pluriennial-hapaxanthic"|OrigValueStr=="strict monocarpic bi-annuals and poly-annuals"|OrigValueStr=="Biennial", "Biennial", 
        ifelse(OriglName=="Plant phenology: Biennial"&OrigValueStr=="no", NA, "Perennial"))))%>%
  filter(!is.na(CleanTraitValue))%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(is.na(Biennial)&is.na(Perennial), "annual", ifelse(is.na(Annual)&is.na(Perennial)|Annual=="Annual"&Biennial=="Biennial"&is.na(Perennial), "biennial", "perennial")))%>%
  select(species_matched, CleanTraitValue)%>%
  mutate(CleanTraitName="lifespan", CleanTraitUnit=NA, source='TRY_59')

table(trait59$CleanTraitValue)


trait42<-dat3%>%
  filter(TraitID==42&OrigValueStr!=""&OrigValueStr!="?"&OriglName!="CONSENSUS")%>%
  mutate(CleanTraitValue=ifelse(OriglName=="aquatic"|OriglName=="carnivory"|OriglName=="Crop"|OriglName=="Ecological type"|OriglName=="final growth form 4 categories"&OrigValueStr=="herb"|OrigValueStr=="amphibiousubmerged"|OrigValueStr=="angiosperm"|OrigValueStr=="annual"|OrigValueStr=="Annual"|OrigValueStr=="aquatic"|OrigValueStr=="Aquatic"|OrigValueStr=="aquatic fresh water"|OrigValueStr=="aquatic, fresh water, floating"|OrigValueStr=="Aquatic_Epiphyte_Shrub_Tree_Vine_Herb"|OrigValueStr=="Aquatic_Herb"|OrigValueStr=="Aquatic_Shrub_Herb"|OrigValueStr=="Aquatic_Tree_Herb"|OriglName=="SGF"&OrigValueStr=="H"|OrigValueStr=="carnivore"|OrigValueStr=="CARNIVORE"|OrigValueStr=="carnivorous"|OrigValueStr=="chasmophyte"|OrigValueStr=="crop"|OrigValueStr=="crops"|OrigValueStr=="Colonizing"|OrigValueStr=="Columnar"|OrigValueStr=="Conical"|OrigValueStr=="Decumbent"|OrigValueStr=="Emergent attached to the substrate"|OrigValueStr=="Epiphyte"|OrigValueStr=="Erect"|OriglName=="shoot growth form"|OriglName=="Shape and Orientation"|OrigValueStr=="free"|OrigValueStr=="geophyte"|OrigValueStr=="gymnosperm"|OriglName=="parasite"|OriglName=="Parasitic"|OrigValueStr=="hydrohalophyte"|OrigValueStr=="Hydrophytes"|OrigValueStr=="Liana_Herb"|OrigValueStr=="HEMI-PARASITE"|OriglName=="Growth form"&OrigValueStr=="herb"|OriglName=="Growth Form (herb,shrub,tree,herbaceous vine,liana/woody vine)"&OrigValueStr=="herb"|OriglName=="Plant growth form"&OrigValueStr=="herb"|OriglName=="Plant growth form"&OrigValueStr=="herbaceous"|OriglName=="Plant growth form"&OrigValueStr=="herb/sub-shrub"|OriglName=="Plant growth form"&OrigValueStr=="herb/shrub"|OriglName=="GrowthForm"&OrigValueStr=="herb"|OriglName=="GROWTHFORM_STD"&OrigValueStr=="Herb"|OriglName=="GROWTHFORM_DIV"&OrigValueStr=="Herb"|OriglName=="GROWTHFORM_DIV"&OrigValueStr=="Shrub_Herb"|OriglName=="Life form"&OrigValueStr=="herbaceous"|OriglName=="Life form"&OrigValueStr=="Herbaceous perennial"|OrigValueStr=="herbaceous monocotyl"|OriglName=="plant growth form"&OrigValueStr=="herb"|OriglName=="GROWTHFORM_ORG"&OrigValueStr=="Herb"|OrigValueStr=="no"|OrigValueStr=="No"|OrigValueStr=="Nano-chamaephyte"|OrigValueStr=="non-woody"|OrigValueStr=="Nontree"|OrigValueStr=="Nonvascular"|OrigValueStr=="Parasite"|OrigValueStr=="Parasite_Herb"|OrigValueStr=="parasitic"|OrigValueStr=="perennial"|OrigValueStr=="Perennial"|OriglName=="Life form: geophyte"|OriglName=="Life form: epiphyte/parasite"|OriglName=="Plant form: Non-distinctive"|OriglName=="Plant form: prostrate"|OriglName=="Plant form: cushion"|OriglName=="Plant form: open"|OrigValueStr=="xerohalophyte"|OrigValueStr=="xerophyte"|OrigValueStr=="weed"|OrigValueStr=="weedy"|OrigValueStr=="Shrub_Vine_Herb"|OrigValueStr=="Tree_Vine_Shrub_Herb"|OrigValueStr=="Vine_Herb"|OrigValueStr=="Vine_Shrub_Tree_Herb"|OrigValueStr=="Vine_Shrub_Tree"|OrigValueStr=="Vine_Shrub_Herb"|OrigValueStr=="terrestrial"|OrigValueStr=="Terrestrial Herb"|OrigValueStr=="therophyte"|OrigValueStr=="Thicket Forming"|OrigValueStr=="moss"|OrigValueStr=="Multiple Stem"|OrigValueStr=="psammophile"|OrigValueStr=="pteridophyte"|OrigValueStr=="Single Crown"|OrigValueStr=="Single Stem"|OrigValueStr=="PS"|OrigValueStr=="rhiz"|OrigValueStr=="rhizomatous"|OrigValueStr=="Rhizomatous"|OrigValueStr=="SC"|OrigValueStr=="Scap"|OrigValueStr=="scrub"|OrigValueStr=="Rhiz"|OrigValueStr=="Herb/Aquatic"|OrigValueStr=="herb/non-woody"|OrigValueStr=="herb/palmoid/non-woody"|OrigValueStr=="Herb/Shrub"|OrigValueStr=="herb/shrub/climber/non-woody/woody"|OrigValueStr=="herb/shrub/non-woody/woody"|OrigValueStr=="herb/shrub/palmoid/non-woody/woody"|OrigValueStr=="Herb/Shrub/Subshrub"|OrigValueStr=="Herb/Shrub/Tree"|OrigValueStr=="Herb/Shrub/Vine"|OrigValueStr=="non-succulent"|DatasetID==79|DatasetID==412|DatasetID==380|OrigValueStr=="creeper"|DatasetID==57|DatasetID==92|DatasetID==50|DatasetID==322|DatasetID==436|DatasetID==422|DatasetID==218|DatasetID==37|DatasetID==236|DatasetID==443|OriglName=="Stem succulent"|DatasetID==200|DatasetID==339|DatasetID==1|DatasetID==460|DatasetID==110|DatasetID==251|DatasetID==4, NA,
    ifelse(OrigValueStr=="b H"|OrigValueStr=="a T"|OrigValueStr=="forb"|OrigValueStr=="herbaceous legume"|OrigValueStr=="annual forb"|OriglName=="GF"&OrigValueStr=="H"|OrigValueStr=="Forb"|OrigValueStr=="forb (herbaceous, with or without woody base)"|OrigValueStr=="h"|OrigValueStr=="D"|OrigValueStr=="HS"|OrigValueStr=="HSA"|OrigValueStr=="HSL"|OrigValueStr=="HSLT"|OrigValueStr=="HST"|OrigValueStr=="Epiphyte_Herb"|OrigValueStr=="Epiphyte_Liana_Tree_Shrub_Herb"|OrigValueStr=="Epiphyte_Tree_Vine_Shrub_Herb"|OrigValueStr=="Forb/herb"|OrigValueStr=="Forbs"|OrigValueStr=='Forb/herb, Subshrub'|OrigValueStr=="Forb/herb, Vine"|OrigValueStr=="Forb/herb, Shrub, Subshrub"|OrigValueStr=="Forb/herb, Shrub, Subshrub, Vine"|OrigValueStr=="H"|OrigValueStr=="herbs"|OrigValueStr=="perennial leguminous herb"|OriglName=="GrowthFormCleartext"&OrigValueStr=="herb"|OriglName=="GrowthformCleartext"&OrigValueStr=="herb"|OriglName=="Life form"&OrigValueStr=="herb"|OriglName=='Life Form'&OrigValueStr=="herb"|OrigValueStr=="herbaceous dicotyl"|OrigValueStr=="Herbaceous Dicot"|OriglName=="type"&OrigValueStr=="herb"|OriglName=="growth_form   TRY"&OrigValueStr=="herb"|OriglName=="Plant Growth Form"&OrigValueStr=="Herb"|OriglName=="Plant Growth Form"&OrigValueStr=="herbaceous"|OriglName=="Plant growth forb"&OrigValueStr=="Herbaceous"|OrigValueStr=="perennial forb"|OriglName=="Life form: forb"|OrigValueStr=="Tree_Vine_Aquatic_Shrub_Herb"|OrigValueStr=="Tree_Vine_Herb"|OrigValueStr=="variable forb"|OrigValueStr=="n hyd"|OrigValueStr=="n Hyd"|OrigValueStr=="m Hel"|OrigValueStr=="rept"|OrigValueStr=="Rept"|OrigValueStr=="herb/aquatic/non-woody"|OrigValueStr=="herb/hemiparasitic/non-woody"|OriglName=="Succulence"|OriglName=="Succulence index"|OriglName=="Succulent"|OrigValueStr=="C"|OrigValueStr=="succulent"|OrigValueStr=="Succulent"|OrigValueStr=="succulent/non-woody"|OrigValueStr=="C"|OrigValueStr=="cactus"|OrigValueStr=="stem-succulent", "Forb", 
     ifelse(OrigValueStr=="Tree"|OrigValueStr=="tree"|OrigValueStr=="Absence"|OrigValueStr=="shrub"|OrigValueStr=="woody plant"|OrigValueStr=="f P"|OrigValueStr=="S"|OrigValueStr=="SH"|OrigValueStr=="ST"|OrigValueStr=="T"|OrigValueStr=="t"|OrigValueStr=="Woody"|OrigValueStr=="W"|OrigValueStr=="Shrub"|OrigValueStr=="subshrub (woody <1m)"|OrigValueStr=="sh"|OrigValueStr=="t"|OrigValueStr=="P"|OrigValueStr=="c C"|OrigValueStr=="Chaemaephyte"|OrigValueStr=="conifer"|OrigValueStr=="Conifers"|OrigValueStr=="d z"|OrigValueStr=="d Z"|OrigValueStr=="Deciduous shrub or tree"|OrigValueStr=="Dwarf shrub"|OrigValueStr=="e N"|OrigValueStr=="Epiphyte_Shrub_Herb"|OrigValueStr=="Epiphyte_Vine_Tree_Shrub_Herb"|OrigValueStr=="erect dwarf shrub"|OrigValueStr=="evergreen shrub or tree"|OrigValueStr=="large shrub"|OrigValueStr=="low to high shrub"|OrigValueStr=="palmoid"|OrigValueStr=="prostrate dwarf shrub"|OriglName=="Life form: erect dwarf shrub"|OriglName=="Life form: prostrate dwarf shrub"|OriglName=="Life form: shrub"|OriglName=="Life form: tree"|OrigValueStr=="Woody Liana"|OrigValueStr=="Woody"|OrigValueStr=="Woody evergreen"|OrigValueStr=="Woody deciduous"|OrigValueStr=="woody at base"|OrigValueStr=="woody"|OrigValueStr=="TREE"|OrigValueStr=="Tree (deciduous)"|OrigValueStr=="Tree (evergreen)"|OrigValueStr=="tree / shrub"|OrigValueStr=="Tree shrub intermediate"|OrigValueStr=="Tree, Shrub"|OrigValueStr=="Tree, Subshrub, Shrub"|OrigValueStr=="tree/palmoid/woody"|OrigValueStr=="Tree/Treelet"|OrigValueStr=="tree/woody"|OrigValueStr=="Tree_Shrub"|OrigValueStr=="trees"|OrigValueStr=="trees/tree"|OrigValueStr=="trees/tree"|OrigValueStr=="Shrub, Subshrub"|OrigValueStr=="Shrub, Subshrub, Tree"|OrigValueStr=="Shrub, Tree"|OrigValueStr=="Shrub,Subshrub"|OrigValueStr=="shrub/palmoid/woody"|OrigValueStr=="Shrub/Subshrub"|OrigValueStr=="shrub/tree"|OrigValueStr=="Shrub/Tree"|OrigValueStr=="Shrub/Tree intermediate"|OrigValueStr=="shrub/tree/palmoid/woody"|OrigValueStr=="Shrub/Tree/Subshrub"|OrigValueStr=="	Shrub/Tree/Treelet"|OrigValueStr=="shrub/tree/woody"|OrigValueStr=="shrub/woody"|OrigValueStr=="Shrub_Tree"|OrigValueStr=="shrub|tree"|OrigValueStr=="shrubs"|OrigValueStr=="small tree"|OrigValueStr=="Small_Tree"|OrigValueStr=="sub-shrub"|OrigValueStr=="Sub-Shrub (Chamaephyte)"|OrigValueStr=="subshrub"|OrigValueStr=="Subshrub"|OrigValueStr=="Subshrub, Shrub"|OrigValueStr=="Subshrub, Shrub, Tree", "Woody", 
     ifelse(OrigValueStr=="grass"|OrigValueStr=="sedge"|OrigValueStr=="G"|OrigValueStr=="annual grass"|OrigValueStr=="Bunch"|OrigValueStr=="grass (Poaceae only)"|OrigValueStr=="g"|OrigValueStr=="se"|OrigValueStr=="rus"|OrigValueStr=="C3 grass"|OrigValueStr=="C4 grass"|OrigValueStr=="Caesp"|OrigValueStr=="caesp"|OrigValueStr=="cereal"|OrigValueStr=="forage grass"|OrigValueStr=="graminoid"|OrigValueStr=="GRAMINOID"|OrigValueStr=="graminoid/aquatic"|OrigValueStr=="graminoid/aquatic/non-woody"|OrigValueStr=="graminoid/non-woody"|OrigValueStr=="Graminoids"|OrigValueStr=="Graminoids Tussock"|OrigValueStr=="Graminoid"|OrigValueStr=="Grass"|OrigValueStr=="grass (clonal)"|OrigValueStr=="grasslike"|OrigValueStr=="Herbaceous Monocot"|OrigValueStr=="pasture grass"|OrigValueStr=="perennial graminoid"|OrigValueStr=="perennial grass"|OrigValueStr=="Perennial grass"|OrigValueStr=="prairie grass"|OriglName=="Life form: graminoid"|OriglName=="Low Growing Grass"|OrigValueStr=="Sedge"|OrigValueStr=="SEDGE","Graminoid",
     ifelse(OriglName=="climber"|OriglName=="ClimbingMode"|OrigValueStr=="climber"|OrigValueStr=="Vine"|OrigValueStr=="V"|OrigValueStr=="twiner/climber."|OrigValueStr=="L"|OrigValueStr=="C+Sc"|OrigValueStr=="climber"|OrigValueStr=="climber or creeper"|OrigValueStr=="climber/non-woody"|OrigValueStr=="climber/parasitic"|OrigValueStr=="climber/vine"|OrigValueStr=="climber/woody"|OrigValueStr=="Climbing"|OrigValueStr=="Climber"|OrigValueStr=="Lianas and climbers"|OrigValueStr=="g L"|OrigValueStr=="liana"|OrigValueStr=="Liana"|OrigValueStr=="lianas"|OrigValueStr=="Lianas (wody climbers)"|OrigValueStr=="lianas/Woody Liana"|OrigValueStr=="Lianna"|OrigValueStr=="Liana_Shrub"|OrigValueStr=="Liana_Vine"|OrigValueStr=="Liana_Vine_Herb"|OrigValueStr=="vine"|OriglName=="Life form: climber"|OriglName=="Life form: liana"|OriglName=="Plant form: climbing"|OrigValueStr=="Vines (non-woody climbers)"|OrigValueStr=="Herb_Liana_Vine"|OrigValueStr=="Shrub_Liana_Vine"|OrigValueStr=="Shrub_Vine"|OrigValueStr=="Tree_Shrub_Liana_Herb_Vine"|OrigValueStr=="Vine"|OrigValueStr=="herb/climber/non-woody"|OrigValueStr=="herb/climber/parasitic/non-woody"|OrigValueStr=="herb/climber/woody"|OrigValueStr=="Herb/Liana", "Vine",
     ifelse(OrigValueStr=="F"|OrigValueStr=="fern"|OrigValueStr=="Fern"|OrigValueStr=="FERN"|OrigValueStr=="FERN ALLY"|OrigValueStr=="Fern or fern ally"|OrigValueStr=="fern/non-woody"|OrigValueStr=="FERNALLY"|OrigValueStr=="Ferns"|OrigValueStr=="Ferns and allies (Lycophytes)"|OrigValueStr=="M"|OriglName=="Life form: fern/fern ally"|OrigValueStr=="club moss"|OrigValueStr=="CLUBMOSS"|OrigValueStr=="Club moss", "Fern",
            NA)))))))%>%
  filter(!is.na(CleanTraitValue))

#can't get these to drop correctly.
#|OriglName=="PlantGrowthFormConsolidated"&OrigValueStr=="shrub/herb"

table(trait42$CleanTraitValue)

trait42_test<-trait42%>%
  select(CleanTraitValue)%>%
  unique()

trait42_fern<-trait42%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Fern=="Fern", "fern", 999))%>%
  filter(CleanTraitValue!=999)%>%
  select(species_matched, CleanTraitValue)

trait42_forb<-trait42%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Forb=="Forb"&is.na(Graminoid)&is.na(Fern)&is.na(Woody)&is.na(Vine), "forb", 999))%>%
  filter(CleanTraitValue!=999)%>%
  select(species_matched, CleanTraitValue)

trait42_gram<-trait42%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Graminoid=="Graminoid"&is.na(Forb)&is.na(Woody)&is.na(Vine)&is.na(Fern),"graminiod", 999))%>%
  filter(CleanTraitValue!=999)%>%
  select(species_matched, CleanTraitValue)

trait42_vine<-trait42%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Vine=="Vine"&is.na(Forb)&is.na(Woody)&is.na(Graminoid)&is.na(Fern),"vine", 999))%>%
  filter(CleanTraitValue!=999)%>%
  select(species_matched, CleanTraitValue)

trait42_woody<-trait42%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Woody=="Woody"&is.na(Forb)&is.na(Vine)&is.na(Graminoid)&is.na(Fern),"woody", 999))%>%
  filter(CleanTraitValue!=999)%>%
  select(species_matched, CleanTraitValue)

trait42_problem<-trait42%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  group_by(species_matched)%>%
  summarise(n=length(CleanTraitValue))%>%
  filter(n!=1)%>%
  select(-n)

trait42_probelm3<-trait42_problem[c(1:40),]%>%
  mutate(CleanTraitValue=ifelse(species_matched=="Amorpha canescens", "woody",
                          ifelse(species_matched=="Artemisia annua", "forb",
                          ifelse(species_matched=="Artemisia frigida", "forb",
                          ifelse(species_matched=="Artemisia gmelinii", "forb",
                          ifelse(species_matched=="Asparagus officinalis", "forb",
                          ifelse(species_matched=="Atriplex canescens", "woody",
                          ifelse(species_matched=="Atriplex patula", "forb",
                          ifelse(species_matched=="Chenopodium glaucum", "forb",
                          ifelse(species_matched=="Chrysocephalum apiculatum", "forb",
                          ifelse(species_matched=="Comandra umbellata", "forb",
                          ifelse(species_matched=="Convolvulus arvensis", "vine",
                          ifelse(species_matched=="Convolvulus erubescens", "vine",
                          ifelse(species_matched=="Coreopsis lanceolata", "forb",
                          ifelse(species_matched=="Crataegus monogyna", "woody",
                          ifelse(species_matched=="Cuscuta glomerata", "vine",
                          ifelse(species_matched=="Cyrilla racemiflora","woody",
                          ifelse(species_matched=="Dalea purpurea", "forb",
                          ifelse(species_matched=="Dryas integrifolia", "woody",
                          ifelse(species_matched=="Dryas octopetala", "woody",
                          ifelse(species_matched=="Dryopteris carthusiana", "fern",
                          ifelse(species_matched=="Elymus repens", "graminoid",
                          ifelse(species_matched=="Equisetum arvense", "fern",
                          ifelse(species_matched=="Erigeron canadensis", "forb",
                          ifelse(species_matched=="Euphorbia corollata", "forb",
                          ifelse(species_matched=="Euphorbia dentata", "forb",
                          ifelse(species_matched=="Fallopia convolvulus", "vine",
                          ifelse(species_matched=="Fallopia scandens", "vine",
                          ifelse(species_matched=="Galium aparine", "forb",
                          ifelse(species_matched=="Galium verum", "forb",
                          ifelse(species_matched=="Gutierrezia sarothrae", "woody",
                          ifelse(species_matched=="Harrimanella hypnoides", "forb",
                          ifelse(species_matched=="Helianthemum nummularium", "woody",
                          ifelse(species_matched=="Hypochaeris radicata", "forb",
                          ifelse(species_matched=="Krascheninnikovia ceratoides", "woody",
                          ifelse(species_matched=="Lathyrus pratensis", "forb",
                          ifelse(species_matched=="Lespedeza capitata", "forb",
                          ifelse(species_matched=="Lespedeza juncea", "woody",
                          ifelse(species_matched=="Linnaea borealis","woody",
                          ifelse(species_matched=="Lonicera japonica", "vine",
                          ifelse(species_matched=="Lonicera periclymenum", "vine", 999)))))))))))))))))))))))))))))))))))))))))

trait42_probelm4<-trait42_problem[c(41:72),]%>%
  mutate(CleanTraitValue=ifelse(species_matched=="Mollugo verticillata", "vine",
                          ifelse(species_matched=="Moneses uniflora", "woody",
                          ifelse(species_matched=="Oenothera biennis", "forb",
                          ifelse(species_matched=="Orthilia secunda", "woody",
                          ifelse(species_matched=="Parthenocissus inserta", "vine",
                          ifelse(species_matched=="Parthenocissus quinquefolia", "vine",
                          ifelse(species_matched=="Phryma leptostachya", "forb",
                          ifelse(species_matched=="Phytolacca americana", "forb",
                          ifelse(species_matched=="Pimelea trichostachya", "woody",
                          ifelse(species_matched=="Plantago coronopus", "forb",
                          ifelse(species_matched=="Portulaca oleracea", "forb",
                          ifelse(species_matched=="Pyrola elliptica", "forb",
                          ifelse(species_matched=="Rosa multiflora", "woody",
                          ifelse(species_matched=="Rubus idaeus", "woody",
                          ifelse(species_matched=="Rubus vestitus", "woody",
                          ifelse(species_matched=="Salix repens", "woody",
                          ifelse(species_matched=="Salsola kali", 'woody',
                          ifelse(species_matched=="Solanum americanum", "forb",
                          ifelse(species_matched=="Solanum dulcamara", "vine",
                          ifelse(species_matched=="Stellaria media", "forb",
                          ifelse(species_matched=="Talinum polygaloides", "forb",
                          ifelse(species_matched=="Thymus praecox", "woody",
                          ifelse(species_matched=="Tofieldia pusilla", "forb",
                          ifelse(species_matched=="Toxicodendron diversilobum", "vine",
                          ifelse(species_matched=="Tribulus terrestris","forb",
                          ifelse(species_matched=="Typha angustifolia", "forb",
                          ifelse(species_matched=="Typha latifolia", "forb",
                          ifelse(species_matched=="Vicia americana", "vine",
                          ifelse(species_matched=="Vicia cracca", "vine",
                          ifelse(species_matched=="Vicia tetrasperma", "vine",
                          ifelse(species_matched=="Vicia villosa", "vine", NA))))))))))))))))))))))))))))))))
              
trait42_probelm2<-trait42%>%
  right_join(trait42_problem)%>%
  select(species_matched, CleanTraitValue, DatasetID)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)

  
trait42_test_messy<-trait42%>%
  filter(TraitID==42&OrigValueStr!=""&OrigValueStr!="?"&OriglName!="CONSENSUS")%>%
  select(OriglName, OrigValueStr, species_matched, CleanTraitValue)%>%
  unique()

trait42_clean<-rbind(trait42_fern, trait42_forb, trait42_gram, trait42_vine, trait42_woody, trait42_probelm3, trait42_probelm4)%>%
  mutate(CleanTraitName="lifeform", CleanTraitUnit=NA, source='TRY_42')


##c3/c4 photosynthesis
trait22<-dat3%>%
  filter(TraitID==22)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=="3"|OrigValueStr=="c3"|OrigValueStr=="C3", "C3",
                         ifelse(OrigValueStr=="C4", "C4",
                         ifelse(OrigValueStr=="CAM", "CAM", NA))))%>%
  filter(!is.na(CleanTraitValue))%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(is.na(C3)&is.na(C4), "CAM", ifelse(is.na(CAM)&is.na(C3), "C4", ifelse(is.na(C4)&is.na(CAM), "C3", NA))))%>%
  filter(!is.na(CleanTraitValue))%>%
  select(species_matched, CleanTraitValue)%>%
  mutate(CleanTraitName="photo_pathway", CleanTraitUnit=NA, source='TRY_22')



##stem support
trait1188_clean <- dat3%>%
  filter(TraitID==1188)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=='self', 'self-supporting', ifelse(OrigValueStr=='creeping', 'prostrate', ifelse(OrigValueStr=='procumbent', 'prostrate', ifelse(OrigValueStr=='twining', 'climbing', ifelse(OrigValueStr=='tendrils', 'climbing', ifelse(OrigValueStr=='scrambling', 'climbing', OrigValueStr)))))))%>%
  mutate(CleanTraitName='stem_support', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)%>%
  unique()%>%
  mutate(CleanTraitValue2=ifelse(CleanTraitValue=='self-supporting', 'self_supporting', CleanTraitValue))%>%
  select(-CleanTraitValue)%>%
  spread(CleanTraitValue2, CleanTraitValue2, fill='fix')%>%
  mutate(CleanTraitValue2=paste(climbing, decumbent, prostrate, self_supporting, sep=','))%>%
  mutate(CleanTraitValue=ifelse(CleanTraitValue2 %in% c('climbing,decumbent,fix,fix','climbing,fix,fix,fix'), 'climbing', ifelse(CleanTraitValue2 %in% c('fix,decumbent,fix,fix', 'fix,decumbent,prostrate,fix'), 'decumbent', ifelse(CleanTraitValue2=='fix,fix,prostrate,fix', 'prostrate', 'self-supporting'))))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)%>%
  mutate(source='TRY_1188')

table(trait1188_clean$CleanTraitValue)



#clonal growth

#208 reproductive type
trait208<-dat3%>%
  filter(TraitID==208)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr %in% c("vegetative & generative", "vegetative", "vegetatively ", "vegetative", "seed_and_vegetative"), 'yes', 'no'))%>%
  filter(!is.na(CleanTraitValue))%>%
  select(species_matched, CleanTraitValue)%>%
  unique()

#341 clonal vs non-clonal
trait341<-dat3%>%
  filter(TraitID==341)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=="no", "no", "yes"))%>%
  select(species_matched, CleanTraitValue)%>%
  filter(!is.na(CleanTraitValue))%>%
  unique()

#329 is very similar to 341, but the low categories too
trait329<-dat3%>%
  filter(TraitID==329)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=="Absence", "no", "yes"))%>%
  select(species_matched, CleanTraitValue)%>%
  filter(!is.na(CleanTraitValue))%>%
  unique()

#334 has a range of spread
trait334<-dat3%>%
  filter(TraitID==334)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr %in% c("", "dispersable"), NA, 'yes'))%>% 
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  filter(!is.na(CleanTraitValue))

#344 has a range of spread
trait344<-dat3%>%
  filter(TraitID==344)%>%
  mutate(CleanTraitValue=ifelse(DatasetID==92|OrigValueStr=="no particuliar mating system", NA,
         ifelse(OrigValueStr=="none"|OrigValueStr=="NotClonal", "no", 'yes')))%>%
  select(species_matched, CleanTraitValue)%>%
  filter(!is.na(CleanTraitValue))%>%
  unique()

##357 role of the clonal organ in growth
trait357 <- dat3%>%
  filter(TraitID==357)%>%
  filter(OrigValueStr!='')%>%
  rename(CleanTraitValue=OrigValueStr)%>%
  mutate(CleanTraitName='clonality', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue2=paste(necessary, additive, regenerative, none, sep=','))%>%
  mutate(CleanTraitValue=ifelse(CleanTraitValue2 %in% c('NA,additive,NA,NA', 'NA,additive,NA,none', 'NA,additive,regenerative,NA', 'NA,additive,regenerative,none', 'necessary,additive,NA,NA', 'necessary,additive,NA,none', 'necessary,additive,regenerative,NA', 'necessary,additive,regenerative,none', 'necessary,NA,NA,NA', 'necessary,NA,NA,none'), 'yes', 'no'))%>%
  select(species_matched, CleanTraitValue)

##609 reproductive type
trait609<-dat3%>%
  filter(TraitID==609)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=="No"|OrigValueStr=="Yes", NA, ifelse(OrigValueStr %in% c("seed and vegetative", "vegetative"), 'yes', 'no')))%>%
  filter(!is.na(CleanTraitValue))%>%
  select(species_matched, CleanTraitValue)%>%
  unique()

##613 vegetative spread rate
trait613 <- dat3%>%
  filter(TraitID==613)%>%
  filter(!is.na(OrigValueStr))%>%
  rename(CleanTraitValue2=OrigValueStr)%>%
  mutate(CleanTraitName='clonality', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue2, CleanTraitUnit)%>%
  unique()%>%
  spread(CleanTraitValue2, CleanTraitValue2)%>%
  mutate(CleanTraitValue=ifelse(is.na(Rapid)&is.na(Moderate)&is.na(Slow), 'no', 'yes'))%>%
  select(species_matched, CleanTraitValue)

clonality <- rbind(trait208, trait329, trait341, trait334, trait344, trait357, trait609, trait613)%>%
  unique()%>%
  spread(key=CleanTraitValue, value=CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(yes=='yes', 'yes', 'no'))%>%
  mutate(CleanTraitName='clonal', CleanTraitUnit=NA, source='TRY_208:329:334:341:344:357:609:613')%>%
  filter(!is.na(CleanTraitValue))%>%
  select(CleanTraitName, species_matched, CleanTraitValue, source, CleanTraitUnit)



##dispersal
trait28<-dat3%>%
  filter(TraitID==28)%>%
  mutate(problem=substr(OrigValueStr,35, 40))%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr==0|OrigValueStr==0.8|OrigValueStr==0.5|OrigValueStr==0.2|OrigValueStr==0.6|OrigValueStr=="Water?"|OrigValueStr=="vegetative dispersule"|OrigValueStr=="unknown"|OrigValueStr=="seed contamination"|OrigValueStr=="other"|OrigValueStr=="one-seeded generative dispersule"|OrigValueStr=="speirochor"|OrigValueStr=="multi-seeded generative dispersule"|OrigValueStr=="hay making machinery"|OrigValueStr=="hay transport"|OrigValueStr=="hemerochor"|OrigValueStr=='germinule'|OrigValueStr=="generative dispersule"|OrigValueStr=="external"|OrigValueStr=="erosion material"|OrigValueStr=="Dispersal prevented"|OrigValueStr=="Dispersal no"|OrigValueStr=="Disp"|OrigValueStr=="car or other vehicle"|OrigValueStr==2|OrigValueStr==3|OrigValueStr=="harvesting"|OrigValueStr=="Combination: animal+unassisted"|OrigValueStr=="Combination: methods originating from parent plant+animal"|OrigValueStr=="Combination: water+wind+animal"|OrigValueStr=="Combination: wind+animal+unassisted"|OrigValueStr=="Combination: wind+unassisted"|OrigValueStr=="Combination: wind+water"|OrigValueStr=="clothes and footwear"|ObsDataID==27123749|ObsDataID==27123751|OrigValueStr=="Wind Animals"|OrigValueStr=="Seeds are produced below ground level"|OrigValueStr=="commerce"|OrigValueStr=="Combination: wind+animal"|problem==", plan", NA,
                         ifelse(OriglName=="wind.disp"&OrigValueStr==1|OrigValueStr=="wind"|OrigValueStr=="Wind"|OriglName=="disp.mode.wind"&OrigValueStr==1|OrigValueStr=="wind/long-distance"|OrigValueStr=="wind-dispersed (with wing, hairs or bristles to provide air-resistance)"|OrigValueStr=="Dispersal wind"|OrigValueStr=="Diaspore is rolled along ground surface by wind"|OrigValueStr=="Diaspore is propelled by action of wind on the plant structure"|OrigValueStr=="Diaspore is blown by wind"|OrigValueStr=="chamaechor"|OrigValueStr=="boleochor"|OrigValueStr=="a-wind"|OrigValueStr=="Anemo"|OrigValueStr=="anemochory"|OrigValueStr=="Anemochory"|OrigValueStr=="Anemochory: Big and round seeds rolling on the ground, pushed by the wind"|OrigValueStr=="Anemochory: Small seeds with pappus or very light seeds (ex,  Crepis sp. or Orchis sp.)"|OrigValueStr=="Anemochory: Stems move with the wind, helping for seed dispersion (ex,  Papaver sp.)"|OrigValueStr=="meteorochor", "Wind",
                      ifelse(OriglName=="disp.mode.animal"&OrigValueStr==1, "Animal",
                         ifelse(OrigValueStr=="water"|OrigValueStr=="Water"|OriglName=="disp.mode.water"&OrigValueStr==1|OrigValueStr=="Wetting by rain or dew"|OrigValueStr=="standing fresh water"|OrigValueStr=="shaken fresh water"|OrigValueStr=="rainwash"|OrigValueStr=="ombrochor"|OrigValueStr=="nautochor "|OrigValueStr=="hydrochory"|OrigValueStr=="Floating in freshwater currents"|OrigValueStr=="nautochor", "Water", 
                         ifelse(OriglName=="disp.mode.gravity"&OrigValueStr==1|OrigValueStr=="unassisted/short-distance"|OrigValueStr=="Unassisted and/or methods originating from parent plant"|OrigValueStr=="unassisted (no morphological structures aiding dispersal)"|OrigValueStr=="unassisted"|OrigValueStr=="Unassisted"|OrigValueStr=="tumbling"|OrigValueStr=="Seeds drop to the ground close to or beneath the parent plant"|OrigValueStr=="non specialized"|OrigValueStr==" Methods originating from parent plant or diaspore"|OrigValueStr=="Gravity"|OrigValueStr=="Explosive mechanism"|OrigValueStr=="explosive mechanism"|OrigValueStr=="explosive"|OrigValueStr=="blastochor"|OrigValueStr=="Barochory"|OrigValueStr=="ballochor "|OrigValueStr=="Autochory"|OrigValueStr=="autochory"|OrigValueStr=="autochor"|OrigValueStr=="unspecialised"|OrigValueStr=="Methods originating from parent plant or diaspore", "Unassisted", "Animal"))))))%>%
  filter(!is.na(CleanTraitValue))%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Animal=="Animal"&is.na(Wind)&is.na(Unassisted)&is.na(Water), "animal",
                         ifelse(Animal=='Animal'&is.na(Water)&is.na(Wind)&Unassisted=="Unassisted", "animal",
                         ifelse(Wind=="Wind"&is.na(Animal)&is.na(Unassisted)&is.na(Water),"wind",
                         ifelse(Wind=='Wind'&is.na(Water)&is.na(Animal)&Unassisted=="Unassisted", "wind",
                         ifelse(Water=="Water"&is.na(Animal)&is.na(Wind)&is.na(Unassisted), "water", 
                         ifelse(Water=='Water'&is.na(Animal)&is.na(Wind)&Unassisted=="Unassisted", "water",
                         ifelse(Unassisted=="Unassisted"&is.na(Animal)&is.na(Wind)&is.na(Water), "unassisted", "combination"))))))))%>%
  select(species_matched, CleanTraitValue)%>%
  mutate(CleanTraitName="dispersal_mode", CleanTraitUnit="NA", source='TRY_28')


# #43 leaf type  -- so many of these are wrong, do we just need people to google each species individually?
# trait43 <- dat3%>%
#   filter(TraitID==43)%>%
#   filter(!is.na(OrigValueStr))%>%
#   rename(CleanTraitValue=OrigValueStr)%>%
#   mutate(CleanTraitName='leaf_type', CleanTraitUnit=NA, source='TRY_43')%>%
#   mutate(CleanTraitValue2=ifelse(CleanTraitValue %in% c('aphyllous', 'photosynthetic stem'), 'none', ifelse(CleanTraitValue %in% c('broad', 'broad-leaved', 'broadleaf', 'broadleaved') | (OriglName=='Leaf type: broad'&CleanTraitValue=='yes') | (OriglName=='Broadleaf/Conifer'&CleanTraitValue=='B'), 'broad', ifelse(CleanTraitValue %in% c('microphylle', 'mikrophylle'), 'microphyll', ifelse(CleanTraitValue=='narrowleaved', 'narrow', ifelse(CleanTraitValue %in% c('needle', 'needle-leaf', 'needle-leaved', 'needleleaf', 'needleleaved'), 'needle', ifelse(CleanTraitValue %in% c('scale', 'scale-leaf', 'scale-like', 'scale-shaped'), 'scale', NA)))))))%>%
#   select(species_matched, CleanTraitName, CleanTraitValue2, CleanTraitUnit)%>%
#   filter(!is.na(CleanTraitValue2))%>%
#   unique()%>%
#   spread(CleanTraitValue2, CleanTraitValue2)%>%
#   mutate(c1=ifelse(is.na(broad), 0, 1), c2=ifelse(is.na(narrow), 0, 1), c3=ifelse(is.na(none), 0, 1), c4=ifelse(is.na(needle), 0, 1), c5=ifelse(is.na(scale), 0, 1), c6=ifelse(is.na(microphyll), 0, 1), drop=c1+c2+c3+c4+c5+c6)%>%
#   filter(drop==1)%>% #dropping species that have multiple values, removes 98 species
#   select(-c1,-c2,-c3,-c4,-c5,-c6,-drop)%>%
#   mutate(CleanTraitValue=ifelse(broad=='broad', 'broad', ifelse(narrow=='narrow', 'narrow', ifelse(needle=='needle', 'needle', ifelse(scale=='scale', 'scale', ifelse(microphyll=='microphyll', 'microphyll', 'none'))))))
#   select(species_matched, CleanTraitValue)

#17 leaf compoundness
trait17 <- dat3%>%
  filter(TraitID==17)%>%
  filter(!is.na(OrigValueStr))%>%
  rename(CleanTraitValue=OrigValueStr)%>%
  mutate(CleanTraitName='leaf_compoundness', CleanTraitUnit=NA, source='TRY_17')%>%
  mutate(CleanTraitValue2=ifelse(CleanTraitValue %in% c('C', 'composite', 'compound', 'palmately compound', 'pinnately compound', 'trifoliolate'), 'compound', ifelse(CleanTraitValue %in% c('S', 'simple'), 'simple', NA)))%>%
  select(CleanTraitName, species_matched, CleanTraitName, CleanTraitValue2, CleanTraitUnit, source)%>%
  filter(!is.na(CleanTraitValue2))%>%
  unique()%>%
  spread(CleanTraitValue2, CleanTraitValue2)%>%
  mutate(c1=ifelse(is.na(compound), 0, 1), c2=ifelse(is.na(simple), 0, 1), drop=c1+c2)%>%
  filter(drop==1)%>% #dropping species that have multiple values, removes 65 species
  select(-c1,-c2,-drop)%>%
  mutate(CleanTraitValue=ifelse(compound=='compound', 'compound', 'simple'))%>%
  mutate(CleanTraitValue=ifelse(is.na(CleanTraitValue), 'simple', CleanTraitValue))%>%
  select(CleanTraitName, species_matched, CleanTraitValue, source, CleanTraitUnit)

#213, 347 mating system; needs update if we want to keep this
trait347 <- dat3%>%
  filter(TraitID==347|TraitID==213)%>%
  filter(!is.na(OrigValueStr))%>%
  rename(CleanTraitValue=OrigValueStr)

#335 phenology; needs update if we want to keep this
trait335 <- dat3%>%
  filter(TraitID==335)%>%
  filter(!is.na(OrigValueStr))%>%
  rename(CleanTraitValue=OrigValueStr)

#29 pollination; needs update if we want to keep this
trait29 <- dat3%>%
  filter(TraitID==29)%>%
  filter(!is.na(OrigValueStr))%>%
  rename(CleanTraitValue=OrigValueStr)%>%
  mutate(CleanTraitName='pollination', CleanTraitUnit=NA, source='TRY_29')%>%
  mutate(CleanTraitValue2=ifelse(CleanTraitValue %in% c('anemogamous/entomogamous', 'Animals Wind', 'autogamous/entomogamous', 'mixed wind/insect pollinated'), 'combination', ifelse(CleanTraitValue %in% c('bee', 'bees', 'bees, bumblebees, wasps, bombylides, syrphids', 'bees, bumble bees, wasps, bombylides, syrphids', 'bees, butterflies', 'beetles, flies, syrphids, wasps, medium tongued bees', 'bumble bees', 'bumblebees, butterfflies', 'butterflies', 'bees, toung < 7 mm', 'bumble bees, lepidoptera', 'bumblebees, lepidoptera', 'butterflies, long tongued bees, syrphids', 'flies', 'entomogamous', 'flies, bees', 'flies, beetles', 'general insect', 'hymenopteres', 'insect', 'insect pollinated', 'insects always', 'insects often', 'insects the rule', 'insects unknown', 'lepidoptera, bees', 'lepidoptera, bumble bees', 'moth/butterfly', 'moths', 'moths, hymenoptera', 'pollination animals', 'short tongued bees, syrphids, flies, beetles', 'Short tongued bees, syrphids, flies, beetles', 'short tounged bees, syrphids, muscids, beetles', 'syrphids', 'syrphids, bees'), 'animal', ifelse(CleanTraitValue %in% c('pollination wind', 'wind', 'Wind', 'wind always', 'wind often', 'wind pollinated', 'wind the rule'), 'wind', ifelse(CleanTraitValue %in% c('hydrogamous'), 'water', ifelse(CleanTraitValue %in% c('cleistogamy the rule', 'selfed', 'selfing always', 'selfing often', 'selfing possible', 'selfing the rule'), 'self', NA))))))%>%
  select(CleanTraitName, species_matched, CleanTraitName, CleanTraitValue2, CleanTraitUnit, source)%>%
  filter(!is.na(CleanTraitValue2))%>%
  unique()%>%
  spread(CleanTraitValue2, CleanTraitValue2)%>%
  mutate(CleanTraitValue=ifelse(animal=="animal"&is.na(wind)&is.na(water)&is.na(self), "animal",
                                ifelse(wind=="wind"&is.na(animal)&is.na(water)&is.na(self),"wind",
                                       ifelse(self=="self"&is.na(animal)&is.na(wind)&is.na(water), "self",
                                              ifelse(water=="water"&is.na(animal)&is.na(wind)&is.na(self), "water", "combination")))))%>%
  select(CleanTraitName, species_matched, CleanTraitValue, source, CleanTraitUnit)



##mycorrhizal traits
mycorr<-read.csv('mycorr_Species_match_for_Kim.csv')%>%
  select(species_matched, Mycorrhizal.type)%>%
  unique()%>%
  mutate(CleanTraitValue=ifelse(Mycorrhizal.type %in% c('AM', 'EcM', 'ErM', 'OM', 'EcM-AM', 'double_AM_EcM', 'NM-AM', 'NM-AM, rarely EcM'), 'yes', ifelse(Mycorrhizal.type=='NM', 'no', ifelse(Mycorrhizal.type=='uncertain', 'uncertain', ''))))%>%
  mutate(CleanTraitName="mycorrhizal", CleanTraitUnit="", source='FungalRoot')%>%
  select(CleanTraitName, species_matched, CleanTraitValue, source, CleanTraitUnit)%>%
  # spread(CleanTraitValue, CleanTraitValue)%>%
  # mutate(CleanTraitValue=ifelse(yes=="yes"&is.na(no)&is.na(uncertain), "yes",
  #                               ifelse(is.na(yes)&no=='no'&is.na(uncertain), "no",
  #                                      ifelse(is.na(yes)&is.na(no)&uncertain=='uncertain', "uncertain",
  #                                             'check'))))%>% #a few were double listed as NA and either yes, no, or uncertain
  filter(!is.na(CleanTraitValue))%>%
  select(CleanTraitName, species_matched, CleanTraitValue, source, CleanTraitUnit)

mycorrType<-read.csv('mycorr_Species_match_for_Kim.csv')%>%
  select(species_matched, Mycorrhizal.type)%>%
  unique()%>%
  mutate(CleanTraitValue=ifelse(Mycorrhizal.type=='AM', 'arbuscular', ifelse(Mycorrhizal.type=='EcM', 'ecto', ifelse(Mycorrhizal.type=='ErM', 'ericaceous', ifelse(Mycorrhizal.type=='OM', 'orchidaceous', ifelse(Mycorrhizal.type=='EcM-AM', 'double_AM_EcM', ifelse(Mycorrhizal.type=='NM', 'none', ifelse(Mycorrhizal.type=='NM-AM', 'facultative_AM', ifelse(Mycorrhizal.type=='NM-AM, rarely EcM', 'facultative_AM_EcM', ifelse(Mycorrhizal.type=='uncertain', 'uncertain', ''))))))))))%>%
  filter(!is.na(CleanTraitValue))%>%
  mutate(CleanTraitName="mycorrhizal_type", CleanTraitUnit="", source='FungalRoot')%>%
  select(CleanTraitName, species_matched, CleanTraitValue, source, CleanTraitUnit)

#n-fixation
nFix <- read.csv('CoRRE_TRY_species_list_N-fixers.csv')%>%
  select(species_matched, fixer)%>%
  mutate(CleanTraitValue=ifelse(fixer==1, 'yes', 'no'))%>%
  mutate(CleanTraitName="n_fixation", CleanTraitUnit="", source='Werner 2014')%>%
  select(CleanTraitName, species_matched, CleanTraitValue, source, CleanTraitUnit)

rhizobial <- read.csv('CoRRE_TRY_species_list_N-fixers.csv')%>%
  select(species_matched, fixer, act)%>%
  mutate(CleanTraitValue=ifelse(fixer==1 & act==0, 'yes', 'no'))%>%
  mutate(CleanTraitName="rhizobial", CleanTraitUnit="", source='Werner 2014')%>%
  select(CleanTraitName, species_matched, CleanTraitValue, source, CleanTraitUnit)

actinorhizal <- read.csv('CoRRE_TRY_species_list_N-fixers.csv')%>%
  select(species_matched, act)%>%
  mutate(CleanTraitValue=ifelse(act==1, 'yes', 'no'))%>%
  mutate(CleanTraitName="actinorhizal", CleanTraitUnit="", source='Werner 2014')%>%
  select(CleanTraitName, species_matched, CleanTraitValue, source, CleanTraitUnit)
  


###combine traits
traitsCat <- rbind(trait59, trait42_clean, trait22, mycorr, mycorrType, trait1188_clean, clonality, trait28, trait17, trait29, nFix, rhizobial, actinorhizal)%>%
  mutate(trait_source=paste(CleanTraitValue, source, sep='::'))%>%
  unique()%>%
  select(-CleanTraitUnit, -source, -CleanTraitValue)%>%
  spread(key=CleanTraitName, value=trait_source)%>%
  separate(lifeform, c('lifeform', 'lifeform_source'), sep='::')%>%
  separate(lifespan, c('lifespan', 'lifespan_source'), sep='::')%>%
  separate(clonal, c('clonal', 'clonal_source'), sep='::')%>%
  separate(dispersal_mode, c('dispersal_mode', 'dispersal_mode_source'), sep='::')%>%
  separate(leaf_compoundness, c('leaf_compoundness', 'leaf_compoundness_source'), sep='::')%>%
  separate(mycorrhizal, c('mycorrhizal', 'mycorrhizal_source'), sep='::')%>%
  separate(mycorrhizal_type, c('mycorrhizal_type', 'mycorrhizal_type_source'), sep='::')%>%
  separate(photo_pathway, c('photo_pathway', 'photo_pathway_source'), sep='::')%>%
  separate(pollination, c('pollination', 'pollination_source'), sep='::')%>%
  separate(stem_support, c('stem_support', 'stem_support_source'), sep='::')%>%
  separate(n_fixation, c('n_fixation', 'n_fixation_source'), sep='::')%>%
  separate(rhizobial, c('rhizobial', 'rhizobial_source'), sep='::')%>%
  separate(actinorhizal, c('actinorhizal', 'actinorhizal_source'), sep='::')%>%
  full_join(splist)%>%
  left_join(read.csv('species_families.csv'))


# write.csv(traitsCat, 'categorical_traits_tofillin.csv', row.names=F)


# write.csv(traits_cat, "C://Users/mavolio2/Dropbox/SDiv_sCoRRE_shared/CoRRE - community and anpp data/TRY_trait_data_categorical.csv", row.names = F)

# summary_cat<-traits_cat%>%
#   group_by(CleanTraitName, TraitCategory, TraitType)%>%
#   summarize(n=length(CleanTraitValue))%>%
#   mutate(per.sp=(n/1954)*100)
# 
# summary_cont<-traits_cont%>%
#   group_by(CleanTraitName, TraitCategory, TraitType)%>%
#   summarize(n=length(CleanTraitValue))%>%
#   mutate(per.sp=(n/1954)*100)
# 
# summary_all<-summary_cat%>%
#   bind_rows(summary_cont)
# 
# ggplot(data=summary_all, aes(x=reorder(CleanTraitName, -per.sp), y=per.sp))+
#   geom_point()+
#   ylab("percent of species")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank())+
#   xlab("")+
#   scale_y_continuous(limits=c(0,100))
# 
# #barplots by trait type
# ggplot(data=subset(summary_all, TraitCategory=="Reproduction"), aes(x=1, y=per.sp))+
#   geom_bar(stat = "identity")+
#  facet_wrap(~CleanTraitName)+
#   scale_y_continuous(limits=c(0,100))+
#   ylab("Percent of Species")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank())+
#   xlab("")