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

##how many traits for sp?
sdivtrt_con<-read.csv("TRY_traits_type_11252019.csv")%>%
  filter(trait_type=="con")

####how many speices have at least 2 observed continouous trait?
sptrait<-dat3%>%
  right_join(sdivtrt_con)%>%
  select(TraitID, species_matched)%>%
  unique()%>%
  group_by(species_matched)%>%
  summarize(n=length(TraitID))%>%
  full_join(splist)
  
#have 2 or more traits
two_or_more<-sptrait%>%
  filter(n>1)

one<-sptrait%>%
  filter(n==1)

none<-sptrait%>%
  filter(is.na(n))
 
write.csv(none, "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/CoRRE data/species_no_traits.csv")

##how many traits for sp
sdivtrt<-read.csv("TRY_traits_type_11252019.csv")

traitnum<-dat3%>%
  select(species_matched, TraitID)%>%
  unique()%>%
  group_by(TraitID)%>%
  summarise(nusp=length(species_matched))%>%
  right_join(sdivtrt)

# write.csv(traitnum, "try_traits_export_nov2019.csv", row.names=F)
trait_test<-dat3%>%
  filter(TraitID==159)


###cleaning life history traits
trait59<-dat3%>%
  filter(TraitID==59&OrigValueStr!="")%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=="1"|OrigValueStr=="always annual"|OrigValueStr=="ann"|OrigValueStr=="annual"|OrigValueStr=="Annual"|OrigValueStr=="annual-winter annual"|OrigValueStr=="annuals"|OrigValueStr=="winter annual"|OrigValueStr=="winter annuals"|OriglName=="Plant phenology: Annual"&OrigValueStr=="yes"|OriglName=="Plant phenology: Perennial"&OrigValueStr=="no"|OrigValueStr=="summer annuals", "Annual", 
        ifelse(OrigValueStr=="2"|OrigValueStr=="1, 2"|OrigValueStr=="1,2"|OrigValueStr=="1-2"|OriglName=="Plant phenology: Biennial"&OrigValueStr=="yes"|OrigValueStr=="always annual, always biennial"|OrigValueStr=="always biennial"|OrigValueStr=="annual-winter annual, biennial"|OrigValueStr=="annual, Biennial"|OrigValueStr=="Annual, Biennial"|OrigValueStr=="annual/bieenial"|OrigValueStr=="annual/biennial"|OrigValueStr=='annual/bisannual'|OrigValueStr=="biannual"|OrigValueStr=="biasannual"|OrigValueStr=="biennial"|OrigValueStr=="sometimes annual, always biennial"|OrigValueStr=="winter annual-biennial"|OrigValueStr=="always annual, always biennial, always pluriennial-hapaxanthic"|OrigValueStr=="strict monocarpic bi-annuals and poly-annuals"|OrigValueStr=="Biennial", "Biennial", 
        ifelse(OriglName=="Plant phenology: Biennial"&OrigValueStr=="no", NA, "Perennial"))))%>%
  filter(!is.na(CleanTraitValue))

table(trait59$CleanTraitValue)

trait59_test<-trait59%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  select(Annual, Biennial, Perennial)%>%
  unique

trait59_clean<-trait59%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(is.na(Biennial)&is.na(Perennial), "Annual", ifelse(is.na(Annual)&is.na(Perennial)|Annual=="Annual"&Biennial=="Biennial"&is.na(Perennial), "Biennial", "Perennial")))%>%
  select(species_matched, CleanTraitValue)%>%
  mutate(CleanTraitName="lifespan", CleanTraitUnit=NA)


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
  mutate(CleanTraitValue=ifelse(Fern=="Fern", "Fern", 999))%>%
  filter(CleanTraitValue!=999)%>%
  select(species_matched, CleanTraitValue)

trait42_forb<-trait42%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Forb=="Forb"&is.na(Graminoid)&is.na(Fern)&is.na(Woody)&is.na(Vine), "Forb", 999))%>%
  filter(CleanTraitValue!=999)%>%
  select(species_matched, CleanTraitValue)

trait42_gram<-trait42%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Graminoid=="Graminoid"&is.na(Forb)&is.na(Woody)&is.na(Vine)&is.na(Fern),"Graminiod", 999))%>%
  filter(CleanTraitValue!=999)%>%
  select(species_matched, CleanTraitValue)

trait42_vine<-trait42%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Vine=="Vine"&is.na(Forb)&is.na(Woody)&is.na(Graminoid)&is.na(Fern),"Vine", 999))%>%
  filter(CleanTraitValue!=999)%>%
  select(species_matched, CleanTraitValue)

trait42_woody<-trait42%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Woody=="Woody"&is.na(Forb)&is.na(Vine)&is.na(Graminoid)&is.na(Fern),"Woody", 999))%>%
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
  mutate(CleanTraitValue=ifelse(species_matched=="Amorpha canescens", "Woody",
                          ifelse(species_matched=="Artemisia annua", "Forb",
                          ifelse(species_matched=="Artemisia frigida", "Forb",
                          ifelse(species_matched=="Artemisia gmelinii", "Forb",
                          ifelse(species_matched=="Asparagus officinalis", "Forb",
                          ifelse(species_matched=="Atriplex canescens", "Woody",
                          ifelse(species_matched=="Atriplex patula", "Forb",
                          ifelse(species_matched=="Chenopodium glaucum", "Forb",
                          ifelse(species_matched=="Chrysocephalum apiculatum", "Forb",
                          ifelse(species_matched=="Comandra umbellata", "Forb",
                          ifelse(species_matched=="Convolvulus arvensis", "Vine",
                          ifelse(species_matched=="Convolvulus erubescens", "Vine",
                          ifelse(species_matched=="Coreopsis lanceolata", "Forb",
                          ifelse(species_matched=="Crataegus monogyna", "Woody",
                          ifelse(species_matched=="Cuscuta glomerata", "Vine",
                          ifelse(species_matched=="Cyrilla racemiflora","Woody",
                          ifelse(species_matched=="Dalea purpurea", "Forb",
                          ifelse(species_matched=="Dryas integrifolia", "Woody",
                          ifelse(species_matched=="Dryas octopetala", "Woody",
                          ifelse(species_matched=="Dryopteris carthusiana", "Fern",
                          ifelse(species_matched=="Elymus repens", "Graminiod",
                          ifelse(species_matched=="Equisetum arvense", "Fern",
                          ifelse(species_matched=="Erigeron canadensis", "Forb",
                          ifelse(species_matched=="Euphorbia corollata", "Forb",
                          ifelse(species_matched=="Euphorbia dentata", "Forb",
                          ifelse(species_matched=="Fallopia convolvulus", "Vine",
                          ifelse(species_matched=="Fallopia scandens", "Vine",
                          ifelse(species_matched=="Galium aparine", "Forb",
                          ifelse(species_matched=="Galium verum", "Forb",
                          ifelse(species_matched=="Gutierrezia sarothrae", "Woody",
                          ifelse(species_matched=="Harrimanella hypnoides", "Forb",
                          ifelse(species_matched=="Helianthemum nummularium", "Woody",
                          ifelse(species_matched=="Hypochaeris radicata", "Forb",
                          ifelse(species_matched=="Krascheninnikovia ceratoides", "Woody",
                          ifelse(species_matched=="Lathyrus pratensis", "Forb",
                          ifelse(species_matched=="Lespedeza capitata", "Forb",
                          ifelse(species_matched=="Lespedeza juncea", "Woody",
                          ifelse(species_matched=="Linnaea borealis","Woody",
                          ifelse(species_matched=="Lonicera japonica", "Vine",
                          ifelse(species_matched=="Lonicera periclymenum", "Vine", 999)))))))))))))))))))))))))))))))))))))))))

trait42_probelm4<-trait42_problem[c(41:72),]%>%
  mutate(CleanTraitValue=ifelse(species_matched=="Mollugo verticillata", "Vine",
                          ifelse(species_matched=="Moneses uniflora", "Woody",
                          ifelse(species_matched=="Oenothera biennis", "Forb",
                          ifelse(species_matched=="Orthilia secunda", "Woody",
                          ifelse(species_matched=="Parthenocissus inserta", "Vine",
                          ifelse(species_matched=="Parthenocissus quinquefolia", "Vine",
                          ifelse(species_matched=="Phryma leptostachya", "Forb",
                          ifelse(species_matched=="Phytolacca americana", "Forb",
                          ifelse(species_matched=="Pimelea trichostachya", "Woody",
                          ifelse(species_matched=="Plantago coronopus", "Forb",
                          ifelse(species_matched=="Portulaca oleracea", "Forb",
                          ifelse(species_matched=="Pyrola elliptica", "Forb",
                          ifelse(species_matched=="Rosa multiflora", "Woody",
                          ifelse(species_matched=="Rubus idaeus", "Woody",
                          ifelse(species_matched=="Rubus vestitus", "Woody",
                          ifelse(species_matched=="Salix repens", "Woody",
                          ifelse(species_matched=="Salsola kali", 'Woody',
                          ifelse(species_matched=="Solanum americanum", "Forb",
                          ifelse(species_matched=="Solanum dulcamara", "Vine",
                          ifelse(species_matched=="Stellaria media", "Forb",
                          ifelse(species_matched=="Talinum polygaloides", "Forb",
                          ifelse(species_matched=="Thymus praecox", "Woody",
                          ifelse(species_matched=="Tofieldia pusilla", "Forb",
                          ifelse(species_matched=="Toxicodendron diversilobum", "Vine",
                          ifelse(species_matched=="Tribulus terrestris","Forb",
                          ifelse(species_matched=="Typha angustifolia", "Forb",
                          ifelse(species_matched=="Typha latifolia", "Forb",
                          ifelse(species_matched=="Vicia americana", "Vine",
                          ifelse(species_matched=="Vicia cracca", "Vine",
                          ifelse(species_matched=="Vicia tetrasperma", "Vine",
                          ifelse(species_matched=="Vicia villosa", "Vine", NA))))))))))))))))))))))))))))))))
              
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
  mutate(CleanTraitName="lifeform", CleanTraitUnit=NA)

##c3/c4 photosynthesis

trait22<-dat3%>%
  filter(TraitID==22)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=="3"|OrigValueStr=="c3"|OrigValueStr=="C3", "C3",
                         ifelse(OrigValueStr=="C4", "C4",
                         ifelse(OrigValueStr=="CAM", "CAM", NA))))%>%
  filter(!is.na(CleanTraitValue))

table(trait22$OrigValueStr)  

trait22_clean<-trait22%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(is.na(C3)&is.na(C4), "CAM", ifelse(is.na(CAM)&is.na(C3), "C4", ifelse(is.na(C4)&is.na(CAM), "C3", NA))))%>%
  filter(!is.na(CleanTraitValue))%>%
  select(species_matched, CleanTraitValue)%>%
  mutate(CleanTraitName="photo_pathway", CleanTraitUnit=NA)

##mycorrhizal traits
trait_1433_1030_clean<-dat3%>%
  filter(TraitID==1433|TraitID==1030)%>%
  group_by(species_matched)%>%
  summarise(CleanTraitValue=mean(StdValue))%>%
  mutate(CleanTraitName="mycorrhizal_percent_colonization", CleanTraitUnit="%")

#all these species are in trait7
# trait_1433_1030_sp<-trait_1433_1030%>%
#   filter(CleanTriatValue!=0)%>%
#   select(species_matched)

trait7<-dat3%>%
  filter(TraitID==7)%>%
  mutate(CleanTraitValue=ifelse(OriglName=="Stable_AM_loss_likelihood"|OriglName=="AM_Stable_likelihood"|OriglName=="AM_retained_likelihood"|OriglName=="Labile_likelihood"|OriglName=="AM_lost_likelihood"|OrigValueStr=="non-ectomycorrhizal"|OrigValueStr=="?va", NA, 
        ifelse(OrigValueStr=="absent"|OrigValueStr=="0"|OrigValueStr=="Absent"|OrigValueStr=="Non"|OrigValueStr=="no"|OrigValueStr=="No"|OrigValueStr=="Non Mycorr"|OrigValueStr=="N", "Non Mycorr", "Mycorr")))%>%
  filter(!is.na(CleanTraitValue))
                        
table(trait7$CleanTraitValue)

trait7_clean<-trait7%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Mycorr=="Mycorr", "Mycorr", NA))%>%
  filter(!is.na(CleanTraitValue))%>%
  select(species_matched, CleanTraitValue)%>%
  mutate(CleanTraitName="mycorrhizal", CleanTraitUnit=NA)


##plant functional type
##this is not immediately usable. there are a lot of acronyms I dont' know.
# trait197<-dat3%>%
#   filter(TraitID==197)
# 
# table(trait197$OrigValueStr)

##leaf area - merging different traits that all correspond to leaf area
#do everything with the StdValue, which is converted to mm2
#DECISION: combine all leaf area data into one clean variable (see below); some of the regressions are very poor, however most are good (>0.7)
#NOTE: this needs some fixing (see pipeline document for information)

#filter outliers
traitLeafAreaGenus <- dat3%>%
  filter(TraitID %in% c(3114, 3108, 3110, 3112, 3109, 3111, 3113))%>% #all data related to leaf areas
  separate(species_matched, into=c('genus', 'species'), sep=' ')%>%
  group_by(genus, TraitID)%>%
  summarise(genus_mean=mean(StdValue), genus_sd=sd(StdValue))%>%
  ungroup()

traitLeafAreaSpp <- dat3%>%
  filter(TraitID %in% c(3114, 3108, 3110, 3112, 3109, 3111, 3113))%>% #all data related to leaf areas
  group_by(species_matched, TraitID)%>%
  summarise(spp_mean=mean(StdValue), spp_sd=sd(StdValue))%>%
  ungroup()

trait3108_3109_3110_3111_3112_3113_3114_clean <- dat3%>%
  filter(TraitID %in% c(3114, 3108, 3110, 3112, 3109, 3111, 3113))%>% #all data related to leaf areas
  separate(species_matched, into=c('genus', 'species'), sep=' ', remove=F)%>%
  left_join(traitLeafAreaGenus)%>%
  left_join(traitLeafAreaSpp)%>%
  mutate(genus_zscore=abs((StdValue-genus_mean)/genus_sd), spp_zscore=abs((StdValue-spp_mean)/spp_sd))%>% #calculate z-scores for genera and species
  filter(genus_zscore<4, spp_zscore<2)%>%
  group_by(DatasetID, species_matched, TraitID)%>%
  summarise(DatasetValue=mean(StdValue))%>% #averaging by ways to measure leaf area and species within each dataset
  ungroup()%>%
  group_by(species_matched, TraitID)%>%
  summarise(SppValue=mean(DatasetValue))%>% #averaging by ways to measure leaf area and species across datasets
  ungroup()%>%
  group_by(species_matched)%>%
  summarise(CleanTraitValue=mean(SppValue))%>% #averaging by species across datasets and ways to measure leaf area
  ungroup()%>%
  mutate(CleanTraitName='leaf_area', CleanTraitUnit='mm2')
  
  # #getting averages within each species for each trait type, are they comparable?
  # traitLeafAreaTest <- trait3108_3109_3110_3111_3112_3113_3114%>% #note that there are a small number of datasets where there is the same number for leaf vs "if leaflet" (so obv it was not a species with compound leaves, but now the data is in there twice); this needs fixing
  #   group_by(DatasetID, species_matched, TraitID2)%>%
  #   summarise(StdValue=mean(StdValue))%>%
  #   spread(key=TraitID2, value=StdValue, fill=NA)%>%
  #   group_by(DatasetID, species_matched, TraitID)%>%
  #   summarise(DatasetValue=mean(StdValue))%>% #averaging by trait and species within each dataset
  #   ungroup()%>%
  #   group_by(species_matched, TraitID)%>%
  #   summarise(SppValue=mean(DatasetValue))%>% #averaging by trait and species across datasets
  #   ungroup()%>%
  #   mutate(TraitID2=paste("Trait",TraitID, sep=''))%>%
  #   select(DatasetID, species_matched, TraitID2, StdValue)
  # 
  # #regressions to check for comparability of different leaf area variables
  # #for whole leaves
  # ggplot(data=subset(traitLeafAreaTest, Trait3108!=Trait3110), aes(x=Trait3108, y=Trait3110)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,2500) + ylim(0,2500) #with vs without petiole for whole leaves
  # cor(traitLeafAreaTest$Trait3108,traitLeafAreaTest$Trait3110, use = "complete.obs") #r=0.86
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3108, y=Trait3112)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,100000) #+ ylim(0,2500) #without vs undefined petiole for whole leaves
  # cor(traitLeafAreaTest$Trait3108,traitLeafAreaTest$Trait3112, use = "complete.obs") #r=0.64
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3108, y=Trait3114)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,100000) #+ ylim(0,2500) #without petiole and whole leaf vs undefined petiole and undefined leaf
  # cor(traitLeafAreaTest$Trait3108,traitLeafAreaTest$Trait3114, use = "complete.obs") #r=0.26
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3110, y=Trait3112)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,100000) #+ ylim(0,2500) #with vs undefined petiole for whole leaves
  # cor(traitLeafAreaTest$Trait3110,traitLeafAreaTest$Trait3112, use = "complete.obs") #r=0.78
  # 
  # #for leaflets
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3109, y=Trait3111)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,50000) + ylim(0,50000) #with vs without petiole for leaflets
  # with(subset(traitLeafAreaTest, Trait3108<50000),cor(Trait3108,Trait3111, use = "complete.obs")) #r=0.49
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3109, y=Trait3113)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,100000) #+ ylim(0,2500) #without vs undefined petiole for leaflets
  # cor(traitLeafAreaTest$Trait3109,traitLeafAreaTest$Trait3113, use = "complete.obs") #r=0.11
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3109, y=Trait3114)) + geom_point() + geom_abline(intercept=0, slope=1) #+ xlim(0,100000) #+ ylim(0,2500) #without petiole and leaflet vs undefined petiole and undefined leaf
  # cor(traitLeafAreaTest$Trait3109,traitLeafAreaTest$Trait3114, use = "complete.obs") #r=0.47
  # 
  # #whole leaves vs leaflets
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3108, y=Trait3109)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,100000) + ylim(0,50000) #whole leaves vs leaflets without petiole
  # with(subset(traitLeafAreaTest, Trait3108<50000),cor(Trait3108,Trait3109, use = "complete.obs")) #r=0.29
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3110, y=Trait3111)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,20000) + ylim(0,20000) #whole leaves vs leaflets with petiole
  # cor(traitLeafAreaTest$Trait3110,traitLeafAreaTest$Trait3111, use = "complete.obs") #r=0.71
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3112, y=Trait3113)) + geom_point() + geom_abline(intercept=0, slope=1) #+ xlim(0,100000) #+ ylim(0,2500) ##whole leaves vs leaflets undefined petiole
  # cor(traitLeafAreaTest$Trait3112,traitLeafAreaTest$Trait3113, use = "complete.obs") #r=0.80

  
  
  
##specific leaf area (SLA) - merging different traits that all correspond to SLA
  #DECISION: combine all SLA data into one clean variable (see below)
  #do everything with the StdValue, which is converted to mm2/mg

#filter outliers
traitSLAGenus <- dat3%>%
  filter(TraitID %in% c(3115, 3116, 3117))%>% #all data related to leaf areas
  separate(species_matched, into=c('genus', 'species'), sep=' ')%>%
  group_by(genus, TraitID)%>%
  summarise(genus_mean=mean(StdValue), genus_sd=sd(StdValue))%>%
  ungroup()

traitSLASpp <- dat3%>%
  filter(TraitID %in% c(3115, 3116, 3117))%>% #all data related to leaf areas
  group_by(species_matched, TraitID)%>%
  summarise(spp_mean=mean(StdValue), spp_sd=sd(StdValue))%>%
  ungroup()

  trait3115_3116_3117_clean<-dat3%>%
    filter(TraitID %in% c(3115, 3116, 3117))%>% #all data related to SLA
    separate(species_matched, into=c('genus', 'species'), sep=' ', remove=F)%>%
    left_join(traitSLAGenus)%>%
    left_join(traitSLASpp)%>%
    mutate(genus_zscore=abs((StdValue-genus_mean)/genus_sd), spp_zscore=abs((StdValue-spp_mean)/spp_sd))%>% #calculate z-scores for genera and species
    filter(genus_zscore<4, spp_zscore<2)%>%
    group_by(DatasetID, species_matched, TraitID)%>%
    summarise(DatasetValue=mean(StdValue))%>% #averaging by ways to measure SLA and species within each dataset
    ungroup()%>%
    group_by(species_matched, TraitID)%>%
    summarise(SppValue=mean(DatasetValue))%>% #averaging by ways to measure SLA and species across datasets
    ungroup()%>%
    group_by(species_matched)%>%
    summarise(CleanTraitValue=mean(SppValue))%>% #averaging by species across datasets and ways to measure SLA
    ungroup()%>%
    mutate(CleanTraitName='SLA', CleanTraitUnit='mm2 mg-1')
  #hist(trait3115_3116_3117$CleanTraitValue)
  
  # #getting averages within each species for each trait type, are they comparable?
  # traitSLAtest <- trait3115_3116_3117%>%
  #   filter(TraitID %in% c(3115, 3116, 3117))%>%
  #   group_by(DatasetID, species_matched, TraitID)%>%
  #   summarise(DatasetValue=mean(StdValue))%>% #averaging by trait and species within each dataset
  #   ungroup()%>%
  #   group_by(species_matched, TraitID)%>%
  #   summarise(SppValue=mean(DatasetValue))%>% #averaging by trait and species across datasets
  #   ungroup()%>%
  #   mutate(TraitID2=paste("Trait",TraitID, sep=''))%>%
  #   select(-TraitID)%>%
  #   spread(key=TraitID2, value=SppValue, fill=NA)
  # 
  # #regressions to check for comparability of different leaf area variables
  # ggplot(data=traitSLAtest, aes(x=Trait3115, y=Trait3116)) + geom_point() + geom_abline(intercept=0, slope=1) #with vs without petiole; close enough, can be combined
  # lm(data=traitSLAtest$Trait3115,traitSLAtest$Trait3116, use = "complete.obs") #r=0.46
  # 
  # ggplot(data=traitSLAtest, aes(x=Trait3115, y=Trait3117)) + geom_point() + geom_abline(intercept=0, slope=1) #with vs undefined petiole; close enough, can be combined
  # cor(traitSLAtest$Trait3115,traitSLAtest$Trait3117, use = "complete.obs") #r=0.73
  # 
  # ggplot(data=traitSLAtest, aes(x=Trait3116, y=Trait3117)) + geom_point() + geom_abline(intercept=0, slope=1) #without vs undefined petiole; close enough, can be combined
  # cor(traitSLAtest$Trait3116,traitSLAtest$Trait3117, use = "complete.obs") #r=0.57

  

##leaf water content - merging different traits that all correspond to leaf water content
#DECISION: combine all leaf water content data into one clean variable (see below)
#do everything with the StdValue, which is converted to mm2
  
  #filter outliers
  traitLeafWaterGenus <- dat3%>%
    filter(TraitID %in% c(3120, 3121, 3122))%>% #all data related to leaf areas
    separate(species_matched, into=c('genus', 'species'), sep=' ')%>%
    group_by(genus, TraitID)%>%
    summarise(genus_mean=mean(StdValue), genus_sd=sd(StdValue))%>%
    ungroup()
  
  traitLeafWaterSpp <- dat3%>%
    filter(TraitID %in% c(3120, 3121, 3122))%>% #all data related to leaf areas
    group_by(species_matched, TraitID)%>%
    summarise(spp_mean=mean(StdValue), spp_sd=sd(StdValue))%>%
    ungroup()
  
  trait3120_3121_3122<-dat3%>%
    filter(TraitID %in% c(3120, 3121, 3122))%>% #all data related to SLA
    separate(species_matched, into=c('genus', 'species'), sep=' ', remove=F)%>%
    left_join(traitLeafWaterGenus)%>%
    left_join(traitLeafWaterSpp)%>%
    mutate(genus_zscore=abs((StdValue-genus_mean)/genus_sd), spp_zscore=abs((StdValue-spp_mean)/spp_sd))%>% #calculate z-scores for genera and species
    filter(genus_zscore<4, spp_zscore<2)
  
  trait3120_3121_3122_clean <- trait3120_3121_3122%>%
    group_by(DatasetID, species_matched, TraitID)%>%
    summarise(DatasetValue=mean(StdValue))%>% #averaging by ways to measure SLA and species within each dataset
    ungroup()%>%
    group_by(species_matched, TraitID)%>%
    summarise(SppValue=mean(DatasetValue))%>% #averaging by ways to measure SLA and species across datasets
    ungroup()%>%
    mutate(AltValue=ifelse(TraitID==3122, (-0.3513 + 1.1673*SppValue), SppValue))%>% #convert saturated to unsaturated (highly correlated)  
    filter(TraitID!=3121)%>% #drop undefined saturation state (not correlated with other two values and least abundant measurement)
    group_by(species_matched)%>%
    summarise(CleanTraitValue=mean(SppValue))%>% #averaging by species across datasets and ways to measure SLA
    ungroup()%>%
    mutate(CleanTraitName='leaf_water_content', CleanTraitUnit='g(W)/g(DM)')
  #hist(trait3115_3116_3117$CleanTraitValue)
  
  # #getting averages within each species for each trait type, are they comparable?
  # traitLeafWatertest <- trait3120_3121_3122%>%
  #   group_by(DatasetID, species_matched, TraitID)%>%
  #   summarise(DatasetValue=mean(StdValue))%>% #averaging by trait and species within each dataset
  #   ungroup()%>%
  #   group_by(species_matched, TraitID)%>%
  #   summarise(SppValue=mean(DatasetValue))%>% #averaging by trait and species across datasets
  #   ungroup()%>%
  #   mutate(TraitID2=paste("Trait",TraitID, sep=''))%>%
  #   select(-TraitID)%>%
  #   spread(key=TraitID2, value=SppValue, fill=NA)
  # 
  # #regressions to check for comparability of different leaf area variables
  # ggplot(data=traitLeafWatertest, aes(x=Trait3120, y=Trait3121)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,5) #undefined vs not saturated
  # summary(lm(data=traitLeafWatertest, Trait3120~Trait3121))
  # 
  # ggplot(data=traitLeafWatertest, aes(x=Trait3120, y=Trait3122)) + geom_point() + geom_abline(intercept=0, slope=1) #saturated vs not saturated
  # summary(lm(data=traitLeafWatertest, Trait3120~Trait3122))
  # 
  # ggplot(data=traitLeafWatertest, aes(x=Trait3121, y=Trait3122)) + geom_point() + geom_abline(intercept=0, slope=1) #undefinted vs saturated
  # summary(lm(data=traitLeafWatertest, Trait3121~Trait3122))
  
  
  
##filtering continuous traits that TRY has already standardized
traitStandardGenus <- dat3%>%
  filter(TraitID %in% c(6,9,12,26,45,46,47,48,55,56,77,80,82,83,84,95,131,145,146,200,363,403,683,1111,2809,3106,3107))%>% #all data related to leaf areas
  separate(species_matched, into=c('genus', 'species'), sep=' ')%>%
  group_by(genus, TraitID)%>%
  summarise(genus_mean=mean(StdValue), genus_sd=sd(StdValue))%>%
  ungroup()

traitStandardSpp <- dat3%>%
  filter(TraitID %in% c(6,9,12,26,45,46,47,48,55,56,77,80,82,83,84,95,131,145,146,200,363,403,683,1111,2809,3106,3107))%>% #all data related to leaf areas
  group_by(species_matched, TraitID)%>%
  summarise(spp_mean=mean(StdValue), spp_sd=sd(StdValue))%>%
  ungroup()

#get list of trait units and types
traitStandardContinuousList <- dat3%>%
  filter(TraitID %in% c(6,9,12,26,45,46,47,48,55,56,77,80,82,83,84,95,131,145,146,200,363,403,683,1111,2809,3106,3107))%>%
  select(TraitID, UnitName)%>%
  unique()%>%
  mutate(remove=ifelse(TraitID==48&UnitName=='', 1, ifelse(TraitID==3107&UnitName=='cm', 1, 0)))%>% #remove two trait unit names that have been fixed (below) or are redundant
  filter(remove==0)%>%
  select(-remove)

traitStandardContinuous_clean <- dat3%>%
  filter(TraitID %in% c(6,9,12,26,45,46,47,48,55,56,77,80,82,83,84,95,131,145,146,200,363,403,683,1111,2809,3106,3107))%>%
  mutate(StdValue2=ifelse(TraitID==3107&UnitName=='cm', StdValue/100, StdValue))%>%  #fix 195 cases where height was cm, but should be standardized to m
  separate(species_matched, into=c('genus', 'species'), sep=' ', remove=F)%>%
  left_join(traitStandardGenus)%>%
  left_join(traitStandardSpp)%>%
  mutate(genus_zscore=abs((StdValue2-genus_mean)/genus_sd), spp_zscore=abs((StdValue2-spp_mean)/spp_sd))%>% #calculate z-scores for genera and species
  filter(genus_zscore<4, spp_zscore<2)%>%
  group_by(DatasetID, species_matched, TraitID)%>%
  summarise(DatasetValue=mean(StdValue2))%>% #averaging by species within each dataset
  ungroup()%>%
  group_by(species_matched, TraitID)%>%
  summarise(CleanTraitValue=mean(DatasetValue))%>% #averaging by across datasets
  ungroup()%>%
  left_join(traitStandardContinuousList)%>%
  mutate(CleanTraitName=ifelse(TraitID==6, 'rooting_depth', ifelse(TraitID==9, 'root:shoot', ifelse(TraitID==12, 'leaf_longevity', ifelse(TraitID==26, 'seed_dry_mass', ifelse(TraitID==45, 'stomata_conductance', ifelse(TraitID==46, 'leaf_thickness', ifelse(TraitID==47, 'LDMC', ifelse(TraitID==48, 'leaf_density', ifelse(TraitID==55, 'leaf_dry_mass', ifelse(TraitID==56, 'leaf_N:P', ifelse(TraitID==77, 'RGR', ifelse(TraitID==80, 'root_N', ifelse(TraitID==82, 'root_density', ifelse(TraitID==83, 'root_diameter', ifelse(TraitID==84, 'root_C', ifelse(TraitID==95, 'germination_efficiency', ifelse(TraitID==131, 'seed_number', ifelse(TraitID==145, 'leaf_width', ifelse(TraitID==146, 'leaf_C:N', ifelse(TraitID==200, 'number_floristic_zones', ifelse(TraitID==363, 'root_dry_mass', ifelse(TraitID==403, 'shoot_dry_mass', ifelse(TraitID==683, 'root_P', ifelse(TraitID==1111, 'seedbank_density', ifelse(TraitID==2809, 'seedbank_duration', ifelse(TraitID==3106, 'plant_height_vegetative', ifelse(TraitID==3107, 'plant_height_regenerative', TraitID))))))))))))))))))))))))))))%>%
  rename(CleanTraitUnit=UnitName)%>%
  #dropping some traits
  filter(CleanTraitName!='germination_efficiency')%>%
  select(-TraitID)

# ggplot(data=traitStandardContinuous_clean, aes(x=CleanTraitValue)) + geom_histogram() + facet_wrap(~CleanTraitName, scales='free')
# #some traits have very skewed distributions, but looking at the data they seem ok (leaf_C:N, leaf_dry_mass, root_dry_mass, seed_dry_mass, seed_number, seedbank_density, seedbank_duration)


##heterotrophy
trait201_clean <- dat3%>%
  filter(TraitID==201)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=='always carnivorous', 'carnivorous', ifelse(OrigValueStr %in% c('always hemiparasitic', 'hemi-parasitic'), 'hemiparasitic', ifelse(OrigValueStr=='mycotrophic', 'mycotrophic', 'autotrophic'))))%>%
  mutate(remove=ifelse(species_matched=='Drosera rotundifolia'&CleanTraitValue=='autotrophic', 1, ifelse(species_matched=='Bartsia alpina'&CleanTraitValue=='autotrophic', 1, ifelse(species_matched=='Botrychium lunaria'&CleanTraitValue=='autotrophic', 1, 0))))%>% #fix three overlapping species (if anything in addition to autotroph, listed as such)
  filter(remove==0)%>%
  mutate(CleanTraitName='heterotrophy', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)%>%
  unique()

##chemical plant defense
trait346 <- dat3%>%
  filter(TraitID==346)%>%
  filter(!is.na(OrigValueStr))%>%
  mutate(CleanTraitValue=OrigValueStr)%>%
  select(species_matched, CleanTraitValue)%>%
  unique()

trait346_clean<-trait346%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Yes=="Yes"&is.na(No), "Yes", 
                         ifelse(No=="No"&is.na(Yes), "No", NA)))%>%
  filter(!is.na(CleanTraitValue))%>%
  select(species_matched, CleanTraitValue)%>%
  mutate(CleanTraitName="allelopathic", CleanTraitUnit=NA)

table(trait346_clean$OrigValueStr)
  
##role of the clonal organ in growth (these are all from one study and categories are not mutually exclusive)
trait357_clean <- dat3%>%
  filter(TraitID==357)%>%
  filter(OrigValueStr!='')%>%
  rename(CleanTraitValue=OrigValueStr)%>%
  mutate(CleanTraitName='clonal_organ_role', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue2=paste(necessary, additive, regenerative, none, sep=','))%>%
  mutate(CleanTraitValue=ifelse(CleanTraitValue2 %in% c('NA,additive,NA,NA', 'NA,additive,NA,none'), 'additive', ifelse(CleanTraitValue2 %in% c('NA,additive,regenerative,NA', 'NA,additive,regenerative,none'), 'additive,regenerative', ifelse(CleanTraitValue2=='NA,NA,NA,none', 'none', ifelse(CleanTraitValue2 %in% c('NA,NA,regenerative,NA', 'NA,NA,regenerative,none'), 'regenerative', ifelse(CleanTraitValue2 %in% c('necessary,additive,NA,NA', 'necessary,additive,NA,none'), 'necessary,additive', ifelse(CleanTraitValue2 %in% c('necessary,additive,regenerative,NA', 'necessary,additive,regenerative,none'), 'necessary,additive,regenerative', ifelse(CleanTraitValue2 %in% c('necessary,NA,NA,NA', 'necessary,NA,NA,none'), 'necessary', 'necessary,regenerative'))))))))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)

##flowering requirement
trait597_clean <- dat3%>%
  filter(TraitID==597)%>%
  filter(!is.na(OrigValueStr))%>%
  rename(CleanTraitValue2=OrigValueStr)%>%
  mutate(CleanTraitName='flowering_requirement', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue2, CleanTraitUnit)%>%
  unique()%>%
  spread(CleanTraitValue2, CleanTraitValue2)%>%
  mutate(CleanTraitValue=ifelse(is.na(High)&is.na(Medium), 'Low', ifelse(is.na(High), 'Medium', 'High')))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)

##vegetative spread rate
trait613_clean <- dat3%>%
  filter(TraitID==613)%>%
  filter(!is.na(OrigValueStr))%>%
  rename(CleanTraitValue2=OrigValueStr)%>%
  mutate(CleanTraitName='vegetative_spread_rate', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue2, CleanTraitUnit)%>%
  unique()%>%
  spread(CleanTraitValue2, CleanTraitValue2)%>%
  mutate(CleanTraitValue=ifelse(is.na(Rapid)&is.na(Moderate)&is.na(Slow), 'None', ifelse(is.na(Rapid)&is.na(Moderate), 'Slow', ifelse(is.na(Rapid), 'Moderate', 'Rapid'))))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)

##stem longevity
trait1187_clean <- dat3%>%
  filter(TraitID==1187)%>%
  rename(CleanTraitValue=OrigValueStr)%>%
  mutate(CleanTraitName='stem_longevity', CleanTraitUnit='year')%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)%>%
  unique()

##stem support
trait1188_clean <- dat3%>%
  filter(TraitID==1188)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=='self', 'self-supporting', OrigValueStr))%>%
  mutate(CleanTraitName='stem_support', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)%>%
  unique()%>%
  mutate(CleanTraitValue2=ifelse(CleanTraitValue=='self-supporting', 'self_supporting', CleanTraitValue))%>%
  select(-CleanTraitValue)%>%
  spread(CleanTraitValue2, CleanTraitValue2, fill='fix')%>%
  mutate(CleanTraitValue2=paste(creeping, decumbent, procumbent, scrambling, self_supporting, tendrils, twining, sep=','))%>%
  mutate(CleanTraitValue=ifelse(CleanTraitValue2 %in% c('creeping,fix,fix,fix,fix,fix,fix','creeping,decumbent,fix,fix,fix,fix,fix'), 'creeping', ifelse(CleanTraitValue2 %in% c('fix,decumbent,fix,fix,fix,fix,fix', 'fix,decumbent,fix,scrambling,fix,fix,fix'), 'decumbent', ifelse(CleanTraitValue2 %in% c('fix,decumbent,procumbent,fix,fix,fix,fix', 'fix,fix,procumbent,fix,fix,fix,fix'), 'procumbent', ifelse(CleanTraitValue2=='fix,fix,fix,fix,fix,fix,twining', 'twining', ifelse(CleanTraitValue2=='fix,fix,fix,fix,fix,tendrils,fix', 'tendrils', ifelse(CleanTraitValue2=='fix,fix,fix,scrambling,fix,fix,fix', 'scrambling', 'self-supporting')))))))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)


#N-fixation capability
trait8_clean <- dat3%>%
  filter(TraitID==8)%>%
  mutate(remove=ifelse(DatasetID==444&OriglName=='Nitrogen fixation capacity', 1, ifelse(DatasetID==444&OriglName=='NFC', 1, 0)))%>% #for these categories, this study lists everything as a yes, which is not true. other categories in this study are more informative
  filter(!is.na(OrigValueStr), remove==0)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr %in% c('0', 'n', 'N', 'no', 'No', 'NO-N-fixer', 'not N2 fixing', 'no, not an N fixer', 'low'), 'no', ifelse(DatasetID==73&OrigValueStr==1, 'no', 'yes')))%>%
  group_by(species_matched, CleanTraitValue)%>%
  summarise(count=length(CleanTraitValue))%>%
  spread(key=CleanTraitValue, value=count, fill=0)%>%
  mutate(CleanTraitValue2=ifelse(no>yes, 'no', 'yes'), check=ifelse(no>0&yes>0, 'check', '0'))%>%
  mutate(CleanTraitValue=ifelse(species_matched %in% c('Dryas integrifolia', 'Senna marilandica'), 'yes', CleanTraitValue2))%>%
  mutate(CleanTraitName='N_fixation', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)


##leaf palatability  --  figure out what categories mean
trait152_clean <- dat3%>%
  filter(TraitID==152)%>%
  mutate(CleanTraitName='leaf_palatability', CleanTraitUnit=NA)%>% 
  rename(CleanTraitValue=OrigValueStr)%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)


#palatability
trait679 <- dat3%>% #actually 5 traits: bloat, toxicity, palatable to graze animals, palatable to browse animals, palatable to humans (then one more study with just "palatable", which I'm adding to the graze animal category based on Leichman study)
  filter(TraitID==679)%>%
  filter(!is.na(OrigValueStr))%>%
  mutate(CleanTraitName=ifelse(OriglName=='Bloat', 'palatability_bloat', ifelse(OriglName=='Palatable Browse Animal', 'palatability_browse', ifelse(OriglName %in% c('Palatable Graze Animal', 'PA, palatability'), 'palatability_graze', ifelse(OriglName=='Palatable Human', 'palatability_human', 'toxicity')))), CleanTraitUnit=NA)%>% 
  rename(CleanTraitValue=OrigValueStr)%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)%>%
  unique()

#need to take out each variable and then make sure the terms are all fine, then bind them all back together in the end
trait679_bloat <- trait679%>%
  filter(CleanTraitName=='palatability_bloat')%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(is.na(High)&is.na(Medium)&is.na(Low), 'None', ifelse(is.na(High)&is.na(Medium), 'Low', ifelse(is.na(High), 'Medium', 'High'))))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)

trait679_browse <- trait679%>%
  filter(CleanTraitName=='palatability_browse')%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(is.na(High)&is.na(Medium), 'Low', ifelse(is.na(High), 'Medium', 'High')))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)

trait679_graze <- trait679%>%
  filter(CleanTraitName=='palatability_graze')%>%
  mutate(CleanTraitValue2=ifelse(CleanTraitValue %in% c('high', 'High'), 'High', ifelse(CleanTraitValue %in% c('Medium', 'moderate', 'variable (e.g. young plants palatable but adult plant not)'), 'Medium', 'Low')))%>%
  select(-CleanTraitValue)%>%
  spread(CleanTraitValue2, CleanTraitValue2)%>%
  mutate(CleanTraitValue=ifelse(is.na(High)&is.na(Medium), 'Low', ifelse(is.na(High), 'Medium', 'High')))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)

trait679_human <- trait679%>%
  filter(CleanTraitName=='palatability_human')%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(is.na(No), 'Yes', 'No'))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)

trait679_toxicity <- trait679%>%
  filter(CleanTraitName=='toxicity')%>%
  mutate(CleanTraitValue2=ifelse(CleanTraitValue %in% c('high', 'Severe'), 'High', ifelse(CleanTraitValue %in% c('medium', 'Moderate'), 'Medium', ifelse(CleanTraitValue %in% c('low', 'Slight'), 'Low', 'None'))))%>%
  select(species_matched, CleanTraitName, CleanTraitValue2, CleanTraitUnit)%>%
  unique()%>%
  spread(CleanTraitValue2, CleanTraitValue2)%>%
  mutate(CleanTraitValue=ifelse(is.na(High)&is.na(Medium)&is.na(Low), 'None', ifelse(is.na(High)&is.na(Medium), 'Low', ifelse(is.na(High), 'Medium', 'High'))))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)

trait679_clean <- rbind(trait679_bloat, trait679_browse, trait679_graze, trait679_human, trait679_toxicity)


#physical defenses
trait345 <- dat3%>% #all comes from one dataset; actually 3 different traits (defenses of stem, leaves, flower/fruit)
  filter(TraitID==345)%>%
  filter(!is.na(OrigValueStr))%>%
  mutate(CleanTraitName=ifelse(OriglName=='Physical defences on flowers/fruits', 'defenses_fruit/flower', ifelse(OriglName=='Physical defences on leaves', 'defenses_leaves', 'defenses_stem')), CleanTraitUnit=NA)%>% 
  mutate(CleanTraitValue=ifelse(OrigValueStr=='Scales', 'scales', ifelse(OrigValueStr %in% c('glandular hairs', 'hairy', 'dense hairs', 'soft hairs', 'stiff hairs'), 'hairs', ifelse(OrigValueStr=='stinging hairs', 'stinging', ifelse(OrigValueStr %in% c('spiny point', 'spines'), 'spines', ifelse(OrigValueStr=='thick cuticle', 'thick_cuticle', OrigValueStr))))))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)

trait345_hair <- trait345%>%
  #making two different traits: (1) hairs or glabrous and (2) other defenses
  mutate(hair=paste(glabrous, subglabrous, hairs, sep=','), CleanTraitName=ifelse(CleanTraitName=='defenses_fruit/flower', 'hairs_fruit/flower', ifelse(CleanTraitName=='defenses_leaves', 'hairs_leaves', 'hairs_stem')))%>%
  mutate(CleanTraitValue=ifelse(hair=='glabrous,NA,NA', 'glabrous', ifelse(hair=='NA,subglabrous,NA', 'subglabrous', 'hairs')))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)

trait345_other <- trait345%>%
  #making two different traits: (1) hairs or glabrous and (2) other defenses
  mutate(other=paste(scales, stinging, spines, thick_cuticle, glandular, prickles, thorns, viscid, sep=','))%>%
  mutate(CleanTraitValue=ifelse(other=='NA,NA,NA,NA,glandular,NA,NA,NA', 'glandular', ifelse(other=='NA,NA,NA,NA,NA,NA,NA,viscid', 'viscid', ifelse(other=='NA,NA,NA,NA,NA,NA,thorns,NA', 'thorns', ifelse(other=='NA,NA,NA,NA,NA,prickles,NA,NA', 'prickles', ifelse(other=='NA,NA,NA,thick_cuticle,NA,NA,NA,NA', 'thick_cuticle', ifelse(other=='NA,NA,spines,NA,NA,NA,NA,NA', 'spines', ifelse(other=='NA,stinging,NA,NA,NA,NA,NA,NA', 'stinging', ifelse(other=='scales,NA,NA,NA,NA,NA,NA,NA', 'scales', 'drop')))))))))%>%
  filter(CleanTraitValue!='drop')%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)



#clonal growth
#341 has to be clonal is it it how it is clonal
trait341<-dat3%>%
  filter(TraitID==341)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=="no", "NotClonal", "Clonal"))%>%
  select(species_matched, CleanTraitValue)%>%
  filter(!is.na(CleanTraitValue))%>%
  unique()

table(trait341$OrigValueStr)

#329 is very similar to 341, but the low categories too
trait329<-dat3%>%
  filter(TraitID==329)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=="Absence", "NotClonal",
                                       ifelse(OrigValueStr=="Little or no vegetative spread", "LowClonality", "Clonal")))%>%
  select(species_matched, CleanTraitValue)%>%
  filter(!is.na(CleanTraitValue))%>%
  unique()

table(trait329$CleanTraitValue)

#344 is very similar to 329/341, but the low categories too
trait334<-dat3%>%
  filter(TraitID==334)%>%
  mutate(CleanTraitValue=ifelse(OrigUnitStr=="cm"|OrigUnitStr==""|OrigValueStr=="dispersable"|OrigValueStr=="", NA, 
                         ifelse(OrigValueStr=="<0.01", "low", 
                         ifelse(OrigValueStr==">0.25", "high", "mid"))))%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  filter(!is.na(CleanTraitValue))

trait334_clean<-trait334%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(low=="low"&is.na(high)&is.na(mid), "<0.01",
                         ifelse(mid=="mid"&is.na(high)&is.na(low)|mid=="mid"&is.na(high)&low=="low","	0.01-0.25", ">0.25" )))%>%
  select(species_matched, CleanTraitValue)%>%
  mutate(CleanTraitName="rate_veg_spread", CleanTraitUnit="m/yr")
  
  
  
table(trait334$CleanTraitValue)

trait344<-dat3%>%
  filter(TraitID==344)%>%
  mutate(CleanTraitValue=ifelse(DatasetID==92|OrigValueStr=="no particuliar mating system", NA,
         ifelse(OrigValueStr=="none"|OrigValueStr=="NotClonal", "NotClonal",
        ifelse(OrigValueStr=="Little or no vegetative spread", "LowClonality", "Clonal"))))%>%
select(species_matched, CleanTraitValue)%>%
  filter(!is.na(CleanTraitValue))%>%
  unique()

table(trait344$OriglName)


trait_341_329_344_clonal<-trait341%>%
  full_join(trait329)%>%
  full_join(trait344)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue) %>% 
  mutate(CleanTraitValue=ifelse(LowClonality=="LowClonality"&is.na(NotClonal)&is.na(Clonal)|LowClonality=="LowClonality"&NotClonal=="NotClonal"&is.na(Clonal), "LowClonality", ifelse(Clonal=="Clonal"&NotClonal=="NotClonal"&is.na(LowClonality)|Clonal=="Clonal"&LowClonality=="LowClonality"&is.na(NotClonal)|Clonal=="Clonal"&is.na(LowClonality)&is.na(NotClonal)|Clonal=="Clonal"&NotClonal=="NotClonal"&LowClonality=="LowClonality", "Clonal", NA)))%>%
  filter(!is.na(CleanTraitValue))%>%
  select(species_matched, CleanTraitValue)

trait_341_329_344_notclonal<-trait341%>%
  full_join(trait329)%>%
  full_join(trait344)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue) %>%
  mutate(CleanTraitValue=ifelse(NotClonal=="NotClonal"&is.na(LowClonality)&is.na(Clonal),"NotClonal", NA))%>%
  filter(!is.na(CleanTraitValue))%>%
  select(species_matched, CleanTraitValue)

trait_341_329_344_clean<-trait_341_329_344_clonal%>%
  bind_rows(trait_341_329_344_notclonal)%>%
  mutate(CleanTraitName="clonality", CleanTraitUnit=NA)


#609 is what I want, how can they reproduce, use 344/341/329 to check i have all clonals 609 and 208 are the same and need to be combined.
trait609<-dat3%>%
  filter(TraitID==609)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=="No"|OrigValueStr=="Yes", NA, ifelse(OrigValueStr=="seed and vegetative", "seed_and_vegetative", OrigValueStr)))%>%
  filter(!is.na(CleanTraitValue))%>%
  select(species_matched, CleanTraitValue)%>%
  unique()

trait208<-dat3%>%
  filter(TraitID==208)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=="generative"|OrigValueStr=="vegetative & generative", NA, 
                                ifelse(OrigValueStr=="vegetative"|OrigValueStr=="vegetatively ", "vegetative",
                                       ifelse(OrigValueStr=="by seed/by spore", "seed", "seed_and_vegetative"))))%>%
  filter(!is.na(CleanTraitValue))%>%
  select(species_matched, CleanTraitValue)%>%
  unique()

trait609_208_clean<-trait609%>%
  bind_rows(trait208)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(seed=="seed"&is.na(vegetative)&is.na(seed_and_vegetative), "seed", ifelse(vegetative=="vegetative"&is.na(seed)&is.na(seed_and_vegetative),"vegetative", "seed_and_vegetative")))%>%
  select(species_matched, CleanTraitValue) %>% 
  mutate(CleanTraitName="reproduction_seed_veg", CleanTraitUnit=NA)

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
  unique()

trait28_clean<-trait28%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Animal=="Animal"&is.na(Wind)&is.na(Unassisted)&is.na(Water), "Animal",
                         ifelse(Wind=="Wind"&is.na(Animal)&is.na(Unassisted)&is.na(Water),"Wind",
                         ifelse(Unassisted=="Unassisted"&is.na(Animal)&is.na(Wind)&is.na(Water), "Unassisted",
                         ifelse(Water=="Water"&is.na(Animal)&is.na(Wind)&is.na(Unassisted), "Water", "Combination")))))%>%
  select(species_matched, CleanTraitValue)%>%
  mutate(CleanTraitName="dispersal_mode", CleanTraitUnit="NA")


table(trait28_clean$CleanTraitValue)

# #woodiness - not very usable. Already have with plant growth form
# trait38<-dat3%>%
#   filter(TraitID==38)
# 
# table(trait38$OrigValueStr)



#Grime triats - not very usable. Already have with plant growth form - will take more thinking.
# trait196<-dat3%>%
#   filter(TraitID==196)
# 
# table(trait196$OrigValueStr)


#combining traits

traits_cont <- trait3115_3116_3117_clean%>%
  bind_rows(trait3108_3109_3110_3111_3112_3113_3114_clean, trait3120_3121_3122_clean, traitStandardContinuous_clean, trait_1433_1030_clean)%>%
  mutate(TraitCategory=ifelse(CleanTraitName=="rooting_depth"|CleanTraitName=="root:shoot"|CleanTraitName=="root_dry_mass"|CleanTraitName=="plant_height_vegetative"|CleanTraitName=="plant_height_regenerative"|CleanTraitName=="RGR"|CleanTraitName=="shoot_dry_mass", "Growth",
                       ifelse(CleanTraitName=="number_floristic_zones", "Habitat",
                       ifelse(CleanTraitName=="mycorrhizal_percent_colonization", "Mutualism",
                       ifelse(CleanTraitName=="leaf_N:P"|CleanTraitName=="leaf_C:N"|CleanTraitName=="root_N"|CleanTraitName=="root_P"|CleanTraitName=="root_C", "Nutrients",
                       ifelse(CleanTraitName=="stomata_conductance", "Physiology", 
                       ifelse(CleanTraitName=="seed_number"|CleanTraitName=="seed_dry_mass"|CleanTraitName=="seedbank_density"|CleanTraitName=="seedbank_duration", "Reproduciton",
                       ifelse(CleanTraitName=="root_diameter"|CleanTraitName=="root_density", "RES", "LES"))))))))%>%
  mutate(TraitType="Continuous")

write.csv(traits_cont, "C://Users/mavolio2/Dropbox/SDiv_sCoRRE_shared/CoRRE - community and anpp data/TRY_trait_data_continuous.csv", row.names = F)


#removed b/c not enough representation trait345_clean,
traits_cat<-trait28_clean%>%
  bind_rows(trait28_clean, trait_341_329_344_clean, trait609_208_clean, trait59_clean, trait201_clean, trait346_clean, trait357_clean, trait597_clean, trait613_clean, trait1187_clean, trait1188_clean, trait679_clean, trait8_clean, trait152_clean, trait42_clean, trait22_clean,  trait7_clean, trait334_clean)%>%
  mutate(TraitCategory=ifelse(CleanTraitName=="palatability_bloat"|CleanTraitName=="palatability_browse"|CleanTraitName=="palatability_graze"|CleanTraitName=="palatability_human"|CleanTraitName=="toxicity"|CleanTraitName=="leaf_palatability"|CleanTraitName=="stem_longevity", "Herbivory",
                       ifelse(CleanTraitName=="lifeform"|CleanTraitName=="lifespan"|CleanTraitName=="heterotrophy"|CleanTraitName=="allelopathic"|CleanTraitName=="stem_support", "Life_History",
                       ifelse(CleanTraitName=="mycorrhizal"|CleanTraitName=="N_fixation", "Mutualism",
                       ifelse(CleanTraitName=="photo_pathway", "Physiology", "Reproduction")))))%>%
  mutate(TraitType="Categorical")


write.csv(traits_cat, "C://Users/mavolio2/Dropbox/SDiv_sCoRRE_shared/CoRRE - community and anpp data/TRY_trait_data_categorical.csv", row.names = F)

summary_cat<-traits_cat%>%
  group_by(CleanTraitName, TraitCategory, TraitType)%>%
  summarize(n=length(CleanTraitValue))%>%
  mutate(per.sp=(n/1954)*100)

summary_cont<-traits_cont%>%
  group_by(CleanTraitName, TraitCategory, TraitType)%>%
  summarize(n=length(CleanTraitValue))%>%
  mutate(per.sp=(n/1954)*100)

summary_all<-summary_cat%>%
  bind_rows(summary_cont)

ggplot(data=summary_all, aes(x=reorder(CleanTraitName, -per.sp), y=per.sp))+
  geom_point()+
  ylab("percent of species")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank())+
  xlab("")+
  scale_y_continuous(limits=c(0,100))

#barplots by trait type
ggplot(data=subset(summary_all, TraitCategory=="Reproduction"), aes(x=1, y=per.sp))+
  geom_bar(stat = "identity")+
 facet_wrap(~CleanTraitName)+
  scale_y_continuous(limits=c(0,100))+
  ylab("Percent of Species")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank())+
  xlab("")

