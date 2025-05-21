################################################################################
##  sCoRRE_dCCA_traits.R: Examining trait spaces occupied by species in various treatments.
##
##  Author: Kimberly Komatsu
##  Date created: August 12, 2022
################################################################################

library(ade4)
library(vegan)
library(FD)
library(psych)
library(FactoMineR)
library(devtools)
# install_github('fawda123/ggord')
library(ggord)
library(ggfortify)
library(tidyverse)

setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\paper 2_PD and FD responses\\data\\')  #kim's laptop

##### functions and themes #####
###standard error function
se <- function(x, na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  

#theme set
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

##### data #####
sppNames <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data\\trait data\\corre2trykey_2021.csv')%>%
  select(genus_species, species_matched)%>%
  unique()

trt <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\RawAbundance.csv')%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id)%>%
  unique()%>%
  left_join(read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CompiledData\\ExperimentInfo.csv'))%>%
  group_by(site_code, project_name, community_type)%>%
  mutate(experiment_length=max(treatment_year))%>%
  ungroup()%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, trt_type, experiment_length, plot_mani, n, p, CO2, precip, temp)%>%
  mutate(alltrts=ifelse(trt_type %in% c("control", "CO2","CO2*temp", "mow_clip","burn","burn*graze","disturbance","burn*mow_clip","drought","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N","mult_nutrient","N*P","P","N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp","temp*mow_clip","drought*temp","irr*temp","irr"),1,0))%>%
  filter(alltrts==1)%>%
  mutate(dist=ifelse(trt_type %in% c("mow_clip","burn","burn*graze","disturbance","burn*mow_clip"), 1, 0),
         # tCO2=ifelse(trt_type %in% c("CO2"), 1, 0),
         drought=ifelse(trt_type %in% c("drought"), 1, 0),
         # therb_removal=ifelse(trt_type %in% c("herb_removal"), 1, 0),
         irg=ifelse(trt_type %in% c("irr"), 1, 0),
         # ttemp=ifelse(trt_type %in% c("temp"), 1, 0),
         # tn=ifelse(trt_type %in% c("N"), 1, 0),
         # tp=ifelse(trt_type %in% c("P"), 1, 0),
         multtrts=ifelse(trt_type %in% c("CO2*temp", "burn*graze","burn*mow_clip","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp*mow_clip","drought*temp","irr*temp","mult_nutrient","N*P"),1,0))%>%
  mutate(trt_type2=ifelse(dist==1, 'disturbance', ifelse(multtrts==1, 'multiple trts', trt_type)))%>%
  mutate(rep=paste(site_code, project_name, community_type, treatment, sep='::'))%>%
  select(-site_code, -project_name, -community_type, -treatment)%>%
  select(rep, n, p, CO2, precip, temp, dist, multtrts, trt_type2)%>%
  unique()

contTraits <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data\\trait data\\Final TRY Traits\\Imputed Continuous_Traits\\data to play with\\imputed_continuous_20220620.csv')%>%
  select(-X.1, -X, -family, -genus, -observation)%>%
  group_by(species_matched)%>%
  summarise_all(funs(mean))%>%
  ungroup()

traits_all <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data\\trait data\\sCoRRE categorical trait data_11302021.csv')%>%
  full_join(contTraits) %>%
  drop_na()

traitsOutliersRemoved <- traits_all %>%
  filter(!leaf_type %in% c("microphyll","frond")) %>%
  filter(!species_matched %in% c("Centrolepis aristata", "Centrolepis strigosa", "Acorus calamus"))

traitsScaled <- traitsOutliersRemoved %>% ## only scales continuous traits
  mutate_at(vars(seed_dry_mass:seed_number), scale)%>%
  select(species_matched, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal, n_fixation, seed_dry_mass, leaf_N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_number)%>%
  rename(genus_species=species_matched)%>%
  arrange(genus_species)%>%
  mutate(keep=ifelse(growth_form=='CHECK' | growth_form=='' | photosynthetic_pathway=='CHECK' | photosynthetic_pathway=='' | lifespan=='CHECK' | lifespan=='' | clonal=='CHECK' | clonal=='' | mycorrhizal=='CHECK' | mycorrhizal=='' | n_fixation=='CHECK' | n_fixation=='', 0, 1))%>%
  filter(keep==1)%>%
  select(-keep)

spp <- read.csv('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data\\CoRRE data\\community composition\\CoRRE_RelativeCover_Dec2021.csv')%>%
  left_join(sppNames)%>%
  filter(!is.na(species_matched))%>%
  select(-genus_species)%>%
  rename(genus_species=species_matched)%>%
  group_by(site_code, project_name, community_type, treatment, genus_species)%>%
  summarise(cover=mean(relcov))%>%
  ungroup()%>%
  mutate(rep=paste(site_code, project_name, community_type, treatment, sep='::'))%>%
  select(-site_code, -project_name, -community_type, -treatment)%>%
  mutate(keep=ifelse(genus_species %in% c(traitsScaled$genus_species), 1, 0))%>%
  filter(keep==1)%>%
  select(-keep)

traitsFinal <- traitsScaled%>%
  mutate(keep=ifelse(genus_species %in% c(spp$genus_species), 1, 0))%>%
  filter(keep==1)%>%
  select(-keep)

sppTraits <- spp%>%
  left_join(traitsFinal)%>%
  mutate(growth_form_ordinal=ifelse(growth_form=='graminoid', 1,
                                    ifelse(growth_form=='forb', 2,
                                           ifelse(growth_form=='vine', 3, 4))),
         photosynthetic_pathway_ordinal=ifelse(photosynthetic_pathway=='C3', 1, 2),
         lifespan_ordinal=ifelse(lifespan=='perennial', 2, 1),
         clonal_ordinal=ifelse(clonal=='yes', 2, 1),
         mycorrhizal_ordinal=ifelse(mycorrhizal=='yes', 2, 1),
         n_fixation_ordinal=ifelse(n_fixation=='yes', 2, 1))

##### case study examples #####

### nitrogen
sppTraitsKUFS <- sppTraits%>%
  filter(rep %in% c('KUFS::E2::0::N1S0H0', 'KUFS::E2::0::N0S0H0'))

PCAmodel <- prcomp(data.matrix(sppTraitsKUFS[,10:22]))
autoplot(PCAmodel, data=sppTraitsKUFS, colour='rep', size='cover', 
         loadings=TRUE, loadings.colour='dark grey',
         loadings.label=F) +
  scale_color_manual(values=c('#595959', '#0000FF')) +
  coord_cartesian(xlim=c(-0.35,0.2)) +
  theme(legend.position='none')

### irrigation
sppTraitsDL <- sppTraits%>%
  filter(rep %in% c('DL::GCME2::0::P', 'DL::GCME2::0::C'))

PCAmodel <- prcomp(data.matrix(sppTraitsDL[,10:22]))
autoplot(PCAmodel, data=sppTraitsDL, colour='rep', size='cover', 
         loadings=T, loadings.colour='dark grey',
         loadings.label=F) +
  scale_color_manual(values=c('#595959', '#0000FF')) + 
  theme(legend.position='none')

### drought
sppTraitsKAEFS <- sppTraits%>%
  filter(rep %in% c('KAEFS::WAPAClip::0::U CC', 'KAEFS::WAPAClip::0::U CH'))

PCAmodel <- prcomp(data.matrix(sppTraitsKAEFS[,10:22]))
autoplot(PCAmodel, data=sppTraitsKAEFS, colour='rep', size='cover', 
         loadings=T, loadings.colour='dark grey',
         loadings.label=T) +
  scale_color_manual(values=c('#FFA300', '#595959')) + 
  theme(legend.position='none')



##### community weighted means #####

communityWeightedMeans <- spp%>%
  left_join(traitsFinal)%>%
  mutate(growth_form_ordinal=ifelse(growth_form=='graminoid', 1,
                                    ifelse(growth_form=='forb', 2,
                                           ifelse(growth_form=='vine', 3, 4))),
         photosynthetic_pathway_ordinal=ifelse(photosynthetic_pathway=='C3', 1, 2),
         lifespan_ordinal=ifelse(lifespan=='perennial', 2, 1),
         clonal_ordinal=ifelse(clonal=='yes', 2, 1),
         mycorrhizal_ordinal=ifelse(mycorrhizal=='yes', 2, 1),
         n_fixation_ordinal=ifelse(n_fixation=='yes', 2, 1))%>%
  mutate(growth_form_weighted=cover*growth_form_ordinal,
         photosynthetic_pathway_weighted=cover*photosynthetic_pathway_ordinal,
         lifespan_weighted=cover*lifespan_ordinal,
         clonal_weighted=cover*clonal_ordinal,
         mycorrhizal_weighted=cover*mycorrhizal_ordinal,
         n_fixation_weighted=cover*n_fixation_ordinal,
         seed_dry_mass_weighted=cover*seed_dry_mass,
         leaf_N_weighted=cover*leaf_N,
         LDMC_weighted=cover*LDMC,
         SLA_weighted=cover*SLA,
         plant_height_vegetative_weighted=cover*plant_height_vegetative,
         rooting_depth_weighted=cover*rooting_depth,
         seed_number_weighted=cover*seed_number)%>%
  group_by(rep)%>%
  summarise_at(vars(growth_form_weighted, photosynthetic_pathway_weighted, lifespan_weighted, clonal_weighted, mycorrhizal_weighted, n_fixation_weighted, seed_dry_mass_weighted, leaf_N_weighted, LDMC_weighted, SLA_weighted, plant_height_vegetative_weighted, rooting_depth_weighted, seed_number_weighted), mean)%>%
  ungroup()%>%
  left_join(trt)%>%
  na.omit()

##### full PCA #####
PCAmodel <- prcomp(data.matrix(communityWeightedMeans[,2:14]))
autoplot(PCAmodel, data=communityWeightedMeans, colour='trt_type2',
         loadings=T, loadings.color='green',
         loadings.label=T, loadings.label.size=3)

### N only
Nexpts <- communityWeightedMeans%>%
  filter(trt_type2=='N')%>%
  separate(rep, into=c('site_code', 'project_name', 'community_type', 'treatment'), sep='::')%>%
  select(site_code, project_name, community_type)%>%
  unique()

communityWeightedMeansN <- communityWeightedMeans%>%
  separate(rep, into=c('site_code', 'project_name', 'community_type', 'treatment'), sep='::')%>%
  filter(site_code %in% c(Nexpts$site_code) & project_name %in% c(Nexpts$project_name) & community_type %in% c(Nexpts$community_type))%>%
  filter(trt_type2 %in% c('control', 'N'))%>%
  mutate(N_factor=ifelse(n>0, 'N', 'control'))
  
PCAmodel <- prcomp(data.matrix(communityWeightedMeansN[,5:17]))
autoplot(PCAmodel, data=communityWeightedMeansN, colour='N_factor',
         loadings=T, loadings.color='green',
         loadings.label=T, loadings.label.size=3)





##### things that don't work #####

# RDAmodel <- rda(data.matrix(sppTraits[,2:14]) ~ n + p + CO2 + precip + temp + dist + multtrts, data = sppTraits)
# RDAmodel
# summary(RDAmodel)
# 
# ggord(RDAmodel, sppTraits$rep, poly = FALSE, ptslab = TRUE, size = 3,veclsz = 0.4, arrow = 0.3, addcol = "grey10", grp_title = "Site", repel = TRUE, alpha = 2) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   scale_linetype_manual(values = c('solid', 'solid'))+
#   labs(color="", linetype = "", shape = "")+
#   scale_colour_manual(values=c("grey60", "grey80"), labels = c("BLM", "Sun Prairie"), name = "")+
#   scale_linetype_manual(values = c("twodash", "solid"), labels = c("BLM", "Sun Prairie"), name = "")+
#   scale_shape_manual(values = c(17, 16),labels = c("BLM", "Sun Prairie"), name = "")
#coord_cartesian(xlim=c(-1.1,1.2), ylim=c(-1.1,1.2))

# pca.cwm <- dudi.pca(sppTraits[,2:14], scannf=FALSE)
# summary(pca.cwm)
# rda.cwm <- pcaiv(pca.cwm, trt[2:8], scannf=FALSE)

# FAMDmodel <- FAMD(sppTraits[,2:14])
# summary(FAMDmodel)
# res <- Factoshiny(sppTraits[,2:14])

# sppMatrix <- spp%>%
#   spread(key=genus_species, value=cover, fill=0)%>%
#   select(-rep)
# 
# test <- sppMatrix%>%
#   gather(key='genus_species_1', value='cover')%>%
#   select(-cover)%>%
#   unique()%>%
#   cbind(traitsFinal)%>%
#   mutate(same=ifelse(genus_species==genus_species_1, 1, 0))
# 
# ##### double CCA #####
# ca1 <- dudi.coa(spp, scannf = F)
# 
# dCCA1 <- dbrda(ca1, trt, traitsScaled, scannf = FALSE)
# 
# 
# traitsTest <- traitsFinal[1:10,]%>%
#   mutate(name=LETTERS[seq(from=1, to=10)])%>%
#   select(genus_species, name, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal, n_fixation, seed_dry_mass, leaf_N, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_number)
# 
# sppTest <- spp%>%
#   filter(genus_species %in% c(traitsTest$genus_species))%>%
#   left_join(traitsTest)%>%
#   select(rep, cover, name)%>%
#   spread(key=name, value=cover, fill=0)%>%
#   select(-rep)
# 
# traitsTest2 <- traitsTest%>%
#   select(-genus_species)
#   
# ##### community weighted means #####
# functcomp(traitsTest2, data.matrix(sppTest), CWM.type = "dom", bin.num = NULL)
