################################################################################
##  sCoRRE_trait_ordination.R: Figuring out what traits to use.
##
##  Author: Kimberly Komatsu
##  Date created: March 10, 2021
################################################################################

library(ecodist)
library(FD)
library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared')

#read in data
contTraits <- read.csv('Trait Data\\TRY Data\\Gap_Filled\\TRY_new.csv')%>%
  rename(species_matched=Species)%>%
  select(-X.1, -X, -Family, -Genus, -ObservationID)%>%
  group_by(species_matched)%>%
  summarise_all(funs(mean))%>%
  ungroup()

traits <- read.csv('CoRRE data\\CoRRE data\\trait data\\sCoRRE categorical trait data - traits_complete_pre spot check_03102021.csv')%>%
  full_join(contTraits)%>%
  drop_na()

traitMatrix <- distance(traits[,15:47], method='gower') #ignoring all categorical traits

traitPCO <- pco(traitMatrix)



plot(traitPCO$vectors[,1:2], col=as.numeric(traits$species_matched),
     pch=as.numeric(traits$species_matched), main="PCO", xlab="PCO 1", ylab="PCO 2")

PCO <- data.frame(traitPCO$vectors[,1:2], species_matched=traits$species_matched)

ggplot(data=PCO, aes(x=X1, y=X2)) +
  geom_point()
