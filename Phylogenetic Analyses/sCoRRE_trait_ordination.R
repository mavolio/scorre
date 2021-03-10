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
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\")

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

# Create Gower trait disimilarity matrix
traitMatrix <- distance(traits[,15:47], method='gower') #ignoring all categorical traits

# Run PCoA
traitPCO <- pco(traitMatrix)

# Create matrix of first two axes from PCoA
PCOOutMatrix <- traitPCO$vectors[,1:2]

# Create vector fitting object (for overlaying trait vectors on top of PCoA plot) -- should probably up the permutations to 1000
trait_vf <- vf(PCOOutMatrix, traits[,15:47], nperm=100)

### Plotting
# plot(PCOOutMatrix, col=1:length(traits$species_matched),
#       pch=1:length(traits$species_matched), main="PCO", xlab="PCO 1", ylab="PCO 2")
plot(PCOOutMatrix, main="PCO", xlab="PCO 1", ylab="PCO 2")
plot(trait_vf)



PCO <- data.frame(traitPCO$vectors[,1:2], species_matched=traits$species_matched)
names(traitPCO)
traitPCO$values

ggplot(data=PCO, aes(x=X1, y=X2, label=species_matched)) +
#  geom_point() +
  geom_text()
