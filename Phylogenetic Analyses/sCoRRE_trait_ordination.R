################################################################################
##  sCoRRE_trait_ordination.R: Figuring out what traits to use.
##
##  Author: Kimberly Komatsu
##  Date created: March 10, 2021
################################################################################

library(ecodist)
library(FD)
library(PerformanceAnalytics)
library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared')
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\")
setwd('C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared')

#our traits to determine how much data we are missing

dat<-read.csv("Trait Data\\TRY Data\\Try Continuous data\\TRY_trait_data_continuous_long.csv")

length(unique(dat$species_matched))

numsp<-dat%>%
  group_by(species_matched, CleanTraitName)%>%
  summarize(ave=mean(StdValue))%>%
  group_by(CleanTraitName)%>%
  summarize(n=length(ave))%>%
  mutate(pct=n/2400)
  
play<-dat%>%
  group_by(species_matched, CleanTraitName)%>%
  summarize(ave=mean(StdValue))%>%
  spread(CleanTraitName, ave, fill=0)%>%
  filter(seed_number<3e+08)

with(subset(play, SRL>00&root_density>0), plot(SRL, root_density))
     
x<-play$seed_number

#read in data
contTraits <- read.csv('Trait Data\\TRY Data\\Gap_Filled\\TRY_new.csv')%>%
  rename(species_matched=Species)%>%
  select(-X.1, -X, -Family, -Genus, -ObservationID)%>%
  group_by(species_matched)%>%
  summarise_all(funs(mean))%>%
  ungroup()

contTraitsSubset <- contTraits%>%
  rename(ssd=X4, rooting_depth=X6, SLA=X11, leaf_C_mass=X13, leaf_N_mass=X14, leaf_P_mass=X15, stem_diameter=X21, seed_mass=X26, seed_length=X27, leaf_thickness=X46, LDMC=X47, leaf_dry_mass=X55, germination_rate=X95, leaf_length=X144, leaf_width=X145, leaf_CN=X146, stem_conduit_density=X169, stem_conduit_diameter=X281, seed_number=X138, SRL=X1080)%>%
  select(-X18, -X50, -X78, -X163, -X223, -X224, -X237, -X282, -X289, -X3112, -X3113, -X3114, -X3120)

traits <- read.csv('CoRRE data\\CoRRE data\\trait data\\sCoRRE categorical trait data - traits_complete_pre spot check_03102021.csv')%>%
  full_join(contTraitsSubset) %>%
  drop_na()%>%
  filter(leaf_P_mass<20, stem_diameter<0.5, seed_mass<50, seed_number<10000, leaf_width<40, stem_conduit_density<1000, stem_conduit_diameter<200)

traitsOutliersRemoved <- traits %>%
  filter(!leaf_type %in% c("microphyll","frond")) %>%
  filter(!species_matched %in% c("Centrolepis aristata", "Centrolepis strigosa", "Acorus calamus"))

traitsScaled <- traitsOutliersRemoved %>% ## only scales continuous traits
  mutate_at(vars(ssd:SRL), scale)


# Create Gower trait disimilarity matrix
traitMatrix <- distance(traitsScaled[,15:34], method='gower') #ignoring all categorical traits

# Run PCoA
traitPCO <- pco(traitMatrix)

# Create matrix of first two axes from PCoA
PCOOutMatrix <- traitPCO$vectors[,1:2]
PCOOutMatrix_3and4 <- traitPCO$vectors[,3:4]

# Create vector fitting object (for overlaying trait vectors on top of PCoA plot) -- should probably up the permutations to 1000
trait_vf <- vf(PCOOutMatrix, traitsScaled[,15:34], nperm=100)
trait_vf_3and4 <- vf(PCOOutMatrix_3and4, traitsScaled[,15:34], nperm=100)

### Plotting
# plot(PCOOutMatrix, col=1:length(traits$species_matched),
#       pch=1:length(traits$species_matched), main="PCO", xlab="PCO 1", ylab="PCO 2")
plot(PCOOutMatrix, main="PCO", xlab="PCO 1", ylab="PCO 2", col="grey")
plot(trait_vf)

plot(PCOOutMatrix_3and4, main="PCO", xlab="PCO 3", ylab="PCO 4", col="grey")
plot(trait_vf_3and4)

# Plot with species names
PCO <- data.frame(traitPCO$vectors[,1:2], species_matched=traitsScaled$species_matched)
PCO_3and4 <- data.frame(traitPCO$vectors[,3:4], species_matched=traitsScaled$species_matched)

ggplot(data=PCO, aes(x=X1, y=X2, label=species_matched)) +
#  geom_point() +
  geom_text()
ggplot(data=PCO_3and4, aes(x=X1, y=X2, label=species_matched)) +
#  geom_point() +
  geom_text()

chart.Correlation(traitsScaled[,15:34], histogram=TRUE, pch=19)
