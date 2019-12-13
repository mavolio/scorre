library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data')

families <- read.csv('CoRRE_tax_rank.csv')

abund <- read.csv('CoRRE_relative_abundance_Nov2019.csv')%>%
  left_join(read.csv('CoRRE_TRY_species_list.csv'))%>%
  select(genus_species, species_matched, site_code, project_name, community_type, calendar_year, treatment_year, treatment, plot_id, relcov)

trait <- read.csv('TRY_trait_data_categorical.csv')%>%
  filter(CleanTraitName=='lifeform'&CleanTraitValue=='Woody')%>%
  left_join(families)%>%
  left_join(abund)

ggplot(data=trait, aes(x=relcov)) +
  geom_histogram() +
  facet_wrap(~genus, scale='free_y')
