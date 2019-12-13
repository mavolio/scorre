### Assess community weighted traits across MAT gradient
###
### Author: Wilcox (kevin.wilcox@uwyo.edu)
### Created Dec 2019

### Set up workspace
# Kevin's laptop
setwd("C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\CoRRE data\\")

library(tidyverse)
library(ggthemes)

### Read in merged TRY and CoRRE data
source("C:\\Users\\wilco\\Dropbox\\Cross_workstation_workspace\\Git projects\\scorre\\Cleaning and Merging to TRY\\merging corre abundance and try trait data.R")
rm(species_name_key, trait_df)

### Calculate relative abundances of only species we have traits for in each plot
unique(abundance_traits_df$CleanTraitName)
unique(abundance_traits_df$type)

## Leaf area
total_abundance_per_plot <- raw_abundance_df %>%
  group_by(site_code, project_name, community_type, calendar_year, block, treatment, plot_id) %>%
  summarize(total_cover = sum(abundance))

leaf_area_abundance_sums <- abundance_traits_df %>%
  filter(CleanTraitName == "leaf_area") %>%
  filter(type=="identified species") %>%
  group_by(site_code, project_name, community_type, calendar_year, block, treatment, plot_id) %>%
  summarize(sum_cover_with_leaf_area = sum(abundance)) %>%
  ungroup() %>%
  full_join(total_abundance_per_plot, by=c("site_code", "project_name", "community_type", "calendar_year", "block","treatment", "plot_id")) %>%
  mutate(rel_cover_with_leaf_area = sum_cover_with_leaf_area/total_cover)

hist(leaf_area_abundance_sums$rel_cover_with_leaf_area)

problem_sets <- leaf_area_abundance_sums %>%
  filter(rel_cover_with_leaf_area > 1.75)

one_problem <- abundance_traits_df %>%
  filter(CleanTraitName == "leaf_area") %>%
  filter(type=="identified species") %>%
  filter(calendar_year==2001 & block==2 & treatment == "1_0_CO" & project_name=="Exp1" & site_code=="ASGA")

another_problem <- abundance_traits_df %>%
  filter(CleanTraitName == "leaf_area") %>%
  filter(type=="identified species") %>%
  filter(calendar_year==2013 & treatment == 9 & project_name=="e001" & site_code=="CDR" & community_type=="A" & plot_id==46)

unique(another_problem$species_matched)
# test <- abundance_traits_df %>%
#   filter(CleanTraitName == "leaf_area") %>%
#   filter(type=="identified species") %>%
#   filter(site_code == "CHY" & project_name == "EDGE" & calendar_year == 2016 & plot_id == 17)






