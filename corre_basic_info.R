library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\CoRRE data')

correProjectInfo <- read.csv('CoRRE_project_summary.csv')
correTrtInfo <- read.csv('CoRRE_treatment_summary.csv')
correRelAbund <- read.csv('CoRRE_relative_abundance_Nov2019.csv')

#number of total plots, species, observations (plot*years)
all <- correRelAbund%>%
  left_join(correTrtInfo)%>%
  left_join(correProjectInfo)
plots <- all%>%
  select(site_code, project_name, community_type, plot_id)%>%
  unique()
years <- all%>%
  select(site_code, project_name, community_type, plot_id, calendar_year)%>%
  unique()


#histogram of experiment lengths
ggplot(data=correProjectInfo, aes(x=experiment_length)) +
  geom_histogram() + xlab('Experiment Length (yrs)') + ylab('Count')
summary(correProjectInfo$experiment_length)

#list of treatment types and number of experiments for which they are present
trtTypes <- correTrtInfo%>%
  select(site_code, project_name, community_type, treatment, n, p, k, CO2, precip, temp, mow_clip, burn, herb_removal, management, other_trt, plant_mani)%>%
  unique()
