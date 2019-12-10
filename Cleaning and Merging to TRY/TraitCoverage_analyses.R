library(tidyverse)

covdat<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/CoRRE data/trait_coverage_byplot.csv")

library(tidyverse)

trt<-read.csv("C:/Users/mavolio2/Dropbox/Sdiv_Scorre_shared/CoRRE Data/CoRRE_treatment_summary.csv")
loc<-read.csv("C:/Users/mavolio2/Dropbox/Sdiv_Scorre_shared/CoRRE Data/CoRRE_experiment_locations.csv")

alldat<-covdat%>%
  left_join(trt)%>%
  left_join(loc)

#get most coverage traits
ave_trait<-covdat%>%
  group_by(trait)%>%
  summarize(mcov=mean(coverage))%>%
  mutate(rank=rank(-mcov))%>%
  filter(rank<11)%>%
  select(trait)

#geographic range
ave_cont<-alldat%>%
  right_join(ave_trait)%>%
  group_by(Continent, project_name, community_type, trait)%>%
  summarize(mcov=mean(coverage))%>%
  ungroup%>%
  group_by(Continent, project_name, trait)%>%
  summarize(mmcov=mean(mcov))%>%
  ungroup%>%
  group_by(Continent, trait)%>%
  summarize(avecov=mean(mmcov))

## USE this df to look at subsets of interest. 

ggplot(aes(x = Continent, y = avecov), data = ave_cont) + 
  geom_bar(stat="identity") + 
  facet_wrap(~trait, ncol=5)+
  theme(axis.text.x = element_text(angle = 90))



#geographic range


