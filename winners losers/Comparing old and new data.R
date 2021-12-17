setwd("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\data\\")

old<-read.table("Species_DCiDiff_newtrts_filtered_Jul2021.csv", sep="")%>%
  select(species_matched, trt_type2, ave_diff)
new<-read.csv("Species_DCiDiff_Dec2021.csv")%>%
  select(species_matched, trt_type2, ave_diff)
new$species_matched<-gsub(" ", "_", new$species_matched)


library(tidyverse)

allmult_old<-old%>%
  filter(trt_type2=="all mult") %>% 
  mutate(old=1)

allmult_new<-new %>% 
  filter(trt_type2=="all mult") %>% 
  mutate(new=1)%>%
  mutate(newdiff=ave_diff) %>% 
  select(-ave_diff)

diff<-allmult_old%>%
  anti_join(allmult_new)

diff<-allmult_old%>%
  full_join(allmult_new)

plot(diff$ave_diff, diff$newdiff)
cor.test(diff$ave_diff, diff$newdiff)
