
setwd('C:\\Users\\mavolio2\\Dropbox\\CoRRE_database\\Data\\CompiledData')

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(paletteer)
library(gridExtra)

theme_set(theme_bw(12))

#mapping sites

sites<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/Species_DCiDiff_formixedmodelsMarch2024.csv") %>% 
  select(site_code, project_name, community_type, trt_type2) %>% 
  unique()

loc<-read.csv("siteLocationClimate.csv")%>%
  left_join(read.csv('SiteBiotic.csv')) %>%
  group_by(site_code) %>% 
  summarize_at(vars(c(Latitude, Longitude, MAP, MAT, rrich, anpp)), mean, na.rm=T) %>% 
  right_join(sites) %>% 
  mutate(trt_type3=factor(trt_type2, levels=c('co2', 'drought', 'irrigation', 'temp', 'n', 'p', 'multnuts', 'all mult')))

loc2<-read.csv("siteLocationClimate.csv")%>%
  left_join(read.csv('SiteBiotic.csv')) %>%
  group_by(site_code) %>% 
  summarize_at(vars(c(Latitude, Longitude, MAP, MAT, rrich, anpp)), mean, na.rm=T) 

labels<-c(
  'all mult'='Interact.',
  'co2'='CO2',
  'drought'='Drought',
  'irrigation'='Irrigation', 
  'multnuts'='Mult. Nut.',
  'n'='N',
  'p'='P',
  'temp'='Temperature')


site_byTreat<-ggplot(data=loc, aes(x=MAP, y=MAT, color=rrich, size=anpp))+
  geom_point()+
  scale_color_gradient(name="Gamma Richness", low='khaki', high='darkturquoise' )+
  scale_size_continuous(name='ANPP', breaks=c(250, 500, 750, 1000))+
  labs(x="MAP (mm)", y='MAT (\u00B0C)')+
  facet_wrap(~trt_type3, ncol=4, labeller=labeller(trt_type3=labels))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill='antiquewhite'), strip.text.x = element_text(face='bold'))
site_byTreat

# This is using code the Meghan Hayden wrote for the IDE community paper
# Load world map 
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# Map sites in each continent
map <- ggplot() + 
  geom_sf(data = world, fill = "antiquewhite") + 
  #geom_sf(data = oz_states, colour = "black", fill = NA) + 
  geom_point(data = loc, mapping = aes(x = Longitude, y = Latitude),
             pch = 21, 
             color = "black", 
             size = 2, 
             fill='black') + 
  scale_fill_paletteer_d(`"dutchmasters::milkmaid"`) +
  geom_jitter(position = "jitter") +
  theme_bw(base_size = 16) +
  labs(x = "Latitude", y = "Longitude") + 
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.25, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) + 
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "aliceblue")) + 
  coord_sf(ylim = c(-80, 80), expand = FALSE)+
  theme(legend.position = 'top')

map

figs1<-grid.arrange(map, site_byTreat)

ggsave('C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\manuscript\\Map2.jpg', figs1, units = 'in', width=7, height=8)


###mapping families
family<-read.csv('C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\CoRRE data\\trait data\\species_families_trees_2021.csv')

familydat<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/Species_DCiDiff_formixedmodelsMarch2024.csv") %>% 
  left_join(family) %>% 
  filter(family %in% c("Poaceae", "Brassicaceae", "Solanaceae", "Cyperaceae", "Polemoniaceae", "Gentianaceae", "Plantaginaceae", "Euphorbiaceae", "Amaranthaceae", "Orchidaceae", "Fabaceae", "Gentianaceae", "Orobanchaceae", "Lamiaceae")) %>% 
  select(site_code, family) %>% 
  unique() %>% 
  left_join(loc2, by='site_code')

site_byFam<-ggplot(data=familydat, aes(x=MAP, y=MAT))+
  geom_point(size=2, color='gray')+
  labs(x="MAP (mm)", y='MAT (\u00B0C)')+
  facet_wrap(~family, ncol=4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill='gray'), strip.text.x = element_text(face='bold'))
site_byFam

ggsave('C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\manuscript\\FamilyMap.jpg', site_byFam, units = 'in', width=7, height=6)


####grid of families by treatment amount
trtinfo<-read.csv('ExperimentInfo_March2024.csv') %>%
  filter(plot_mani!=0&resource_mani==1) %>% 
 select(site_code, n, p, precip, temp) %>% 
  unique()

familydatTrt<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/Species_DCiDiff_formixedmodelsMarch2024.csv") %>% 
  left_join(family) %>% 
  filter(family %in% c("Poaceae", "Brassicaceae", "Solanaceae", "Cyperaceae", "Polemoniaceae", "Gentianaceae", "Plantaginaceae", "Euphorbiaceae", "Amaranthaceae", "Orchidaceae", "Fabaceae", "Gentianaceae", "Orobanchaceae", "Lamiaceae")) %>% 
  select(site_code, family) %>% 
  unique() %>% 
  left_join(loc2, by='site_code') %>% 
  left_join(trtinfo) %>% 
  select(site_code, family, MAP, MAT, n, p, precip, temp) %>% 
  pivot_longer(n:temp, names_to = 'treatment', values_to = 'amount') %>% 
  filter(amount!=0)

site_byFamTrt<-ggplot(data=familydatTrt, aes(x=MAP, y=MAT, color=treatment))+
  geom_point(size=2)+
  geom_jitter()+
  labs(x="MAP (mm)", y='MAT (\u00B0C)')+
  facet_wrap(~family)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill='gray'), strip.text.x = element_text(face='bold'))

site_byFamTrt

Family_byTrt<-ggplot(data=familydatTrt, aes(x=MAP, y=amount, color=family))+
  geom_point(size=2)+
  #labs(x="MAP (mm)", y='MAT (\u00B0C)')+
  facet_wrap(~treatment, scales='free')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill='gray'), strip.text.x = element_text(face='bold'))
Family_byTrt

