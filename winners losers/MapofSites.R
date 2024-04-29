
setwd('C:\\Users\\mavolio2\\Dropbox\\CoRRE_database\\Data\\CompiledData')

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(paletteer)
library(gridExtra)

theme_set(theme_bw(12))

sites<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/Species_DCiDiff_formixedmodelsNov22.csv") %>% 
  select(site_code, project_name, community_type, trt_type2) %>% 
  unique()

loc<-read.csv("siteLocationClimate.csv")%>%
  left_join(read.csv('SiteBiotic.csv')) %>%
  group_by(site_code) %>% 
  summarize_at(vars(c(Latitude, Longitude, MAP, MAT, rrich, anpp)), mean, na.rm=T) %>% 
  right_join(sites) %>% 
  mutate(trt_type3=factor(trt_type2, levels=c('co2', 'drought', 'irrigation', 'temp', 'n', 'p', 'multnuts', 'all mult')))

labels<-c(
  'all mult'='Interact.',
  'co2'='CO2',
  'drought'='Drt',
  'irrigation'='Irg.', 
  'multnuts'='Mult. Nut.',
  'n'='N',
  'p'='P',
  'temp'='Temp.')


all<-ggplot(data=loc, aes(x=MAP, y=MAT, color=rrich, size=anpp))+
  geom_point()+
  scale_color_gradient(name="Gamma Richness", low='khaki', high='darkturquoise' )+
  scale_size_continuous(name='ANPP', breaks=c(250, 500, 750, 1000))+
  labs(x="MAP (mm)", y='MAT (\u00B0C)')+
  facet_wrap(~trt_type3, ncol=4, labeller=labeller(trt_type3=labels))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill='antiquewhite'), strip.text.x = element_text(face='bold'))
all


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

figs1<-grid.arrange(map, all)

ggsave('C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\WinnersLosers paper\\manuscript\\Map.jpg', figs1, units = 'in', width=7, height=8)
