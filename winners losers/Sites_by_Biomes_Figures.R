# Goal: Plot SCoRRE sites by Biome and Experimental Type 

library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(maps)
library(remotes)
library(raster)
library(maptools)
library(cowplot)
devtools::install_github("valentinitnelav/plotbiomes")


setwd("~/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/")


location_data<-read.csv("environmental data/CoRRE_experiment_locations.csv")
experimental_data<-read.csv("community composition/CoRRE_ExperimentInfoMar2021.csv") 

experimental_data_unique<-experimental_data[!duplicated(experimental_data[c("site_code", "project_name")]),]

# Create "experimental type" columns and then go from wide to long?
View(experimental_data_unique)

experimental_data_unique$Exp_Type<-NA
experimental_data_unique$Exp_Type<-paste(experimental_data_unique$nutrients, experimental_data_unique$light, experimental_data_unique$carbon, experimental_data_unique$water, sep = "_" )

experimental_data_unique$Exp_Type_2<-ifelse(experimental_data_unique$Exp_Type == "1_0_0_0", "nutrients", 
                                     ifelse(experimental_data_unique$Exp_Type == "0_1_0_0", "light",
                                     ifelse(experimental_data_unique$Exp_Type == "0_0_1_0", "carbon",
                                     ifelse(experimental_data_unique$Exp_Type == "0_0_0_1", "water", "multiple"))))

# Need to merge lat long into dataset 

experimental_data_unique_wLatLong<-merge(experimental_data_unique, location_data, "site_code", all.y = TRUE)
experimental_data_unique_wLatLong<-experimental_data_unique_wLatLong[!is.na(experimental_data_unique_wLatLong["Latitude"]),]

# Plot studies onto Whittaker biome plot

# Get temperature and precipitation data for each study site 
r <- getData("worldclim",var="bio",res=10)
r <- r[[c(1,12)]] # annual temperature and precipitation 
names(r) <- c("Temp","Prec")

# Lat and Long for Whittacker

Lat_long_fun<- function(df, type){
  lats_ind <- df[df$Exp_Type_2 == type,]$Latitude
  longs_ind <- df[df$Exp_Type_2  == type,]$Longitude
  coords_ind <- data.frame(x=longs_ind,y=lats_ind)
  points_ind <- SpatialPoints(coords_ind, proj4string = r@crs)
  values_ind <- extract(r,points_ind)
  values_ind <- values_ind/10 # to get into mm cm and temp C 
  df_ind <- cbind.data.frame(coordinates(points_ind),values_ind)
  df_ind
}


# Get points for Whittacker plots
Lat_long_nutrients<-Lat_long_fun(experimental_data_unique_wLatLong, "nutrients")
Lat_long_carbon<-Lat_long_fun(experimental_data_unique_wLatLong, "carbon")
Lat_long_water<-Lat_long_fun(experimental_data_unique_wLatLong, "water")
Lat_long_multi<-Lat_long_fun(experimental_data_unique_wLatLong, "multiple")

# Make plot 

plot_1 <- ggplot() +
  geom_polygon(data = Whittaker_biomes,
               aes(x    = temp_c,
                   y    = precp_cm,
                   fill = biome),
               # adjust polygon borders
               colour = "gray98",
               size   = 1) +
  labs(y = "Precipitaiton (cm)", x = "Temperature (°C)")+
  theme_bw() +
  scale_fill_manual(name   = "Whittaker biomes",
                    breaks = names(Ricklefs_colors),
                    labels = names(Ricklefs_colors),
                    values = Ricklefs_colors)


Whittacker_plot_fun<-function(df1){
  plot_1 +
    geom_point(data = df1, 
               aes(x = Temp, 
                   y = Prec), 
               size   = 2.5,
               shape  = 21,
               colour = "gray90", 
               fill   = "gray45") +
    theme_bw()+ 
    theme(legend.position = "none", 
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size =14))+
    guides(shape = guide_legend(override.aes = list(size = 3)))  
}

#whit_legend<-get_legend(Whittacker_nutrients)

Whittacker_nutrients  <- Whittacker_plot_fun(Lat_long_nutrients) + ggtitle("Nutrients")
Whittacker_carbon <- Whittacker_plot_fun(Lat_long_carbon) + ggtitle("Carbon")
Whittacker_water  <- Whittacker_plot_fun(Lat_long_water) + ggtitle("Water")
Whittacker_multi <- Whittacker_plot_fun(Lat_long_multi) + ggtitle("Multiple")

whit_grid<-plot_grid(Whittacker_nutrients, Whittacker_carbon, Whittacker_water, Whittacker_multi,whit_legend, nrow=1)

pdf("~/Downloads/WhittackerPlots_SCoRRE.pdf", height = 4, width = 16)
whit_grid
dev.off()
