library(tidyverse)
library(plyr)
library(hypervolume)
library(BAT)
library(reshape2)
library(lmerTest)
library(ggpubr)




#Read in data
traits <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/trait data/Final Cleaned Traits/Continuous_Traits/Backtrans_GapFilled_sCorre.csv")

cover <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Dec2021.csv")

experimentinfo <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfo_Dec2021.csv")


#create average trait values per species
traits <- ddply(traits,.(species_matched),
      function(x)data.frame(
        seed_dry_mass = mean(x$seed_dry_mass),
        stem_spec_density = mean(x$stem_spec_density),
        leaf_N = mean(x$leaf_N),
        leaf_P = mean(x$leaf_P),
        LDMC = mean(x$LDMC),
        leaf_C = mean(x$leaf_C),
        leaf_dry_mass = mean(x$leaf_dry_mass),
        plant_height_vegetative = mean(x$plant_height_vegetative),
        leaf_C.N = mean(x$leaf_C.N),
        SLA = mean(x$SLA),
        water_content = mean(x$water_content),
        rooting_depth = mean(x$rooting_depth),
        SRL = mean(x$SRL),
        seed_number = mean(x$seed_number)
      ))



#standardize the scale of all the traits
cols <- c( "seed_dry_mass", 
            "stem_spec_density",
           "leaf_N",
           "leaf_P",
           "LDMC",
           "leaf_C",
           "leaf_dry_mass",
           "plant_height_vegetative",
           "leaf_C.N",
           "SLA",
           "water_content",
           "rooting_depth",
           "SRL",
           "seed_number")

traits[cols] <- scale(traits[cols])

#make the species column merge-able
traits$species_matched <- tolower(traits$species_matched)




#Reduce cover data to focal data using a series of merges


    #minimum number of replicates (I guess I still haven't figured this part out yet)
repnum <- cover%>%
  dplyr::select(site_code, project_name, community_type, treatment, plot_id)%>%
  unique()%>%
  #dplyr::mutate(present = 1)%>%
  dplyr::group_by(site_code, project_name, community_type, treatment)%>%
  dplyr::summarise(rep_num = length(plot_id))%>%
  dplyr::ungroup()


  #only certain manipulations
foctrt <- experimentinfo[c("site_code", "project_name", "community_type", "trt_type")]%>%
          unique()%>%
          subset( trt_type == "N" | trt_type == "P" | trt_type == "irr" | trt_type == "drought")

foctrt$treats_wanted <- foctrt$trt_type
foctrt <- foctrt[c("site_code", "project_name", "community_type", "treats_wanted")]



  #Last year of experiment
lastyear <- ddply(experimentinfo, .(site_code, project_name, community_type),
      function(x)data.frame(
        last_trt_yr = max(x$calendar_year)
      ))


  #minimum treatment length
nyear <- experimentinfo[c("site_code", "project_name", "community_type", "treatment_year")] %>%
        unique()%>%
        subset(treatment_year != 0)%>%
        ddply(.(site_code, project_name, community_type),
              function(x)data.frame(
                n.trt.yrs = max(x$treatment_year)
              ))



#Merge all the datasets above to create columns so subset by
crest <- cover %>%
  merge(nyear,by = c("site_code", "project_name", "community_type"))%>%
  merge(lastyear, by = c("site_code", "project_name", "community_type"))%>%
  merge( experimentinfo, by = c("site_code", "project_name", "community_type", "treatment", "calendar_year", "treatment_year"), all.x = TRUE)%>%
  merge(foctrt, by = c("site_code", "project_name", "community_type"), all.x = TRUE)%>%
  merge(repnum, by = c("site_code", "project_name", "community_type"), all.x = TRUE)


#subset by criteria
test <- crest %>%
  subset(treats_wanted != "NA")%>%
        subset( n.trt.yrs >=5)%>%
        subset(last_trt_yr == calendar_year)%>%
  subset( trt_type == "control" | trt_type == "N" | trt_type == "P" | trt_type == "irr" | 
            trt_type == "drought")%>%
  subset(rep_num >= 5)%>%
  subset(site_code == "KNZ" & project_name == "change" 
    #site_code == "NWT" | 
      #site_code == "SEV" &  project_name == "EDGE" & community_type == "blue_gramma"
    #site_code == "KBS" & project_name == "T7" & community_type == "0"
    )##THIS IS ONLY TO GET THE CODE TO RUN FASTER WHILE MAKING THE WORKFLOW

test <- test[c("site_code", "project_name", "community_type", "treatment_year", "plot_id", "genus_species", "relcov", "trt_type")]%>%
        unique()


plot.treatment <- test[c("site_code", "project_name", "community_type", "plot_id", "trt_type")]%>%
                  unique()
plot.treatment <- tidyr::unite(plot.treatment, "rep", c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = FALSE)

df <- merge(test, traits, by.x = "genus_species", by.y = "species_matched", all.x = TRUE)

df <- unite(df, rep, c("site_code", "project_name", "community_type", "plot_id"), sep = "::")

df$ok <- complete.cases(df[,c("seed_dry_mass", 
                              #"LDMC",
                              "plant_height_vegetative",
                              "rooting_depth"
                              )])
df <- subset(df, ok == TRUE)

hv_split <- base::split(df[,c("seed_dry_mass", 
                              #"LDMC",
                              "plant_height_vegetative",
                              "rooting_depth",
                              "relcov",
                              "rep")], df$rep)
hv_split <- subset(hv_split, lapply(hv_split, nrow) >1)





hv_func <- function(x) {
  hypervolume_gaussian(data = x[1:3], name = unique(x$rep), weight = x$relcov, 
                       verbose = FALSE) 
  }

#hv_volumes <- lapply(hv_split,  hv_func)



#same thing but in parallel
library(parallel)
numCores <- detectCores()
cl <- makeCluster(numCores)


clusterEvalQ(cl, {library(plyr)
  library(hypervolume)
  library(lmerTest)
  library(visreg)
  library(emmeans)
  library(tidyverse)
  library(ggthemes)
  library(MuMIn)
  library(dplyr)
  library(BAT)
})

hv_volumes <- parLapply(cl,  hv_split, hv_func)


#names(hv_volumes) <- names(hv_split)                        


hvs_joined = hypervolume_join(hv_volumes)


#make metrics?
#alpha <- sapply(hv_volumes, kernel.alpha)



#alpha <- kernel.alpha(hvs_joined)
#evenness <- kernel.evenness(hvs_joined)
#dispersion <- kernel.dispersion(hvs_joined)


#rep <- names(hv_volumes)
#hv_div <- data.frame(rep, alpha, evenness, dispersion)

#hv_div <- separate(hv_div, rep, c("site_code", "project_name", "community_type", "plot_id"), sep = "-")

#hv_div <- unite(hv_div, expgroup, c("site_code", "project_name", "community_type"), sep = "-", remove = FALSE)

#hv_div <- merge(hv_div, plot.treatment, by = c("site_code", "project_name", "community_type", "plot_id"))


#ggplot(hv_div, aes(x=trt_type, y = alpha))+
#  facet_wrap(~expgroup)+
#  geom_boxplot()+
#  ylim(0,3)+
#  theme_bw()


#dist_func <- function(x) {
#  hypervolume_distance(hv1 = x, hv2 = (x+1)) 
#}

#lapply(hv_volumes, dist_func)


distance.measures <- kernel.similarity(hvs_joined) #this takes a long time, but once it's done I hope to be able to subset dataframe by sites and treatments which is the last step for comparison

#distance.centroids <- distance$Distance_centroids


distance.centroids.dat <- data.frame(as.table(as.matrix(distance.measures$Distance_centroids)))[lower.tri(as.matrix(distance.measures$Distance_centroids), diag = FALSE), ]
  
distance.centroids.dat$rep1 <- distance.centroids.dat$Var1
distance.centroids.dat$rep2 <- distance.centroids.dat$Var2
distance.centroids.dat$value <- distance.centroids.dat$Freq
#  as.matrix(distance.measures$Distance_centroids)%>%
#  reshape2::melt( varnames = c("rep1", "rep2"))%>%
#    unique()


scooby <- distance.centroids.dat%>%
          merge(plot.treatment, by.x = "rep1", by.y = "rep")

dooby <- scooby[c("rep1", "rep2", "value", "trt_type", "site_code", "project_name", "community_type")]%>%
        unite( expgroup1, c("site_code", "project_name", "community_type"), sep = "::")

dooby$trt_rep1 <- dooby$trt_type

doo <- dooby%>%
      merge( plot.treatment, by.x = "rep2", by.y = "rep")%>%
      unite( expgroup2, c("site_code", "project_name", "community_type"), sep = "::")

doo$trt_rep2 <- doo$trt_type.y

doo <- doo[c("rep1", "rep2", "value", "trt_rep1", "trt_rep2", "expgroup1", "expgroup2")]


last <- doo%>%
          unite( trt_comp, c("trt_rep1", "trt_rep2"), sep = ".", remove = FALSE)%>%
          subset( value != 0)%>%
          subset( trt_comp == "control.control" | trt_comp == "drought.drought" | trt_comp == "N.N" | trt_comp == "P.P" | trt_comp == "irr.irr" )%>%
          mutate( match = ifelse(expgroup1 == expgroup2 , TRUE, FALSE))

withinsite <- subset(last, match == TRUE)

ggplot(withinsite, aes(trt_comp, value))+
  facet_wrap(~expgroup1)+
  geom_boxplot()+
  ylab("Distance between centroids")+
  stat_compare_means(method = "t.test")+
  theme_bw()





##Figure out how to do inter-site comparisons. All the distances are already calculated. Have to use site pairs as a random effect?


last    <- last%>%
        unite(exppair, c("expgroup1", "expgroup2" ), remove = FALSE, sep = "||")


mod <- lmer(value~trt_comp + (1|exppair), data = last)
summary(mod)

library(visreg)
visreg(mod, ylab = "Distance between centroids")







##############################
####
library(vegan)



#'test' dataframe has all the cover data but only with focal sites and treatments and such

#For each treatment at each site, pull the treatment and control data, spread, calculate distance matrix, then betadisper 


wide <- test%>%
    pivot_wider(names_from = genus_species, values_from = relcov, values_fill = 0)

distances <- dist(wide[7:ncol(wide)])
mod <- betadisper(distances, group = wide$trt_type, type = "centroid")
permutest(mod, permutations = 99)
anova(mod)
#plot(mod)
boxplot(mod)
