library(tidyverse)
library(plyr)
library(hypervolume)
library(BAT)
library(reshape2)
library(lmerTest)
library(ggpubr)

library(vegan)


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


    #minimum number of replicates
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
          subset( trt_type == "N" | trt_type == "P" | trt_type == "irr" | trt_type == "drought"| trt_type == "N*P" | trt_type == "mult_nutrient")

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
        subset( n.trt.yrs >=4)%>%
        subset(last_trt_yr == calendar_year)%>%
  subset( trt_type == "control" | trt_type == "N" | trt_type == "P" | trt_type == "irr" | 
        trt_type == "drought" | trt_type == "N*P" | trt_type == "mult_nutrient"
          )#%>%

test$trt_type <- revalue(test$trt_type, c("N*P" = "mult_nutrient"))

N <-  subset(test[test$trt_type %in% "N",], rep_num >= 5)
P <-  subset(test[test$trt_type %in% "P",], rep_num >= 5)
irr <-  subset(test[test$trt_type %in% "irr",], rep_num >= 5)
mult_nutrient <-  subset(test[test$trt_type %in% "mult_nutrient",], rep_num >= 5)
drought <-  subset(test[test$trt_type %in% "drought",], rep_num >= 4)
control <-  subset(test[test$trt_type %in% "control",], rep_num >= 4)

test <- bind_rows(N, P, irr, mult_nutrient, drought, control)
  #subset(rep_num >= 5)#%>%
#  subset(#site_code == "KNZ" |# & project_name == "change" 
    #site_code == "NWT" | 
#      site_code == "SEV" #|#&  project_name == "EDGE" & community_type == "blue_gramma"
    #site_code == "KBS" #& project_name == "T7" & community_type == "0"
    #site_code == "SEV" &  project_name == "WENNDEx"
#    )##THIS IS ONLY TO GET THE CODE TO RUN FASTER WHILE MAKING THE WORKFLOW

#test$trt_type <- revalue(test$trt_type, c("N*P" = "mult_nutrient"))

test <- test[c("site_code", "project_name", "community_type", "treatment_year", "plot_id", "genus_species", "relcov", "trt_type", "plot_mani", "treatment.x")]%>%
        unique()


plot.treatment <- test[c("site_code", "project_name", "community_type", "plot_id", "trt_type", "treatment.x")]%>%
                  unique()
plot.treatment <- tidyr::unite(plot.treatment, "rep", c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = FALSE)

df <- merge(test, traits, by.x = "genus_species", by.y = "species_matched", all.x = TRUE)

df <- unite(df, rep, c("site_code", "project_name", "community_type", "plot_id"), sep = "::", remove = FALSE)

df <- unite(df, expgroup, c("site_code", "project_name", "community_type"), sep = "::")

df$ok <- complete.cases(df[,c(#"seed_dry_mass",
  "stem_spec_density",
  #"leaf_N",
  #"leaf_P",
  "LDMC",
  "plant_height_vegetative",
  "SLA"
  #"rooting_depth"
                              )])
df <- subset(df, ok == TRUE)

#hv_split <- base::split(df[,c("seed_dry_mass", 
#                              "LDMC",
#                              "plant_height_vegetative",
#                              "rooting_depth",
#                              "relcov",
#                              "rep")], df$rep)
#hv_split <- subset(hv_split, lapply(hv_split, nrow) >1)





#hv_func <- function(x) {
#  hypervolume_gaussian(data = x[1:3], name = unique(x$rep), weight = x$relcov, 
#                       verbose = FALSE) 
#  }

#hv_volumes <- lapply(hv_split,  hv_func)



#same thing but in parallel
#library(parallel)
#numCores <- detectCores()
#cl <- makeCluster(numCores)


#clusterEvalQ(cl, {library(plyr)
#  library(hypervolume)
#  library(lmerTest)
#  library(visreg)
#  library(emmeans)
#  library(tidyverse)
#  library(ggthemes)
#  library(MuMIn)
#  library(dplyr)
#  library(BAT)
#})

#hv_volumes <- parLapply(cl,  hv_split, hv_func)


#names(hv_volumes) <- names(hv_split)                        


#hvs_joined = hypervolume_join(hv_volumes)


#make metrics?
#alpha <- sapply(hv_volumes, kernel.alpha)



#alpha <- kernel.alpha(hvs_joined)
#evenness <- kernel.evenness(hvs_joined)
#dispersion <- kernel.dispersion(hvs_joined)


#rep <- names(hv_volumes)
#hv_div <- data.frame(rep, alpha, evenness, dispersion)






##############################
####



#'test' dataframe has all the cover data but only with focal sites and treatments and such

#For each treatment at each site, pull the treatment and control data, spread, calculate distance matrix, then betadisper 

kevin <- unite(test, expgroup, c("site_code", "project_name", "community_type"), remove = FALSE, sep = "::" )

expgroup_vector <- unique(kevin$expgroup)

distances_master <- {}

for(i in 1:length(expgroup_vector)) {
  temp.df <- subset(kevin, expgroup == expgroup_vector[i])

  temp.wide <- temp.df%>%
                  pivot_wider(names_from = genus_species, values_from = relcov, values_fill = 0)
  temp.distances <- vegdist(temp.wide[10:ncol(temp.wide)], method = "bray")
  temp.mod <- betadisper(temp.distances, group = temp.wide$trt_type, type = "centroid")
  distances_temp <- data.frame(expgroup = expgroup_vector[i], trt_type = temp.wide$trt_type, treatment = temp.wide$treatment.x, plot_mani = temp.wide$plot_mani, dist = temp.mod$dist)
  distances_master <- rbind(distances_master, distances_temp )
  rm(temp.df, temp.wide, temp.distances, temp.mod, distances_temp)
}


mean.dist.df <- ddply(distances_master,.(expgroup, trt_type, treatment, plot_mani), function(x)data.frame( mean_dist = mean(x$dist)))

trt.df <- subset(mean.dist.df, plot_mani >= 1)%>%
            dplyr::rename(dist.trt = mean_dist)
con.df <- subset(mean.dist.df, plot_mani == 0)%>%
            dplyr::rename(dist.con = mean_dist)%>%
            dplyr::select(expgroup, dist.con)

lrr.df <- merge(trt.df, con.df, by = "expgroup", all.x = TRUE)%>%
    mutate(lrr = log(dist.trt/dist.con))%>%
    mutate(con_minus_trt = dist.trt/dist.con)

lrr.df.conf <- lrr.df%>%
                ddply(.(trt_type), function(x)data.frame(
                  lrr.mean = mean(x$lrr),
                  lrr.error = qt(0.975, df=length(x$trt_type)-1)*sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type)-1)
                ))

ggplot(lrr.df.conf, aes(trt_type, lrr.mean))+
  #geom_point()+
  #geom_errorbar(aes(ymin = lrr.mean-lrr.error, ymax = lrr.mean+lrr.error))+
  geom_pointrange(aes(ymin = lrr.mean-lrr.error, ymax = lrr.mean+lrr.error), size = 1.5)+
  geom_hline(yintercept = 0, size = 1)+
  xlab("")+
  ylab("Bray-Curtis LRR Distance between plots within sites")+
  theme_bw()



################################################
##########################
#########Trying that trait stuff again but with loops


expgroup_vector <- unique(df$expgroup)

tdistances_master <- {}

hv_func <- function(x) {
  hypervolume_gaussian(data = x[1:4], name = unique(x$rep), weight = x$relcov, #changing number of traits included scales time exponentially
                       chunk.size = 1000000,
                       verbose = FALSE) 
}


for(i in 1:length(expgroup_vector)) {
  temp.df <- subset(df, expgroup == expgroup_vector[i])
  temp.hv_split <- base::split(temp.df[,c(#"seed_dry_mass",
    "stem_spec_density",
    #"leaf_N",
    #"leaf_P",
    "LDMC",
    "plant_height_vegetative",
    "SLA",
    #"rooting_depth",
                                "relcov",
                                "rep")], temp.df$rep)
  temp.hv_split <- subset(temp.hv_split, lapply(temp.hv_split, nrow) >1)
  
  temp.hv_volumes <- lapply(temp.hv_split,  hv_func)
  
  temp.hvs_joined = hypervolume_join(temp.hv_volumes)
  
  temp.distance.measures <- kernel.similarity(temp.hvs_joined) #this takes a long time,
  
  
  temp.distance.centroids.dat <- data.frame(as.table(as.matrix(temp.distance.measures$Distance_centroids)))[lower.tri(as.matrix(temp.distance.measures$Distance_centroids), diag = FALSE), ]
  
  tdistances_master <- rbind(tdistances_master, temp.distance.centroids.dat )

  rm(temp.df, temp.hv_split, temp.hv_volumes)
}


tdistances_full        <-   tdistances_master%>%
  merge(plot.treatment, by.x = c("Var1"), by.y = c("rep"))%>%
  dplyr::select(Var1, Var2, Freq, trt_type, treatment.x)%>%
  dplyr::rename(c(trt_type.1 = trt_type, treatment.1 = treatment.x))%>%
  merge(plot.treatment, by.x = "Var2", by.y = "rep")%>%
    dplyr::select(Var2, Var2, Freq, trt_type.1, treatment.1, trt_type, treatment.x, site_code, project_name, community_type )%>%
  dplyr::rename(c(trt_type.2 = trt_type, treatment.2 = treatment.x))%>%
  unite( expgroup, c("site_code", "project_name", "community_type"), sep = "::", remove = TRUE)%>%
  dplyr::filter(trt_type.1 == trt_type.2)

  


tmean.dist.df <- ddply(tdistances_full,.(expgroup, trt_type.1, treatment.1), function(x)data.frame( mean_dist = mean(x$Freq)))

trt.df <- subset(tmean.dist.df, trt_type.1 != "control")%>%
  dplyr::rename(dist.trt = mean_dist)
con.df <- subset(tmean.dist.df, trt_type.1 == "control")%>%
  dplyr::rename(dist.con = mean_dist)%>%
  dplyr::select(expgroup, dist.con)

lrr.df <- merge(trt.df, con.df, by = "expgroup", all.x = TRUE)%>%
  mutate(lrr = log(dist.trt/dist.con))%>%
  mutate(con_minus_trt = dist.trt/dist.con)

#lrr.df <- read.csv("~/lrr.df_traits.csv")

lrr.df.conf <- lrr.df%>%
  ddply(.(trt_type.1), function(x)data.frame(
    lrr.mean = mean(x$lrr),
    lrr.error = qt(0.975, df=length(x$trt_type.1)-1)*sd(x$lrr, na.rm=TRUE)/sqrt(length(x$trt_type.1)-1)
  ))


ggplot(lrr.df.conf, aes(trt_type.1, lrr.mean))+
  #geom_point()+
  #geom_errorbar(aes(ymin = lrr.mean-lrr.error, ymax = lrr.mean+lrr.error))+
  geom_pointrange(aes(ymin = lrr.mean-lrr.error, ymax = lrr.mean+lrr.error), size = 1.5)+
  geom_hline(yintercept = 0, size = 1)+
  xlab("")+
  ylab("Trait space LRR Distance between plots within sites")+
  theme_bw()


lrr.df_traits <- lrr.df

write.csv(lrr.df_traits, "C:/Users/ohler/Documents/lrr.df_traits.csv")

#################################################################################
#################################################################################
##Cross site comparison

##Need some some sort of site df with cover values averaged across all the replicates

test <- unite(test, expgroup, c("site_code", "project_name", "community_type"), remove = FALSE, sep = "::" )

expgroup_vector <- unique(test$expgroup)

test_w0_master <- {}


for(i in 1:length(expgroup_vector)) {
  temp.df <- subset(test, expgroup == expgroup_vector[i])
  temp.wide <- pivot_wider(temp.df, names_from = genus_species, values_from = relcov, values_fill = 0)
  temp.test_w0 <- pivot_longer(temp.wide, cols = 10:ncol(temp.wide), names_to = "genus_species", values_to = "cover")
  
  
  test_w0_master <- rbind(test_w0_master, temp.test_w0 )
  
  rm(temp.df, temp.wide,temp.test_w0)

  
}

site.df <- test_w0_master%>%
            ddply(.(expgroup, site_code, project_name, community_type, treatment_year, trt_type, plot_mani, treatment.x, genus_species), function(x)data.frame( cover = mean(x$cover)))

#start editing here
site.traits <- merge(site.df, traits, by.x = "genus_species", by.y = "species_matched", all.x = TRUE)

site.traits <- unite(site.traits, rep, c("site_code", "project_name", "community_type", "trt_type", "treatment.x"), sep = "::", remove = FALSE)

#df <- unite(df, expgroup, c("site_code", "project_name", "community_type"), sep = "::")

site.traits$ok <- stats::complete.cases(site.traits[,c(#"seed_dry_mass",
  "stem_spec_density",
  #"leaf_N",
  #"leaf_P",
  "LDMC",
  "plant_height_vegetative",
  "SLA"
  #"rooting_depth"
                                  )])

site.traits <- subset(site.traits, ok == TRUE)



##loop to get distances

trt_type_vector <- unique(site.traits$trt_type)


tdistances_master <- {}

hv_func <- function(x) {
  hypervolume_gaussian(data = x[1:4], name = unique(x$rep), weight = x$cover, #changing number of traits included scales time exponentially
                       chunk.size = 2000,
                       verbose = FALSE) 
}


for(i in 1:length(trt_type_vector)) {
  temp.df <- subset(site.traits, trt_type == trt_type_vector[i])
  temp.hv_split <- base::split(temp.df[,c(#"seed_dry_mass",
                                          "stem_spec_density",
                                          #"leaf_N",
                                          #"leaf_P",
                                          "LDMC",
                                    "plant_height_vegetative",
                                          "SLA",
                                          #"rooting_depth",
                                          "cover",
                                          "rep")], temp.df$rep)
  temp.hv_split <- subset(temp.hv_split, lapply(temp.hv_split, nrow) >1)
  
  temp.hv_volumes <- lapply(temp.hv_split,  hv_func)
  
  temp.hvs_joined = hypervolume_join(temp.hv_volumes)
  
  temp.distance.measures <- kernel.similarity(temp.hvs_joined, num.points.max =100) #this takes a long time,
  
  
  temp.distance.centroids.dat <- data.frame(as.table(as.matrix(temp.distance.measures$Distance_centroids)))[lower.tri(as.matrix(temp.distance.measures$Distance_centroids), diag = FALSE), ]
  
  temp.distance.centroids.dat$trt_type <- trt_type_vector[i]
  
  tdistances_master <- rbind(tdistances_master, temp.distance.centroids.dat )
  
  rm(temp.df, temp.hv_split, temp.hv_volumes)
}


tdistances_full        <-   tdistances_master%>%
  separate(Var1, c("site_code", "project_name", "community_type", "trt_type", "treatment"), sep = "::", remove = FALSE)%>%
  unite(expgroup.1, c("site_code", "project_name", "community_type"), remove = TRUE, sep = "::" )%>%
  separate(Var2,c("site_code", "project_name", "community_type", "trt_type", "treatment"), sep = "::", remove = FALSE)%>%
  unite(expgroup.2, c("site_code", "project_name", "community_type"), remove = TRUE, sep = "::" )%>%
  unite(exp_pair, c("expgroup.1", "expgroup.2"))%>%
  dplyr::select(exp_pair, Freq, trt_type)


#write.csv(tdistances_full, "C:/Users/ohler/Documents/tdistances_full2.csv")
#tdistances_full <- read.csv("C:/Users/ohler/Documents/tdistances_full2.csv")

explist.mult_nutrient <- data.frame(exp_pair = unique(tdistances_full$exp_pair[tdistances_full$trt_type %in% "mult_nutrient"]))

explist.drought <- data.frame(exp_pair = unique(tdistances_full$exp_pair[tdistances_full$trt_type %in% "drought"]))

explist.P <- data.frame(exp_pair = unique(tdistances_full$exp_pair[tdistances_full$trt_type %in% "P"]))

explist.N <- data.frame(exp_pair = unique(tdistances_full$exp_pair[tdistances_full$trt_type %in% "N"]))

explist.irr <- data.frame(exp_pair = unique(tdistances_full$exp_pair[tdistances_full$trt_type %in% "irr"]))

  
finalframe.mult_nutrient <- explist.mult_nutrient%>%
                merge(tdistances_full, by = "exp_pair", all.x = TRUE )%>%
                subset(trt_type == "mult_nutrient" | trt_type == "control")

finalframe.drought <- explist.drought%>%
  merge(tdistances_full, by = "exp_pair", all.x = TRUE )%>%
  subset(trt_type == "drought" | trt_type == "control")

finalframe.P <- explist.P%>%
  merge(tdistances_full, by = "exp_pair", all.x = TRUE )%>%
  subset(trt_type == "P" | trt_type == "control")

finalframe.N <- explist.N%>%
  merge(tdistances_full, by = "exp_pair", all.x = TRUE )%>%
  subset(trt_type == "N" | trt_type == "control")

finalframe.irr <- explist.irr%>%
  merge(tdistances_full, by = "exp_pair", all.x = TRUE )%>%
  subset(trt_type == "irr" | trt_type == "control")




  
ggplot(finalframe.mult_nutrient, aes(trt_type, Freq))+
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  ylab("Trait distance between sites")+
  ggtitle("Multiple nutrients ")+
  xlab("")+
    theme_bw()

hist(subset(finalframe.mult_nutrient, trt_type == "mult_nutrient")$Freq)
  
  
ggplot(finalframe.drought, aes(trt_type, Freq))+
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  ylab("Trait distance between sites")+
  ggtitle("Drought")+
  xlab("")+
    theme_bw()

hist(subset(finalframe.drought, trt_type == "drought")$Freq)


ggplot(finalframe.P, aes(trt_type, Freq))+
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  ylab("Trait distance between sites")+
  ggtitle("Phosphorous")+
  xlab("")+
    theme_bw()

hist(subset(finalframe.P, trt_type == "P")$Freq)



ggplot(finalframe.N, aes(trt_type, Freq))+
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  ylab("Trait distance between sites")+
  ggtitle("Nitrogen")+
  xlab("")+
    theme_bw()

hist(subset(finalframe.N, trt_type == "N")$Freq)




ggplot(finalframe.irr, aes(trt_type, Freq))+
  geom_boxplot()+
  stat_compare_means(method = "t.test")+
  ylab("Trait distance between sites")+
  ggtitle("Irrigation")+
  xlab("")+
  theme_bw()

hist(subset(finalframe.irr, trt_type == "irr")$Freq)






  


