library(tidyverse)
library(plyr)
library(hypervolume)
library(BAT)




traits <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/trait data/Final Cleaned Traits/Continuous_Traits/Backtrans_GapFilled_sCorre.csv")

cover <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Dec2021.csv")

experimentinfo <- read.csv("C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfo_Dec2021.csv")



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

traits$species_matched <- tolower(traits$species_matched)




#Reduce cover data to focal data using a series of right merges

    #Last year of experiment
    #only certain manipulations
    #minimum treatment length
    #minimum number of replicates


foctrt <- experimentinfo[c("site_code", "project_name", "community_type", "trt_type")]%>%
          unique()%>%
          subset( trt_type == "N" | trt_type == "P" | trt_type == "irr" | trt_type == "drought")

foctrt$treats_wanted <- foctrt$trt_type
foctrt <- foctrt[c("site_code", "project_name", "community_type", "treats_wanted")]


lastyear <- ddply(experimentinfo, .(site_code, project_name, community_type),
      function(x)data.frame(
        last_trt_yr = max(x$calendar_year)
      ))


#temp <- experimentinfo[c("site_code", "project_name", "community_type", "treatment_year")]
#unique(temp)

nyear <- experimentinfo[c("site_code", "project_name", "community_type", "treatment_year")] %>%
        unique()%>%
        subset(treatment_year != 0)%>%
        ddply(.(site_code, project_name, community_type),
              function(x)data.frame(
                n.trt.yrs = length(x$treatment_year)
              ))



temp <- experimentinfo[c("site_code", "project_name", "community_type", "trt_type")] %>%
        unique()

#temp <- cover[c("site_code", "project_name", "community_type", "plot_id")] %>%
#        merge(experimentinfo, by = )
#  unique()%>%
 

test <- cover %>%
  #merge( temp, by = c("site_code", "project_name", "community_type"), all.y = TRUE)%>%
  merge(nyear,by = c("site_code", "project_name", "community_type"))%>%
  merge(lastyear, by = c("site_code", "project_name", "community_type"))%>%
  merge( experimentinfo, by = c("site_code", "project_name", "community_type", "treatment", "calendar_year", "treatment_year"), all.x = TRUE)%>%
  merge(foctrt, by = c("site_code", "project_name", "community_type"), all.x = TRUE)%>%
  subset(treats_wanted != "NA")
#%>%
#  subset(trt_type == "N" | trt_type == "P" | trt_type == "irr" | trt_type == "drought")






test <- test %>%
        subset( n.trt.yrs >=10)%>%
        subset(last_trt_yr == calendar_year)%>%
  subset( trt_type == "control" | trt_type == "N" | trt_type == "P" | trt_type == "irr" | trt_type == "drought")

test <- test[c("site_code", "project_name", "community_type", "treatment_year", "plot_id", "genus_species", "relcov", "trt_type")]%>%
        unique()





df <- merge(test, traits, by.x = "genus_species", by.y = "species_matched", all.x = TRUE)

df <- unite(df, rep, c("site_code", "project_name", "community_type", "plot_id"), sep = "-")

df$ok <- complete.cases(df[,c("seed_dry_mass", 
                             "stem_spec_density",
                             "leaf_N",
                             "leaf_P",
                             "LDMC")])
df <- subset(df, ok == TRUE)

hv_split <- base::split(df[,c("seed_dry_mass", 
                              "stem_spec_density",
                              "leaf_N",
                              "leaf_P",
                              "LDMC",
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



alpha <- kernel.alpha(hvs_joined)
evenness <- kernel.evenness(hvs_joined)
dispersion <- kernel.dispersion(hvs_joined)


rep <- names(hv_volumes)
hv_div <- data.frame(rep, alpha)

hv_div <- separate(hv_div, rep, c("site_code", "project_name", "community_type", "plot_id"), sep = "-")







