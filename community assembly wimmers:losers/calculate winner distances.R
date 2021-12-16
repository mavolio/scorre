
#load packages
library(gawdis)

#set directory
my.wd <- "~/Dropbox/sCoRRE/sDiv_sCoRRE_shared"
my.wd <- "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/" #Padu's wd


#load data:
corre<-read.table(paste(my.wd, "CoRRE_RelativeCover_Dec2021.csv", sep=""), header=T, sep=",", fill = TRUE)
trait<-read.table(paste(my.wd, "Backtrans_GapFilled_sCorre.csv", sep=""), header=T, sep=",", fill = TRUE)

#get mean trait values:
trait <- trait %>%  group_by(species_matched) %>% summarise_at(vars("seed_dry_mass","stem_spec_density", "leaf_N",                 
                                                                  "leaf_P", "LDMC", "leaf_C","leaf_dry_mass",          
                                                                   "plant_height_vegetative", "leaf_C.N", "SLA","water_content",          
                                                                   "rooting_depth", "SRL"), median, na.rm=T) %>% as.data.frame() 
rownames(trait)<-trait$species_matched
trait$species_matched<-NULL

#calculate functional dissimilarities between species using GAWDIS:
#I haven't log-transformed traits. Pending!!
#trait.dis<-gawdis(trait, w.type="optimized") #This took forever and I had to stop it. Can remove "w.type="optiimized" but a warning pops up)
#trait.dis<-gawdis(trait) #with "analytical" method
trait.dis<-cluster::daisy(trait, metric="gower") #or with gower distance

#convert distance object into phylotree using the Ward's algorithm (https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust):
ftree<-as.phylo(hclust(trait.dis, method="average"))

#save output:
write.tree(ftree, paste(my.wd, "ftree.scorre.gowdis.log.upgma.tre", sep=""))





# only experimennts that are at least 5 years

##############################################
# identify winners and losers per experiment #
##############################################
# DCi though time for control and treatment and then slopes in Adam's framework + extinct and immigration species
# important!remove too rare species; maybe identify via CSi values, do histograms to identify threshlds but it could be around 0.1


# Emily developed code to generate DCi through time: “emily ambient v2.R”







##############################################
# calculate mdns, dnns, wmdns
##############################################

# first step in one experiment
ab_data_nutnet <- ab_data[(ab_data$project_name %in% c("NutNet") & ab_data$treatment %in% c("Control", "NPK")),]
nutnet_sp_list <- unique(ab_data_nutnet[c("genus_species")])
nutnet_traits <- trait_data_cont[trait_data_cont$species_matched %in% c(nutnet_sp_list$genus_species),]    

nutnet_traits_ss <- nutnet_traits[c("species_matched", "LDMC", "leaf_N", "SLA", "seed_dry_mass", "plant_height_generative",
                                    "plant_height_vegetative", "root_density", "rooting_depth")]

# take average trait values 
df_traits_avg <- function(df_traits_ss){
  aggregate(list(LDMC = df_traits_ss$LDMC, 
                 leaf_N = df_traits_ss$leaf_N, 
                 SLA = df_traits_ss$SLA, 
                 seed_dry_mass = df_traits_ss$seed_dry_mass, 
                 plant_height_generative = df_traits_ss$plant_height_generative, 
                 plant_height_vegetative = df_traits_ss$plant_height_vegetative,
                 root_density = df_traits_ss$root_density, 
                 rooting_depth = df_traits_ss$rooting_depth), 
            by = list(species_matched = df_traits_ss$species_matched), FUN = mean, na.rm = TRUE)
}

nutnet_traits_avg <- df_traits_avg(nutnet_traits_ss)

# merge with abundance values 
nutnet_comm_traits <- merge(ab_data_nutnet, nutnet_traits_avg, by.x = "genus_species", by.y = "species_matched", all.x = TRUE)













#