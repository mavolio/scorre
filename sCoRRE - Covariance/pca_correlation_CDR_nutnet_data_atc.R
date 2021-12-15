## get data and run PCA
require(vegan)

# clean up data 
ab_data <- read.csv("~/Dropbox/SharedFolders/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_RawAbundance_Dec2021.csv")
trait_data_cont <- read.csv("~/Dropbox/SharedFolders/sDiv_sCoRRE_shared/Trait Data/TRY Data/TRY Continuous data/TRY_trait_data_continuous_Nov2021.csv")
trait_data_cont$species_matched <- tolower(trait_data_cont$species_matched)

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

# pull the last year of the experiment from each site 
nutnet_splits <- split(nutnet_comm_traits, nutnet_comm_traits$site_code)
nutnet_comm_traits_lastyr <-  lapply(nutnet_splits, function (x) {
  x = x; x$yr_group = ifelse(x$treatment_year == max(x$treatment_year),"last", "other");
  y = x[x$yr_group =="last",]; return (y)})

# get unique species by traits at *site/treatment* level 
sps_function <- function(df){
  df_comm_traits_lastyr_sp <- unique(df[c("site_code", "genus_species", "treatment","LDMC", "leaf_N", "SLA",
                                          "plant_height_generative", "plant_height_vegetative", "root_density", "rooting_depth")])
  return(df_comm_traits_lastyr_sp)
}

site_sps_out_nutnet <- lapply(nutnet_comm_traits_lastyr, sps_function)
site_sps_out_nutnet_all <- do.call(rbind, site_sps_out_nutnet)

# work with CDR ?

CDR_nutnet <- site_sps_out_nutnet$CDR

# standardise to zero mean and unit sd before analysis
D1 = apply(site_sps_out_nutnet_all[,4:10],2,function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))
D1 = apply(D1, 2, function(x) {x[!is.finite(x)] = 0; x})
treatment = site_sps_out_nutnet_all$treatment
sites = site_sps_out_nutnet_all$site_code


# drop rows without any values 
iris.pca1 <- rda(D1)

## pca1
iris.pca1$CA$v # the PCA vectors - these are the formula for relationships among traits

# overall explanatory power of each axis
iris.pca1$CA$eig/sum(iris.pca1$CA$eig) # explanatory power of axes based on eigenvalues
cumsum(iris.pca1$CA$eig/sum(iris.pca1$CA$eig))

# predict SLA traits from PC1
pca_axis_to_use = 1
sla_from_pc1 = iris.pca1$CA$v[3,pca_axis_to_use]*iris.pca1$CA$u[,pca_axis_to_use]

plot(sla_from_pc1, D1[,"SLA"])
summary(lm(D1[,"SLA"]~sla_from_pc1))


plot(sla_from_pc1[treatment=="Control"], D1[treatment=="Control","SLA"])
points(sla_from_pc1[treatment=="NPK"], D1[treatment=="NPK","SLA"], col = 2)

summary(lm(D1[treatment=="Control","SLA"]~sla_from_pc1[treatment=="Control"]))
summary(lm(D1[treatment=="NPK","SLA"]~sla_from_pc1[treatment=="NPK"]))


# predict SLA traits from multiple PCA axis
pca_axis_to_use = c(1:7)
sla_from_pc14 = iris.pca1$CA$u[,pca_axis_to_use]%*%t(iris.pca1$CA$v[3,pca_axis_to_use,drop=FALSE])

plot(sla_from_pc14[treatment=="Control"], D1[treatment=="Control","SLA"])
summary(lm(D1[treatment=="Control","SLA"]~sla_from_pc14[treatment=="Control"]))
plot(sla_from_pc14[treatment=="NPK"], D1[treatment=="NPK","SLA"])
summary(lm(D1[treatment=="NPK","SLA"]~sla_from_pc14[treatment=="NPK"]))


plot(sla_from_pc14[sites=="Bt"], D1[sites=="Bt","SLA"])
summary(lm(D1[sites=="Bt","SLA"]~sla_from_pc14[sites=="Bt"]))
plot(sla_from_pc14[sites=="CDR"], D1[sites=="CDR","SLA"])
summary(lm(D1[sites=="CDR","SLA"]~sla_from_pc14[sites=="CDR"]))



# now work on two traits at once from one pca axis
pca_axis_to_use = 1
traits_to_use = c(1,3)
ldmc_from_pc1 = iris.pca1$CA$v[traits_to_use,pca_axis_to_use][1]*iris.pca1$CA$u[,pca_axis_to_use]
sla_from_pc1 = iris.pca1$CA$v[traits_to_use,pca_axis_to_use][2]*iris.pca1$CA$u[,pca_axis_to_use]

plot(ldmc_from_pc1, sla_from_pc1)
(lm(ldmc_from_pc1~sla_from_pc1))
iris.pca1$CA$v[traits_to_use,pca_axis_to_use][1]/iris.pca1$CA$v[traits_to_use,pca_axis_to_use][2]



# predict SLA traits from multiple PCA axis
pca_axis_to_use = c(1:6)
traits_to_use = c(1,3)
ldmc_sla_from_pc123 = iris.pca1$CA$u[,pca_axis_to_use]%*%t(iris.pca1$CA$v[traits_to_use,pca_axis_to_use])

plot(ldmc_sla_from_pc123[,1], D1[,traits_to_use[1]])
plot(ldmc_sla_from_pc123[,2], D1[,traits_to_use[2]])

plot(ldmc_sla_from_pc123[,1], ldmc_sla_from_pc123[,2])
summary(lm(ldmc_sla_from_pc123[,1]~ldmc_sla_from_pc123[,2]))

plot(D1[,traits_to_use[1]], D1[,traits_to_use[2]])



plot(ldmc_sla_from_pc123[sites=="Bt",1], ldmc_sla_from_pc123[sites=="Bt",2])
summary(lm(ldmc_sla_from_pc123[sites=="Bt",2]~ldmc_sla_from_pc123[sites=="Bt",1]))
plot(ldmc_sla_from_pc123[sites=="CDR",1], ldmc_sla_from_pc123[sites=="CDR",2])
summary(lm(ldmc_sla_from_pc123[sites=="CDR",2]~ldmc_sla_from_pc123[sites=="CDR",1]))


plot(ldmc_sla_from_pc123[treatment=="Control",1], ldmc_sla_from_pc123[treatment=="Control",2], xlab="LDMC", ylab="SLA")
abline(lm(ldmc_sla_from_pc123[treatment=="Control",2]~ldmc_sla_from_pc123[treatment=="Control",1]))
points(ldmc_sla_from_pc123[treatment=="NPK",1], ldmc_sla_from_pc123[treatment=="NPK",2],col=2)
abline(lm(ldmc_sla_from_pc123[treatment=="NPK",2]~ldmc_sla_from_pc123[treatment=="NPK",1]),col=2)



get_pca_exp = function(pca, data, dimensions, observations = 1:nrow(data)) {
  # r2-like index for PCA
  cor(c(pca$CA$u[observations,dimensions,drop=FALSE] %*% t(pca$CA$v[,dimensions,drop=FALSE])), c(data[observations,]))^2
}

# run across all data in D1
get_pca_exp(iris.pca1, D1, 1) # dimension 1
get_pca_exp(iris.pca1, D1, 2) # dimension 2
get_pca_exp(iris.pca1, D1, 3) # dimension 3
get_pca_exp(iris.pca1, D1, 4) # dimension 4
# matches explanatory power from eigenvalues

# run on subsets of data (just first 10 observations)
get_pca_exp(iris.pca1, D1, 1, observations = 1:10) # dimension 1
get_pca_exp(iris.pca1, D1, 2, observations = 1:10) # dimension 2
get_pca_exp(iris.pca1, D1, 3, observations = 1:10) # dimension 3
get_pca_exp(iris.pca1, D1, 4, observations = 1:10) # dimension 4


## calculate correspondence between pca1 and pca2
require(MatrixCorrelation)
PSI(iris.pca1$CA$v, iris.pca2$CA$v, center = TRUE) # very high correlation
