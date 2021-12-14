d = read.csv("~/Dropbox/SharedFolders/sDiv_sCoRRE_shared/Trait Data/TRY Data/TRY Continuous data/TRY_trait_data_continuous_long_Nov2021.csv")
tmp = d[grep("Larrea", d$species_matched),]
hist(tmp[tmp$CleanTraitName=="leaf_N",]$StdValue)
hist(d[d$CleanTraitName=="leaf_N",]$StdValue)
abline(v = mean(tmp[tmp$CleanTraitName=="leaf_N",]$StdValue))


tmp_N = with(d[d$CleanTraitName=="leaf_N",], tapply(StdValue, species_matched, mean))
tmp_SLA = with(d[d$CleanTraitName=="SLA",], tapply(StdValue, species_matched, mean))
tmp_RD = with(d[d$CleanTraitName=="rooting_depth",], tapply(StdValue, species_matched, mean))


plot_data = data.frame(species = sort(unique(c(names(tmp_N), names(tmp_SLA), names(tmp_RD)))), N = NA, SLA = NA, RD = NA)
plot_data$N = tmp_N[match(plot_data$species, names(tmp_N))]
plot_data$SLA = tmp_SLA[match(plot_data$species, names(tmp_SLA))]
plot_data$RD = tmp_RD[match(plot_data$species, names(tmp_RD))]
plot_data = plot_data[!is.na(rowSums(plot_data[,-1])),]

require(rgl)
plot3d(plot_data$RD, plot_data$N, plot_data$SLA, type = "h", zlim = c(0, 50))





d2 = read.csv("~/Dropbox/SharedFolders/sDiv_sCoRRE_shared/CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Dec2021.csv")
head(d2)

species_drought = sort(unique(d2[d2$treatment=="drought",]$genus_species))
species_control = sort(unique(d2[d2$treatment=="control",]$genus_species))
species_control = species_control[!species_control%in%species_drought]



plot_data$indrought = match(tolower(plot_data$species), tolower(species_drought))
plot_data$incontrol = match(tolower(plot_data$species), tolower(species_control))


psd = !is.na(plot_data$indrought)
plot3d(plot_data$RD[psd], plot_data$N[psd], plot_data$SLA[psd], type = "h", zlim = c(0, 50))

psc = !is.na(plot_data$incontrol)
plot3d(plot_data$RD[psc], plot_data$N[psc], plot_data$SLA[psc], type = "h", zlim = c(0, 50))


cor(plot_data[psd,2:4])
cor(plot_data[psc,2:4])




require(vegan)
# standardise to zero mean and unit sd before analysis
D1 = apply(plot_data[c(which(psc), which(psd)),2:4],2,function(x) (x-mean(x))/sd(x))

iris.pca1 <- rda(D1)

## pca1
iris.pca1$CA$v # the PCA vectors - these are the formula for relationships among traits

# overall explanatory power of each axis
iris.pca1$CA$eig/sum(iris.pca1$CA$eig) # explanatory power of axes based on eigenvalues

get_pca_exp = function(pca, data, dimensions, observations = 1:nrow(data)) {
  # r2-like index for PCA
  cor(c(pca$CA$u[observations,dimensions,drop=FALSE] %*% t(pca$CA$v[,dimensions,drop=FALSE])), c(data[observations,]))^2
}

# run across all data in D1
get_pca_exp(iris.pca1, D1, 1) # dimension 1
get_pca_exp(iris.pca1, D1, 2) # dimension 2
get_pca_exp(iris.pca1, D1, 3) # dimension 3
# matches explanatory power from eigenvalues

# run on controls
get_pca_exp(iris.pca1, D1, 1, observations = 1:60) # dimension 1
get_pca_exp(iris.pca1, D1, 2, observations = 1:60) # dimension 2
get_pca_exp(iris.pca1, D1, 3, observations = 1:60) # dimension 3



# run on drought
get_pca_exp(iris.pca1, D1, 1, observations = 61:nrow(D1)) # dimension 1
get_pca_exp(iris.pca1, D1, 2, observations = 61:nrow(D1)) # dimension 2
get_pca_exp(iris.pca1, D1, 3, observations = 61:nrow(D1)) # dimension 3






D2 = apply(plot_data[c(which(psc)),2:4],2,function(x) (x-mean(x))/sd(x))
D3 = apply(plot_data[c(which(psd)),2:4],2,function(x) (x-mean(x))/sd(x))
iris.pca2 <- rda(D2)
iris.pca3 <- rda(D3)

## pca2
iris.pca2$CA$v # the PCA vectors - these are the formula for relationships among traits
iris.pca3$CA$v # the PCA vectors - these are the formula for relationships among traits


# overall explanatory power of each axis
iris.pca2$CA$eig/sum(iris.pca2$CA$eig) # explanatory power of axes based on eigenvalues
iris.pca3$CA$eig/sum(iris.pca3$CA$eig) # explanatory power of axes based on eigenvalues

require(MatrixCorrelation)
cor(iris.pca2$CA$v, iris.pca3$CA$v)




