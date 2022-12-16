################################################################################
##  sCoRRE_phylogeneticDiversity_causalModels.R: Examining differences in phylogenetic and functional diversity with causal modeling of the CoRRE database.
##
##  Author: Kimberly Komatsu
##  Date created: December 4, 2022
##  Edited 13 December 2022 (A Clark)
################################################################################

##### mixed effects model #####
#NOTE: these models do not account for biotic or abiotic env drivers at a site or for trt magnitude (but do include a random effect of site)
#library(lme4)
require(nlme)
#library(car)
#library(ggeffects)
#library(grid)
require(viridis)

allDivRR = read.csv("sCoRRE_funcphyl_RR.csv")
allDivRR$treatment_year2 = allDivRR$treatment_year^2

#options(contrasts=c('contr.sum','contr.poly')) 
data_subset = na.omit(subset(allDivRR, FDis_RR_mean<5 & trt_type2!='herb_removal' & treatment_year<50))
data_subset = data_subset[,c("FDis_RR_mean", "trt_type2", "treatment_year", "treatment_year2", "site_proj_comm")]
data_subset$trt_type2 = as.factor(data_subset$trt_type2)
data_subset$site_proj_comm = as.factor(data_subset$site_proj_comm)

FDisModel <- lme(FDis_RR_mean ~ -1 + trt_type2 + (treatment_year+treatment_year2):trt_type2,
                         data=data_subset,
                         na.action = na.omit,
                         random = ~1|site_proj_comm/trt_type2,
                         weights = varExp(form = ~fitted(.)))
summary(FDisModel)
plot(FDisModel) # weighting seems to have fixed the residual problem

# get coefficients
rand_coef = random.effects(FDisModel)
rand_names1 = row.names(rand_coef[[1]])
rand_names2 = row.names(rand_coef[[2]])

fix_coef = fixed.effects(FDisModel)
fix_names = names(fix_coef)

# make plot
category_list = unique(data_subset[,c("trt_type2", "site_proj_comm")])
trt_list = unique(data_subset[,c("trt_type2")])
category_colors = viridis(length(unique(category_list$trt_type2)))

pdf("FDis_RR_mean_temporal.pdf", width = 6, height = 4)
par(mar = c(4,4,2,2))
plot(FDis_RR_mean~jitter(treatment_year), data = data_subset,
     col = adjustcolor(category_colors[as.factor(data_subset$trt_type2)], alpha.f = 0.1),
     cex = 0.5, ylim=c(-1,1),
     ylab = "FDis RR", xlab = "Treatment Year")

# plot random effects
for(i in 1:nrow(category_list)) {
  ps = data_subset$site_proj_comm == category_list$site_proj_comm[i] &
    data_subset$trt_type2 == category_list$trt_type2[i]
  xrng = range(data_subset$treatment_year[ps], na.rm = TRUE)
  xsq = seq(xrng[1], xrng[2], length=100)
  
  # find the right model coefficients
  coef_ps1 = which(rand_names1 == as.character(category_list$site_proj_comm[i]))
  coef_ps2 = which(rand_names2 == paste(category_list$site_proj_comm[i], category_list$trt_type2[i], sep = "/"))
  rand_int = rand_coef[[1]][coef_ps1,] + rand_coef[[2]][coef_ps2,]
  
  fixef_ps = fix_coef[grep(category_list$trt_type2[i], fix_names, fixed = TRUE)]
  
  pred_ps = rand_int + fixef_ps[1] +
    fixef_ps[2]*xsq + fixef_ps[3]*xsq^2
  
  lines(xsq, pred_ps, lwd = 0.3,
        col = adjustcolor(category_colors[data_subset$trt_type2[i]], alpha.f = 0.4))
}

# plot fixed effects
for(i in 1:length(trt_list)) {
  ps = data_subset$trt_type2 == trt_list[i]
  xrng = range(data_subset$treatment_year[ps], na.rm = TRUE)
  xsq = seq(xrng[1], xrng[2], length=100)
  
  # find the right model coefficients
  fixef_ps = fix_coef[grep(trt_list[i], fix_names, fixed = TRUE)]
  
  pred_ps = fixef_ps[1] +
    fixef_ps[2]*xsq + fixef_ps[3]*xsq^2
  
  lines(xsq, pred_ps, lwd = 2,
        col = adjustcolor(category_colors[trt_list[i]], alpha.f = 0.8))
}

legend("topright", legend = as.character(trt_list), fill = category_colors, bty = "n", ncol=2)
dev.off()



Anova(FDisModel, type='III')
plot(FDisModel)
qqnorm(resid(FDisModel)) #not a great line


summary(MNTDModel <- lmer(mntd_diff_mean ~ poly(treatment_year,2)*as.factor(trt_type2) + (1 + trt_type2|site_proj_comm),
                          data=na.omit(subset(allDivRR, FDis_RR_mean<5 & trt_type2!='herb_removal' & treatment_year<50))))