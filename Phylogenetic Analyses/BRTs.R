
#####
# Boosted Regression Trees for Phylogenetic diversity
####

#load packages:
library(mgcv)
library(gratia)
library(gbm)
library(dismo)
library(magrittr)
library(data.table)

#set directory:
#my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/"
my.wd <- "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"
#my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"

#load data:
dat<-fread(paste(my.wd, "PD_mean_effectsizes.csv", sep=""))
names(dat)

#remove rows with NAs in the response variable:
dat<-dat[complete.cases(dat[ , 19:21]),]



#####
#run BRTs:

#PD:
set.seed(123) #for reproducibility
mod <-  gbm.step(data=dat, gbm.x = c("MAP", "MAT", "rrich", "anpp", "n", "p", "k", "CO2", "precip", "temp", "herb_removal"), 
                     gbm.y = "pd_diff_avg", family = "gaussian",
                     tree.complexity = 5, learning.rate = 0.001, 
                     bag.fraction = 0.5, step.size=100, max.trees = 20000)
#saveRDS(mod, file = "mod.PD.rds")

#get some metrics from the model:
mod$gbm.call$best.trees
hist(mod$residuals)
mod$contributions
cor(dat$pd_diff_avg, mod$fitted)^2 #pseudo-R2
mean(mod$residuals * mod$residuals) #MSE
mod$cv.statistics$deviance.mean #minimum cv deviance
mod$cv.statistics$deviance.se #cv deviance se

#Partial dependence plots:
png("PD_pdp.png",res=600,height=7,width=7,units="in"); 
gbm.plot(mod, smooth=T, write.title=F, y.label="Average diff PD")
dev.off()



#MPD:
set.seed(123) #for reproducibility
mod <-  gbm.step(data=dat, gbm.x = c("MAP", "MAT", "rrich", "anpp", "n", "p", "k", "CO2", "precip", "temp", "herb_removal"), 
                  gbm.y = "mpd_diff_avg", family = "gaussian",
                  tree.complexity = 5, learning.rate = 0.001, 
                  bag.fraction = 0.5, step.size=100, max.trees = 20000)
#saveRDS(mod, file = "mod.MPD.rds")

#get some metrics from the model:
mod$gbm.call$best.trees
hist(mod$residuals)
mod$contributions
cor(dat$pd_diff_avg, mod$fitted)^2 #pseudo-R2
mean(mod$residuals * mod$residuals) #MSE
mod$cv.statistics$deviance.mean #minimum cv deviance
mod$cv.statistics$deviance.se #cv deviance se

#Partial dependence plots:
png("MPD_pdp.png",res=600,height=7,width=7,units="in"); 
gbm.plot(mod, smooth=T, write.title=F, y.label="Average diff MPD")
dev.off()



#MNTD:
set.seed(123) #for reproducibility
mod <-  gbm.step(data=dat, gbm.x = c("MAP", "MAT", "rrich", "anpp", "n", "p", "k", "CO2", "precip", "temp", "herb_removal"), 
                 gbm.y = "mntd_diff_avg", family = "gaussian",
                 tree.complexity = 5, learning.rate = 0.001, 
                 bag.fraction = 0.5, step.size=100, max.trees = 20000)
#saveRDS(mod, file = "mod.MNTD.rds")

#get some metrics from the model:
mod$gbm.call$best.trees
hist(mod$residuals)
mod$contributions
cor(dat$pd_diff_avg, mod$fitted)^2 #pseudo-R2
mean(mod$residuals * mod$residuals) #MSE
mod$cv.statistics$deviance.mean #minimum cv deviance
mod$cv.statistics$deviance.se #cv deviance se

#Partial dependence plots:
png("MNTD_pdp.png",res=600,height=7,width=7,units="in"); 
gbm.plot(mod, smooth=T, write.title=F, y.label="Average diff MNTD")
dev.off()

#clean-up:
rm(list = ls())
