
#####
# Boosted Regression Trees for Phylogenetic diversity
####

#load packages:
library(ggpubr)
library(mgcv)
library(gratia)
library(gbm)
library(dismo)
library(magrittr)
library(data.table)
library(rsq)
library(tidyverse)

#set directory:
#my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/paper 2_PD and FD responses/data/"
my.wd <- "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"
#my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- 'C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\sDiv\\sDiv_sCoRRE_shared\\paper 2_PD and FD responses\\data\\' #kim's laptop

#load data:
dat<-fread(paste(my.wd, "PD_mean_effectsizes.csv", sep=""))
names(dat)

#remove rows with NAs in the response variable:
dat<-dat[complete.cases(dat[ , 6:21]),]


###N
nTrts <- dat%>%
  filter(tn==1)

mpdN <-  gbm.step(data=nTrts, gbm.x = c("MAP", "MAT", "rrich", "anpp", "n"), 
                         gbm.y = "mpd_diff_avg", family = "gaussian", #can consider other distributions
                         tree.complexity = 5, learning.rate = 0.001, #these values punish the models for over fitting
                         bag.fraction = 0.5, step.size=100, max.trees = 50000) #learning rate (0.5 is high) and step size are important
#saveRDS(mod, file = "mod.PD.rds")

#get some metrics from the model:
mpdN$gbm.call$best.trees
hist(mpdN$residuals)
contN <- mpdN$contributions
cor(nTrts$pd_diff_avg, mpdN$fitted)^2 #pseudo-R2 = 0.537
mean(mpdN$residuals * mpdN$residuals) #MSE
mpdN$cv.statistics$deviance.mean #minimum cv deviance
mpdN$cv.statistics$deviance.se #cv deviance se

nContFig <- ggplot(data=contN, aes(x=var, y=rel.inf)) + 
  geom_bar(stat='identity') +
  coord_flip() + 
  scale_x_discrete(limits = c('n', 'rrich', 'anpp', 'MAT', 'MAP')) +
  xlab('Variable') + ylab('Relative Influence')

#Partial dependence plots:
png(paste(my.wd,"N_MPD_pdp.png"),res=600,height=7,width=7,units="in"); 
gbm.plot(mpdN, smooth=T, write.title=F, y.label="Average diff MPD")
dev.off()

find.int <- gbm.interactions(mpdN) #this is the interactions
find.int$interactions
find.int$rank.list
# The returned object, here named test.int, is a list. The first 2 components summarise the results, first as a ranked list of the 5 most important pairwise interactions, and the second tabulating all pairwise interactions. The variable index numbers in $rank.list can be used for plotting.
# You can plot pairwise interactions like this:
gbm.perspec(mpdN,1,4, y.range = c(-0,20), z.range=c(0,0.6))

nDriver <- ggplot(data=nTrts, aes(x=MAP, y=mpd_diff_avg)) + geom_point() + ylab('Avg MPD Difference (trt-ctl)')


###P - can't do boosted tree because not enough data
pTrts <- dat%>%
  filter(tp==1)

summary(pMPD<-lm(mpd_diff_avg ~ p + MAP + rrich + anpp + MAT, data=pTrts))
pPartialR <- rsq.partial(pMPD)
contP <- data.frame(pPartialR$variable, pPartialR$partial.rsq)

pContFig <- ggplot(data=contP, aes(x=pPartialR.variable, y=pPartialR.partial.rsq)) + 
  geom_bar(stat='identity') +
  coord_flip() + 
  scale_x_discrete(limits = c('p', 'rrich', 'anpp', 'MAT', 'MAP')) +
  xlab('Variable') + ylab('Relative Influence')

pDriver <- ggplot(data=pTrts, aes(x=rrich, y=mpd_diff_avg)) + geom_point() + ylab('Avg MPD Difference (trt-ctl)')


# ###CO2 - can't do boosted tree because not enough data, probably shouldn't do this one at all, only 7 obs
# co2Trts <- dat%>%
#   filter(tCO2==1)
# 
# summary(co2MPD<-lm(mpd_diff_avg ~ CO2 + MAP + rrich + anpp + MAT, data=co2Trts))
# co2PartialR <- rsq.partial(co2MPD)
# contco2 <- data.frame(co2PartialR$variable, co2PartialR$partial.rsq)
# 
# co2ContFig <- ggplot(data=contco2, aes(x=co2PartialR.variable, y=co2PartialR.partial.rsq)) + 
#   geom_bar(stat='identity') +
#   coord_flip() + 
#   scale_x_discrete(limits = c('CO2', 'rrich', 'anpp', 'MAT', 'MAP')) +
#   xlab('Variable') + ylab('Relative Influence')
# 
# ggplot(data=co2Trts, aes(x=MAP, y=mpd_diff_avg)) + geom_point()


###drought - can't do boosted tree because not enough data
droTrts <- dat%>%
  filter(drought==1)

summary(droMPD<-lm(mpd_diff_avg ~ precip + MAP + rrich + anpp + MAT, data=droTrts))
droPartialR <- rsq.partial(droMPD)
contdro <- data.frame(droPartialR$variable, droPartialR$partial.rsq)

droContFig <- ggplot(data=contdro, aes(x=droPartialR.variable, y=droPartialR.partial.rsq)) + 
  geom_bar(stat='identity') +
  coord_flip() + 
  scale_x_discrete(limits = c('precip', 'rrich', 'anpp', 'MAT', 'MAP')) +
  xlab('Variable') + ylab('Relative Influence')

droDriver <- ggplot(data=droTrts, aes(x=rrich, y=mpd_diff_avg)) + geom_point() + ylab('Avg MPD Difference (trt-ctl)')



###irrigation - can't do boosted tree because not enough data
irrTrts <- dat%>%
  filter(irg==1)

summary(irrMPD<-lm(mpd_diff_avg ~ precip + MAP + rrich + anpp + MAT, data=irrTrts))
irrPartialR <- rsq.partial(irrMPD)
contirr <- data.frame(irrPartialR$variable, irrPartialR$partial.rsq)

irrContFig <- ggplot(data=contirr, aes(x=irrPartialR.variable, y=irrPartialR.partial.rsq)) + 
  geom_bar(stat='identity') +
  coord_flip() + 
  scale_x_discrete(limits = c('precip', 'rrich', 'anpp', 'MAT', 'MAP')) +
  xlab('Variable') + ylab('Relative Influence')

irrDriver <- ggplot(data=irrTrts, aes(x=anpp, y=mpd_diff_avg)) + geom_point() + ylab('Avg MPD Difference (trt-ctl)')



###disturbance - can't do boosted tree because not enough data
distTrts <- dat%>%
  filter(dist==1)

summary(distMPD<-lm(mpd_diff_avg ~ MAP + rrich + anpp + MAT, data=distTrts))
distPartialR <- rsq.partial(distMPD)
contdist <- data.frame(distPartialR$variable, distPartialR$partial.rsq)

distContFig <- ggplot(data=contdist, aes(x=distPartialR.variable, y=distPartialR.partial.rsq)) + 
  geom_bar(stat='identity') +
  coord_flip() + 
  scale_x_discrete(limits = c('rrich', 'anpp', 'MAT', 'MAP')) +
  xlab('Variable') + ylab('Relative Influence')

distDriver <- ggplot(data=distTrts, aes(x=rrich, y=mpd_diff_avg)) + geom_point() + ylab('Avg MPD Difference (trt-ctl)')



# ###herb removal - can't do boosted tree because not enough data
# herbTrts <- dat%>%
#   filter(therb_removal==1)
# 
# summary(herbMPD<-lm(mpd_diff_avg ~ MAP + rrich + anpp + MAT, data=herbTrts))
# distPartialR <- rsq.partial(distMPD)
# contherb <- data.frame(herbPartialR$variable, herbPartialR$partial.rsq)
# 
# distContFig <- ggplot(data=contdist, aes(x=distPartialR.variable, y=distPartialR.partial.rsq)) + 
#   geom_bar(stat='identity') +
#   coord_flip() + 
#   scale_x_discrete(limits = c('rrich', 'anpp', 'MAT', 'MAP')) +
#   xlab('Variable') + ylab('Relative Influence')
# 
# distDriver <- ggplot(data=distTrts, aes(x=rrich, y=mpd_diff_avg)) + geom_point() + ylab('Avg MPD Difference (trt-ctl)')



###interactions
interactionTrts <- dat%>%
  filter(multtrts==1)

mpdInteractions <-  gbm.step(data=interactionTrts, gbm.x = c("MAP", "MAT", "rrich", "anpp"), 
                  gbm.y = "mpd_diff_avg", family = "gaussian", #can consider other distributions
                  tree.complexity = 5, learning.rate = 0.001, #these values punish the models for over fitting
                  bag.fraction = 0.5, step.size=100, max.trees = 50000) #learning rate (0.5 is high) and step size are important
#saveRDS(mod, file = "mod.PD.rds")

#get some metrics from the model:
mpdInteractions$gbm.call$best.trees
hist(mpdInteractions$residuals)
mpdInteractions$contributions
contInt <- mpdInteractions$contributions
mean(mpdInteractions$residuals * mpdInteractions$residuals) #MSE
mpdInteractions$cv.statistics$deviance.mean #minimum cv deviance
mpdInteractions$cv.statistics$deviance.se #cv deviance se

interactionsContFig <- ggplot(data=contInt, aes(x=var, y=rel.inf)) + 
  geom_bar(stat='identity') +
  coord_flip() + 
  scale_x_discrete(limits = c('rrich', 'anpp', 'MAT', 'MAP')) +
  xlab('Variable') + ylab('Relative Influence')

#Partial dependence plots:
png(paste(my.wd,"interactions_MPD_pdp.png"),res=600,height=7,width=7,units="in"); 
gbm.plot(mpdInteractions, smooth=T, write.title=F, y.label="Average diff MPD")
dev.off()

find.int <- gbm.interactions(mpdInteractions) #this is the interactions
find.int$interactions
# The returned object, here named test.int, is a list. The first 2 components summarise the results, first as a ranked list of the 5 most important pairwise interactions, and the second tabulating all pairwise interactions. The variable index numbers in $rank.list can be used for plotting.
# You can plot pairwise interactions like this:
gbm.perspec(mpdN,3,5, y.range = c(-1000,20), z.range=c(0,0.6))

intDriver <- ggplot(data=interactionTrts, aes(x=anpp, y=mpd_diff_avg)) + geom_point() + ylab('Avg MPD Difference (trt-ctl)')

ggarrange(nContFig, pContFig, droContFig, irrContFig, distContFig, interactionsContFig,
          ncol = 6, nrow = 1)
#export at 2000x400

ggarrange(nDriver, pDriver, droDriver, irrDriver, distDriver, intDriver,
          ncol = 6, nrow = 1)
#export at 2000x400







# #####
# #run BRTs:
# 
# #PD:
# set.seed(123) #for reproducibility
# 
# mod <-  gbm.step(data=dat, gbm.x = c("MAP", "MAT", "rrich", "anpp", "n", "p", "k", "CO2", "precip", "temp", "herb_removal"), 
#                  gbm.y = "pd_diff_avg", family = "gaussian", #can consider other distributions
#                  tree.complexity = 5, learning.rate = 0.001, #these values punish the models for over fitting
#                  bag.fraction = 0.5, step.size=100, max.trees = 20000) #learning rate (0.5 is high) and step size are important
# #saveRDS(mod, file = "mod.PD.rds")
# 
# #get some metrics from the model:
# mod$gbm.call$best.trees
# hist(mod$residuals)
# mod$contributions
# cor(dat$pd_diff_avg, mod$fitted)^2 #pseudo-R2
# mean(mod$residuals * mod$residuals) #MSE
# mod$cv.statistics$deviance.mean #minimum cv deviance
# mod$cv.statistics$deviance.se #cv deviance se
# 
# #Partial dependence plots:
# png("PD_pdp.png",res=600,height=7,width=7,units="in"); 
# gbm.plot(mod, smooth=T, write.title=F, y.label="Average diff PD")
# dev.off()
# 
# 
# 
# 
# #adding treatment in the model to look at relative importance
# dat$treatment <- as.factor(dat$treatment)
# mod2 <-  gbm.step(data=dat, gbm.x = c("treatment", "MAP", "MAT", "rrich", "anpp", "n", "p", "k", "CO2", "precip", "temp", "herb_removal"), 
#                  gbm.y = "pd_diff_avg", family = "gaussian",
#                  tree.complexity = 5, learning.rate = 0.0001, 
#                  bag.fraction = 0.5, step.size=100, max.trees = 20000)
# 
# mod2$contributions
# cor(dat$pd_diff_avg, mod2$fitted)^2 #pseudo-R2
# 
# ##################### from Elith paper: Interrogate and plot the interactions
# # This code assesses the extent to which pairwise interactions exist in the data.
# 
# find.int <- gbm.interactions(mod2) #this is the interactions
# 
# 
# # The returned object, here named test.int, is a list. The first 2 components summarise the results, first as a ranked list of the 5 most important pairwise interactions, and the second tabulating all pairwise interactions. The variable index numbers in $rank.list can be used for plotting.
# # You can plot pairwise interactions like this:
# 
# gbm.perspec(mod,2,1, y.range = c(15,20), z.range=c(0,0.6))
# 
# #only treatment & MAP
# 
# mod3 <-  gbm.step(data=dat, gbm.x = c("treatment", "MAP"), 
#                   gbm.y = "pd_diff_avg", family = "gaussian",
#                   tree.complexity = 5, learning.rate = 0.0001, 
#                   bag.fraction = 0.5, step.size=100, max.trees = 20000)
# 
# mod3$contributions
# 
# cor(dat$pd_diff_avg, mod3$fitted)^2 #pseudo-R2
# 
# find.int <- gbm.interactions(mod3)
# 
# 
# 
# 
# 
# #MPD:
# set.seed(123) #for reproducibility
# mod <-  gbm.step(data=dat, gbm.x = c("MAP", "MAT", "rrich", "anpp", "n", "p", "k", "CO2", "precip", "temp", "herb_removal"), 
#                   gbm.y = "mpd_diff_avg", family = "gaussian",
#                   tree.complexity = 5, learning.rate = 0.001, 
#                   bag.fraction = 0.5, step.size=100, max.trees = 20000)
# #saveRDS(mod, file = "mod.MPD.rds")
# 
# #get some metrics from the model:
# mod$gbm.call$best.trees
# hist(mod$residuals)
# mod$contributions
# cor(dat$pd_diff_avg, mod$fitted)^2 #pseudo-R2
# mean(mod$residuals * mod$residuals) #MSE
# mod$cv.statistics$deviance.mean #minimum cv deviance
# mod$cv.statistics$deviance.se #cv deviance se
# 
# #Partial dependence plots:
# png("MPD_pdp.png",res=600,height=7,width=7,units="in"); 
# gbm.plot(mod, smooth=T, write.title=F, y.label="Average diff MPD")
# dev.off()
# 
# 
# 
# #MNTD:
# set.seed(123) #for reproducibility
# mod <-  gbm.step(data=dat, gbm.x = c("MAP", "MAT", "rrich", "anpp", "n", "p", "k", "CO2", "precip", "temp", "herb_removal"), 
#                  gbm.y = "mntd_diff_avg", family = "gaussian",
#                  tree.complexity = 5, learning.rate = 0.001, 
#                  bag.fraction = 0.5, step.size=100, max.trees = 20000)
# #saveRDS(mod, file = "mod.MNTD.rds")
# 
# #get some metrics from the model:
# mod$gbm.call$best.trees
# hist(mod$residuals)
# mod$contributions
# cor(dat$pd_diff_avg, mod$fitted)^2 #pseudo-R2
# mean(mod$residuals * mod$residuals) #MSE
# mod$cv.statistics$deviance.mean #minimum cv deviance
# mod$cv.statistics$deviance.se #cv deviance se
# 
# #Partial dependence plots:
# png("MNTD_pdp.png",res=600,height=7,width=7,units="in"); 
# gbm.plot(mod, smooth=T, write.title=F, y.label="Average diff MNTD")
# dev.off()
# 
# 
# #plot(mpd_diff_avg ~ treatment, data=dat[dat$treatment == "N",])
# #abline(h=0)
# 
# #clean-up:
# rm(list = ls())
