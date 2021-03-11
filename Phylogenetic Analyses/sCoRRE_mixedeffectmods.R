#mixed effects models for everything?

library(tidyverse)
library(gridExtra)
library(nlme)
library(lmerTest)
library(lme4)
library(ggplot2)

###read in data

my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/"
#my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"
#my.wd <- "C:/Users/ohler/Dropbox/sDiv_sCoRRE_shared/"

#read in the data

#raw abundance data
#dat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RawAbundance_Feb2021.csv",sep=""))
#relative abundance data
#reldat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RelativeAbundance_Feb2021.csv",sep=""))
# PD data
pd <- read.csv(paste(my.wd, "paper 2_PD and FD responses/data/CoRRE_pd_metrics_nonweighted_mg.csv", sep = ""), fileEncoding = "UTF-8-BOM")

trt <- read.csv(paste(my.wd,'CoRRE data/CoRRE data/community composition/CoRRE_RawAbundanceMar2021.csv', sep = ""))%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_id, block)%>%
  unique()

trts<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfoMar2021.csv", sep=""))%>%
  select(-X)%>%
  select(site_code, project_name, community_type, treatment, trt_type, pulse, plot_mani,resource_mani)%>%
  unique()

trt_analysis<-trts%>%
  mutate(control = ifelse(trt_type == "control", 1,0))%>%
  mutate(alltrts=ifelse(trt_type %in% c("CO2","CO2*temp", "mow_clip","burn","burn*graze","disturbance","burn*mow_clip","drought","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N","mult_nutrient","N*P","P","N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp","temp*mow_clip","drought*temp","irr*temp","irr", "control"),1,0))%>%
  filter(alltrts==1)%>%
  mutate(dist=ifelse(trt_type %in% c("mow_clip","burn","burn*graze","disturbance","burn*mow_clip"), 1, 0), 
         dist_other=ifelse(trt_type %in% c("drought*mow_clip","drought*temp*mow_clip", "irr*mow_clip","irr*temp*mow_clip","N*irr*mow_clip","N*P*burn*graze","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip", "mult_nutrient*mow_clip","N*burn*mow_clip", "N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp*mow_clip"), 1, 0),
         CO2=ifelse(trt_type %in% c("CO2"), 1, 0),
         CO2_other=ifelse(trt_type %in% c("CO2*temp","drought*CO2*temp","irr*CO2","irr*CO2*temp","N*CO2*temp","N*irr*CO2","N*irr*CO2*temp","N*CO2"), 1, 0),
         drought=ifelse(trt_type %in% c("drought"), 1, 0),
         drought_other=ifelse(trt_type %in% c("drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","N*drought","drought*temp"), 1, 0),
         herb_removal=ifelse(trt_type %in% c("herb_removal"), 1, 0),
         herb_removal_other=ifelse(trt_type %in% c("herb_removal*mow_clip","irr*herb_removal","N*herb_removal","P*herb_removal","mult_nutrient*herb_removal*mow_clip","mult_nutrient*herb_removal"), 1, 0),
         irg=ifelse(trt_type %in% c("irr"), 1, 0),
         irg_other=ifelse(trt_type %in% c("irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*irr*CO2","N*irr*mow_clip","mult_nutrient*irr","N*irr*CO2*temp","N*irr","N*irr*temp","irr*temp"), 1, 0),
         temp=ifelse(trt_type %in% c("temp"), 1, 0),
         temp_other=ifelse(trt_type %in% c("control","CO2*temp","drought*CO2*temp","drought*temp*mow_clip","irr*CO2*temp","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2*temp","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","temp*mow_clip","drought*temp","irr*temp"), 1, 0),
         nuts=ifelse(trt_type %in% c("N","mult_nutrient","N*P","P"), 1, 0),
         nuts_other=ifelse(trt_type %in% c("control","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze","mult_nutrient*irr","N*irr*CO2*tempN*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip"), 1, 0),
         multtrts=ifelse(trt_type %in% c("CO2*temp", "burn*graze","burn*mow_clip","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp*mow_clip","drought*temp","irr*temp"),1,0))

dat2<-do.call(rbind, strsplit(pd$plot.id, "\\."))
dat3 <- as.data.frame(dat2)
dat3$plot.id <- pd$plot.id
dat4 <- dat3[which(dat3$V1 != dat3$V6),]
nutnet <- dat4[which(dat4$V3 == "NutNet"),]
nutnet$site_code <- paste(nutnet$V1, nutnet$V2, sep = ".")
nutnet <- nutnet[,-c(1,2)]
names(nutnet) <- c("project_name", "community_type", "treatment_year", "plot_id", "plot.id", "site_code")
dat4 <- dat4[-which(dat4$V3 == "NutNet"),]
dat4$plot_id <- paste(dat4$V5, dat4$V6, sep = ".")
dat4 <- dat4[,-c(5,6)]
names(dat4) <- c("site_code", "project_name", "community_type", "treatment_year","plot.id", "plot_id")
dat3 <- dat3[which(dat3$V1 == dat3$V6),-6]
names(dat3) <- c("site_code", "project_name", "community_type", "treatment_year", "plot_id", "plot.id")
sites <- rbind(dat3, dat4, nutnet)

newdf <- merge(pd, sites)

newdf1 <- merge(newdf, trt)
newdf2 <- merge(newdf1,trt_analysis)
newdf2$experimentid <- paste(newdf2$site_code, newdf2$project_name, newdf2$community_type, sep = ",")
newdf2$treatment_year<-as.numeric(newdf2$treatment_year)

create_exp_list<-function(df, trt_type_name){
  test <- newdf2[which(newdf2[,trt_type_name] ==1),]
  expvector <- unique(test$experimentid)
  test2 <- df%>%
    subset(trt_type == "control" & experimentid %in% expvector)
  test3 <- rbind(test, test2)
  test3
}

#####################
## Drought data ####
###################
Drought_data<-create_exp_list(newdf2, "drought")

Drought_data_exp_list<-split(Drought_data, Drought_data$experimentid)

mod_function<- function(df,category){
  mod1<-lmer(pd.raw ~ df[,category] * treatment_year + (1|plot_id), df)
  mod2<-lmer(pd.ses ~ df[,category] * treatment_year + (1|plot_id), df)
  mod3<-lmer(mpd.raw ~ df[,category] * treatment_year + (1|plot_id), df)
  mod4<-lmer(mpd.ses ~ df[,category] * treatment_year + (1|plot_id), df)
  mod5<-lmer(mntd.raw ~ df[,category] * treatment_year + (1|plot_id), df)
  mod6<-lmer(mntd.ses ~ df[,category] * treatment_year + (1|plot_id), df)
  summary_out<-summary(mod1)
  summary_out2<-summary(mod2) 
  summary_out3<-summary(mod3) 
  summary_out4<-summary(mod4) 
  summary_out5<-summary(mod5) 
  summary_out6<-summary(mod6) 
  effects_out<-as.data.frame(summary_out$coefficients[c(2,3,4),c(1,5)])
  effects_out$estimate_type<-rownames(summary_out$coefficients)[c(2,3,4)]
  effects_out$experimentid<-paste(df[1,40]) # this is sketchy
  effects_out$measure <- "pd.raw"
  effects_out2<-as.data.frame(summary_out2$coefficients[c(2,3,4),c(1,5)])
  effects_out2$estimate_type<-rownames(summary_out2$coefficients)[c(2,3,4)]
  effects_out2$experimentid<-paste(df[1,40]) # this is sketchy
  effects_out2$measure <- "pd.ses"
  effects_out3<-as.data.frame(summary_out3$coefficients[c(2,3,4),c(1,5)])
  effects_out3$estimate_type<-rownames(summary_out3$coefficients)[c(2,3,4)]
  effects_out3$experimentid<-paste(df[1,40]) # this is sketchy
  effects_out3$measure <- "mpd.raw"
  effects_out4<-as.data.frame(summary_out4$coefficients[c(2,3,4),c(1,5)])
  effects_out4$estimate_type<-rownames(summary_out4$coefficients)[c(2,3,4)]
  effects_out4$experimentid<-paste(df[1,40]) # this is sketchy
  effects_out4$measure <- "mpd.ses"
  effects_out5<-as.data.frame(summary_out5$coefficients[c(2,3,4),c(1,5)])
  effects_out5$estimate_type<-rownames(summary_out5$coefficients)[c(2,3,4)]
  effects_out5$experimentid<-paste(df[1,40]) # this is sketchy
  effects_out5$measure <- "mntd.raw"
  effects_out6<-as.data.frame(summary_out6$coefficients[c(2,3,4),c(1,5)])
  effects_out6$estimate_type<-rownames(summary_out6$coefficients)[c(2,3,4)]
  effects_out6$experimentid<-paste(df[1,40]) # this is sketchy
  effects_out6$measure <- "mntd.ses"
  effects_out <- rbind(effects_out, effects_out2, effects_out3, effects_out4, effects_out5, effects_out6)
  effects_out
  }
mod_function(Drought_data_exp_list$`BUX,PQ,0`, "drought")

out<-lapply(Drought_data_exp_list, mod_function, category = "drought")
out_2<-do.call(rbind, out)
out_2[,2]<-as.numeric(format(out_2[,2], scientific = FALSE))
out_2[,c(1,2)]<-round(out_2[,c(1,2)], 4)
names(out_2)[2] <- "p.val"
count.drought <- out_2 %>% mutate(cutoff = ifelse(p.val <0.05, "sig", "non")) %>% 
  count(estimate_type, measure, cutoff)

# think about double counting 

###################
## nutrients #####
#################

nutrient_data<-create_exp_list(newdf2, )

# Figure 

g1 <- ggplot(aes(x = treatment_year, y = mntd.raw), data = Drought_data_exp_list$`SEV,EDGE,blue_gramma`) +
  geom_point(aes (color = trt_type))+
  geom_smooth(aes (color = trt_type), method = "lm")
g2 <- ggplot(aes(x = treatment_year, y = mntd.ses), data = Drought_data_exp_list$`SEV,EDGE,blue_gramma`) +
  geom_point(aes (color = trt_type))+
  geom_smooth(aes (color = trt_type), method = "lm")
g3 <- ggplot(aes(x = treatment_year, y = mpd.raw), data = Drought_data_exp_list$`SEV,EDGE,blue_gramma`) +
  geom_point(aes (color = trt_type))+
  geom_smooth(aes (color = trt_type), method = "lm")
g4 <- ggplot(aes(x = treatment_year, y = mpd.ses), data = Drought_data_exp_list$`SEV,EDGE,blue_gramma`) +
  geom_point(aes (color = trt_type))+
  geom_smooth(aes (color = trt_type), method = "lm")
g5 <- ggplot(aes(x = treatment_year, y = pd.raw), data = Drought_data_exp_list$`SEV,EDGE,blue_gramma`) +
  geom_point(aes (color = trt_type))+
  geom_smooth(aes (color = trt_type), method = "lm")
g6 <- ggplot(aes(x = treatment_year, y = pd.ses), data = Drought_data_exp_list$`SEV,EDGE,blue_gramma`) +
  geom_point(aes (color = trt_type))+
  geom_smooth(aes (color = trt_type), method = "lm")

png(paste(my.wd, "paper 2_PD and FD responses/sample_fig.png", sep = ""), width = 6, height = 9, units = "in", res = 300)
ggpubr::ggarrange(plotlist = list(g1,g2,g3,g4,g5,g6),ncol = 2,nrow =3, legend = "top", common.legend = TRUE)
dev.off()

# Of the 18 models for "drought" three are "boudary (singular) fit" 
# check out https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#singular-models-random-effect-variances-estimated-as-zero-or-correlations-estimated-as---1 
