#mixed effects models for everything?

library(tidyverse)
library(gridExtra)
library(nlme)
library(lmerTest)
library(lme4)
library(ggplot2)
library(colorspace)
library(fixest)

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
  effects_out$pval.cor <- p.adjust(p = effects_out$`Pr(>|t|)`, method = "fdr")
  effects_out
  }

#mod_function(Drought_data_exp_list$`BUX,PQ,0`, "drought")

out<-lapply(Drought_data_exp_list, mod_function, category = "drought")
out_2<-do.call(rbind, out)
out_2[,2]<-as.numeric(format(out_2[,2], scientific = FALSE))
out_2[,c(1,2)]<-round(out_2[,c(1,2)], 4)
names(out_2)[2] <- "p.val"

count.drought <- out_2 %>% mutate(mainef = ifelse(pval.cor <0.05 & estimate_type == "df[, category]", "main.sig", "non"),
                                  time = ifelse(pval.cor <0.05 & estimate_type == "treatment_year", "time.sig", "non"),
                                  main.time = ifelse(pval.cor <0.05 & estimate_type == "df[, category]:treatment_year", "int.sig","non")) %>%
  mutate (double = ifelse(mainef == "sig" & main.time == "int.sig", 1, 0)) %>%
  count(mainef, time, main.time,double, measure) %>%
  mutate(trt_type = "drought")
count.drought$prop <- 2*count.drought$n/(length(Drought_data_exp_list)*18)


# think about double counting 

###################
## nutrients #####
#################

nutrient_data<-create_exp_list(newdf2, "nuts")
nutrient_data <- nutrient_data[which(nutrient_data$site_code != "YMN"),]

nutrient_data_exp_list<-split(nutrient_data, nutrient_data$experimentid)
out<-lapply(nutrient_data_exp_list, mod_function, category = "nuts")
out_2<-do.call(rbind, out)
out_2[,2]<-as.numeric(format(out_2[,2], scientific = FALSE))
out_2[,c(1,2)]<-round(out_2[,c(1,2)], 4)
names(out_2)[2] <- "p.val"
count.nut <- out_2 %>% mutate(mainef = ifelse(pval.cor <0.05 & estimate_type == "df[, category]", "main.sig", "non"),
                              time = ifelse(pval.cor <0.05 & estimate_type == "treatment_year", "time.sig", "non"),
                              main.time = ifelse(pval.cor <0.05 & estimate_type == "df[, category]:treatment_year", "int.sig","non")) %>%
  mutate (double = ifelse(mainef == "sig" & main.time == "int.sig", 1, 0)) %>%
  count(mainef, time, main.time,double, measure) %>%
  mutate(trt_type = "nuts")
count.nut$prop <- 2*count.nut$n/(length(nutrient_data_exp_list)*18)

# Of the 18 models for "drought" three are "boudary (singular) fit" 
# check out https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#singular-models-random-effect-variances-estimated-as-zero-or-correlations-estimated-as---1 




#####################
## CO2 data ####
###################
CO2_data<-create_exp_list(newdf2, "CO2")

CO2_data_exp_list<-split(CO2_data, CO2_data$experimentid)


out<-lapply(CO2_data_exp_list, mod_function, category = "CO2")
out_2<-do.call(rbind, out)
out_2[,2]<-as.numeric(format(out_2[,2], scientific = FALSE))
out_2[,c(1,2)]<-round(out_2[,c(1,2)], 4)
names(out_2)[2] <- "p.val"
count.CO2 <- out_2 %>% mutate(mainef = ifelse(pval.cor <0.05 & estimate_type == "df[, category]", "main.sig", "non"),
                              time = ifelse(pval.cor <0.05 & estimate_type == "treatment_year", "time.sig", "non"),
                              main.time = ifelse(pval.cor <0.05 & estimate_type == "df[, category]:treatment_year", "int.sig","non")) %>%
  mutate (double = ifelse(mainef == "sig" & main.time == "int.sig", 1, 0)) %>%
  count(mainef, time, main.time,double, measure) %>%
  mutate(trt_type = "CO2")
count.CO2$prop <- 2*count.CO2$n/(length(CO2_data_exp_list)*18)





#####################
## dist data #### 
###################
dist_data<-create_exp_list(newdf2, "dist")

dist_data_exp_list<-split(dist_data, dist_data$experimentid)


out<-lapply(dist_data_exp_list, mod_function, category = "dist")
out_2<-do.call(rbind, out)
out_2[,2]<-as.numeric(format(out_2[,2], scientific = FALSE))
out_2[,c(1,2)]<-round(out_2[,c(1,2)], 4)
names(out_2)[2] <- "p.val"
count.dist <- out_2 %>% mutate(mainef = ifelse(pval.cor <0.05 & estimate_type == "df[, category]", "main.sig", "non"),
                               time = ifelse(pval.cor <0.05 & estimate_type == "treatment_year", "time.sig", "non"),
                               main.time = ifelse(pval.cor <0.05 & estimate_type == "df[, category]:treatment_year", "int.sig","non")) %>%
  mutate (double = ifelse(mainef == "sig" & main.time == "int.sig", 1, 0)) %>%
  count(mainef, time, main.time,double, measure) %>%
  mutate(trt_type = "dist")
count.dist$prop <- 2*count.dist$n/(length(dist_data_exp_list)*18)





#####################
## herb_removal data ####
###################
herb_data<-create_exp_list(newdf2, "herb_removal")

herb_data_exp_list<-split(herb_data, herb_data$experimentid)


out<-lapply(herb_data_exp_list, mod_function, category = "herb_removal")
out_2<-do.call(rbind, out)
out_2[,2]<-as.numeric(format(out_2[,2], scientific = FALSE))
out_2[,c(1,2)]<-round(out_2[,c(1,2)], 4)
names(out_2)[2] <- "p.val"
count.herb <- out_2 %>% mutate(mainef = ifelse(pval.cor <0.05 & estimate_type == "df[, category]", "main.sig", "non"),
                               time = ifelse(pval.cor <0.05 & estimate_type == "treatment_year", "time.sig", "non"),
                               main.time = ifelse(pval.cor <0.05 & estimate_type == "df[, category]:treatment_year", "int.sig","non")) %>%
  mutate (double = ifelse(mainef == "sig" & main.time == "int.sig", 1, 0)) %>%
  count(mainef, time, main.time,double, measure) %>%
  mutate(trt_type = "herb_removal")

count.herb$prop <- 2*count.herb$n/(length(herb_data_exp_list)*18)


#####################
## irg data #### 
###################
irg_data<-create_exp_list(newdf2, "irg")

irg_data_exp_list<-split(irg_data, irg_data$experimentid)


out<-lapply(irg_data_exp_list, mod_function, category = "irg")
out_2<-do.call(rbind, out)
out_2[,2]<-as.numeric(format(out_2[,2], scientific = FALSE))
out_2[,c(1,2)]<-round(out_2[,c(1,2)], 4)
names(out_2)[2] <- "p.val"
count.irg <- out_2 %>% mutate(mainef = ifelse(pval.cor <0.05 & estimate_type == "df[, category]", "main.sig", "non"),
                              time = ifelse(pval.cor <0.05 & estimate_type == "treatment_year", "time.sig", "non"),
                              main.time = ifelse(pval.cor <0.05 & estimate_type == "df[, category]:treatment_year", "int.sig","non")) %>%
  mutate (double = ifelse(mainef == "sig" & main.time == "int.sig", 1, 0)) %>%
  count(mainef, time, main.time,double, measure) %>%
  mutate(trt_type = "irg")

count.irg$prop <- 2*count.irg$n/(length(irg_data_exp_list)*18)

#####################
## temp data ####
###################
temp_data<-create_exp_list(newdf2, "temp")

temp_data_exp_list<-split(temp_data, temp_data$experimentid)


out<-lapply(temp_data_exp_list, mod_function, category = "temp")
out_2<-do.call(rbind, out)
out_2[,2]<-as.numeric(format(out_2[,2], scientific = FALSE))
out_2[,c(1,2)]<-round(out_2[,c(1,2)], 4)
names(out_2)[2] <- "p.val"
count.temp <- out_2 %>% mutate(mainef = ifelse(pval.cor <0.05 & estimate_type == "df[, category]", "main.sig", "non"),
                               time = ifelse(pval.cor <0.05 & estimate_type == "treatment_year", "time.sig", "non"),
                               main.time = ifelse(pval.cor <0.05 & estimate_type == "df[, category]:treatment_year", "int.sig","non")) %>%
  mutate (double = ifelse(mainef == "sig" & main.time == "int.sig", 1, 0)) %>%
  count(mainef, time, main.time,double, measure) %>%
  mutate(trt_type = "temp")
count.temp$prop <- 2*count.temp$n/(length(temp_data_exp_list)*18)



#####################
## multtrts data ####  
###################
multtrts_data<-create_exp_list(newdf2, "multtrts")

multtrts_data_exp_list<-split(multtrts_data, multtrts_data$experimentid)


out<-lapply(multtrts_data_exp_list, mod_function, category = "multtrts")
out_2<-do.call(rbind, out)
out_2[,2]<-as.numeric(format(out_2[,2], scientific = FALSE))
out_2[,c(1,2)]<-round(out_2[,c(1,2)], 4)
names(out_2)[2] <- "p.val"
count.multtrts <- out_2 %>% mutate(mainef = ifelse(pval.cor <0.05 & estimate_type == "df[, category]", "main.sig", "non"),
                                   time = ifelse(pval.cor <0.05 & estimate_type == "treatment_year", "time.sig", "non"),
                                   main.time = ifelse(pval.cor <0.05 & estimate_type == "df[, category]:treatment_year", "int.sig","non")) %>%
  mutate (double = ifelse(mainef == "sig" & main.time == "int.sig", 1, 0)) %>%
  count(mainef, time, main.time,double, measure) %>%
  mutate(trt_type = "multtrts")
count.multtrts$prop <- 2*count.multtrts$n/nrow(out_2) # this is weird. Not sure why i need to multiply by 2.


##Bind all the vote counting together
bigdf <- rbind(count.CO2, count.dist, count.drought, count.herb, count.irg, count.multtrts, count.nut, count.temp)
test <- bigdf %>% gather(effect, sig, mainef:main.time)
test$sig <- factor(test$sig,levels = c("main.sig", "int.sig", "time.sig", "non"))

ggplot(data = test) +
  geom_bar(aes(x = prop, y = trt_type, fill = sig), stat = "identity") +
  facet_wrap( ~ measure) +
  scale_fill_discrete_sequential(palette = "Plasma")+
  theme_classic()




##############################

mod <- lme(pd.raw ~ experimentid + nuts*treatment_year, random = ~1|experimentid/plot_id, data = nutrient_data) 
mod1 <- lme(pd.ses ~ experimentid + nuts*treatment_year, random = ~1|experimentid/plot_id, data = nutrient_data)
mod2 <- lme(mpd.raw ~ experimentid + nuts*treatment_year, random = ~1|experimentid/plot_id, data = nutrient_data)  
mod3 <- lme(mpd.ses ~ experimentid + nuts*treatment_year, random = ~1|experimentid/plot_id, data = nutrient_data)
mod4 <- lme(mntd.raw ~ experimentid + nuts*treatment_year, random = ~1|experimentid/plot_id, data = nutrient_data)  
mod5 <- lme(mntd.ses ~ experimentid + nuts*treatment_year, random = ~1|experimentid/plot_id, data = nutrient_data)  

testmod <- feols(pd.raw ~ nuts *treatment_year  | experimentid, nutrient_data) 
testmod <- feols(pd.ses ~ nuts *treatment_year  | experimentid, nutrient_data)

mod6 <- lme(pd.raw ~ experimentid + drought*treatment_year, random = ~1|experimentid/plot_id, data = Drought_data) 
mod7 <- lme(pd.ses ~ experimentid + drought*treatment_year, random = ~1|experimentid/plot_id, data = Drought_data)
mod8 <- lme(mpd.raw ~ experimentid + drought*treatment_year, random = ~1|experimentid/plot_id, data = Drought_data)  
mod9 <- lme(mpd.ses ~ experimentid + drought*treatment_year, random = ~1|experimentid/plot_id, data = Drought_data)
mod10 <- lme(mntd.raw ~ experimentid + drought*treatment_year, random = ~1|experimentid/plot_id, data = Drought_data)  
mod11 <- lme(mntd.ses ~ experimentid + drought*treatment_year, random = ~1|experimentid/plot_id, data = Drought_data)
mod12 <- lme(pd.raw ~ experimentid + CO2*treatment_year, random = ~1|experimentid/plot_id, data = CO2_data)
mod13 <- lme(pd.ses ~ experimentid + CO2*treatment_year, random = ~1|experimentid/plot_id, data = CO2_data)
mod14 <- lme(mpd.raw ~ experimentid + CO2*treatment_year, random = ~1|experimentid/plot_id, data = CO2_data)
mod15 <- lme(mpd.ses ~ experimentid + CO2*treatment_year, random = ~1|experimentid/plot_id, data = CO2_data)
mod16 <- lme(mntd.raw ~ experimentid + CO2*treatment_year, random = ~1|experimentid/plot_id, data = CO2_data)
mod17 <- lme(mntd.ses ~ experimentid + CO2*treatment_year, random = ~1|experimentid/plot_id, data = CO2_data)
View(dist_data)
mod18 <- lme(pd.raw ~ experimentid + dist*treatment_year, random = ~1|experimentid/plot_id, data = dist_data )
mod19 <- lme(pd.ses ~ experimentid + dist*treatment_year, random = ~1|experimentid/plot_id, data = dist_data)
mod20 <- lme(mpd.raw ~ experimentid + dist*treatment_year, random = ~1|experimentid/plot_id, data = dist_data)
mod21 <- lme(mpd.ses ~ experimentid + dist*treatment_year, random = ~1|experimentid/plot_id, data = dist_data)
mod22 <- lme(mntd.raw ~ experimentid + dist*treatment_year, random = ~1|experimentid/plot_id, data = dist_data)
mod23 <- lme(mntd.ses ~ experimentid + dist*treatment_year, random = ~1|experimentid/plot_id, data = dist_data)

mod24 <- lme(pd.raw ~ experimentid + irg*treatment_year, random = ~1|experimentid/plot_id, data = irg_data) 
mod25 <- lme(pd.ses ~ experimentid + irg*treatment_year, random = ~1|experimentid/plot_id, data = irg_data)
mod26 <- lme(mpd.raw ~ experimentid + irg*treatment_year, random = ~1|experimentid/plot_id, data = irg_data)  
mod27 <- lme(mpd.ses ~ experimentid + irg*treatment_year, random = ~1|experimentid/plot_id, data = irg_data)
mod28 <- lme(mntd.raw ~ experimentid + irg*treatment_year, random = ~1|experimentid/plot_id, data = irg_data)  
mod29 <- lme(mntd.ses ~ experimentid + irg*treatment_year, random = ~1|experimentid/plot_id, data = irg_data)

mod30 <- lme(pd.raw ~ experimentid + herb_removal*treatment_year, random = ~1|experimentid/plot_id, data = herb_data)
mod31 <- lme(pd.ses ~ experimentid + herb_removal*treatment_year, random = ~1|experimentid/plot_id, data = herb_data)
mod32 <- lme(mpd.raw ~ experimentid + herb_removal*treatment_year, random = ~1|experimentid/plot_id, data = herb_data)
mod33 <- lme(mpd.ses ~ experimentid + herb_removal*treatment_year, random = ~1|experimentid/plot_id, data = herb_data)
mod34 <- lme(mntd.raw ~ experimentid + herb_removal*treatment_year, random = ~1|experimentid/plot_id, data = herb_data)
mod35 <- lme(mntd.ses ~ experimentid + herb_removal*treatment_year, random = ~1|experimentid/plot_id, data = herb_data)
mod36 <- lme(pd.raw ~ experimentid + temp*treatment_year, random = ~1|experimentid/plot_id, data = temp_data)
mod37 <- lme(pd.ses ~ experimentid + temp*treatment_year, random = ~1|experimentid/plot_id, data = temp_data)
mod38 <- lme(mpd.raw ~ experimentid + temp*treatment_year, random = ~1|experimentid/plot_id, data = temp_data)
mod39 <- lme(mpd.ses ~ experimentid + temp*treatment_year, random = ~1|experimentid/plot_id, data = temp_data)
mod40 <- lme(mntd.raw ~ experimentid + temp*treatment_year, random = ~1|experimentid/plot_id, data = temp_data)
mod41 <- lme(mntd.ses ~ experimentid + temp*treatment_year, random = ~1|experimentid/plot_id, data = temp_data)
modlist <- list(mod, mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12, mod13, mod14, mod15, mod16, mod17, mod18, mod19, mod20, mod21, mod22, mod23, mod24, mod25, mod26, mod27, mod28, mod29, mod30, mod31, mod32, mod33, mod34, mod35,mod36, mod37,mod38,mod39,mod40,mod41)

coef.fun <- function(mod.no){
  x <- summary(mod.no)
  df1 <- as.data.frame(x$tTable[c(nrow(x$tTable)-2, nrow(x$tTable)-1, nrow(x$tTable)),c(1,2,5)])
  rownames(df1) <- NULL
  df1$effect <- c("main", "time", "interaction")
  df1
  }

test <- lapply(modlist,coef.fun)
out_2<-as.data.frame(do.call(rbind, test))
out_2$trt_type<- c(rep("nut",6),rep("drought",6), rep("CO2",6), rep("dist",6), rep("irg",6), rep("herb_removal",6), rep("temp",6))
out_2$pd.metric <- rep(c("pd.raw", "pd.ses", "mdp.raw", "mpd.ses", "mntd.raw", "mntd.ses"),7)




