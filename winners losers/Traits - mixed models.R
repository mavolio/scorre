######
#####
##### code to study how species traits affect response to GCDs
##### code by M.Avolio with help from Adam Clark and Tamara Munkemuller
#### created summer 2021, updated with new trait dataset June 20, 2022
######

library(tidyverse)
library(lme4)
library(emmeans)
#library(relaimpo)

theme_set(theme_bw(12))

### Trait data

# #barGraphStats(data=, variable="", byFactorNames=c(""))
# barGraphStats <- function(data, variable, byFactorNames) {
#   count <- length(byFactorNames)
#   N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
#   names(N)[1:count] <- byFactorNames
#   names(N) <- sub("^x$", "N", names(N))
#   mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
#   names(mean)[1:count] <- byFactorNames
#   names(mean) <- sub("^x$", "mean", names(mean))
#   sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
#   names(sd)[1:count] <- byFactorNames
#   names(sd) <- sub("^x$", "sd", names(sd))
#   preSummaryStats <- merge(N, mean, by=byFactorNames)
#   finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
#   finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
#   return(finalSummaryStats)
# }  

#read in data
contTraits1 <- read.csv('C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\CoRRE data\\trait data\\Final TRY Traits\\Imputed Continuous_Traits\\data to play with\\imputed_continuous_20220620.csv')

contTraits<-contTraits1%>%
  select(-X.1, -X, -family, -genus, -observation)%>%
  select(species_matched, LDMC, SLA, plant_height_vegetative, rooting_depth, seed_dry_mass) %>% 
  group_by(species_matched)%>%
  summarise_all(funs(mean))%>%
  ungroup()  %>% 
  filter(seed_dry_mass<30, plant_height_vegetative<10, rooting_depth<3, SLA<75)

# # ##code to drop outliers
# contLong<-contTraits %>%
#   pivot_longer(LDMC:seed_dry_mass, names_to = "trait", values_to = "value")
# ggplot(data=contLong, aes(x=value))+
#   geom_histogram()+
#   facet_wrap(~trait, scales = "free")

#why do I need to scale? Kim, Kevin and I discussed this on 8/10/22 and decided I did not have to do this.
# traitsScaled <- contTraits%>%
#   mutate_at(vars(LDMC:seed_dry_mass), scale)



catTraits <- read.csv('C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\CoRRE data\\trait data\\Final TRY Traits\\sCoRRE categorical trait data_final_20211209.csv') %>% 
  select(species_matched, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation) %>% 
  filter(growth_form!="moss",species_matched!="", growth_form!="lycophyte") %>% 
  mutate(mycorrhizal=ifelse(mycorrhizal_type=="none", 'no', ifelse(mycorrhizal_type=="facultative_AM"|mycorrhizal_type=="facultative_AM_EcM", "facultative", ifelse(mycorrhizal_type=="uncertain", "unk", "yes")))) %>% 
  select(-mycorrhizal_type) %>% 
  mutate(photo_path=ifelse(photosynthetic_pathway=="possible C4"|photosynthetic_pathway=="possible C4/CAM", "C4", ifelse(photosynthetic_pathway=="possible CAM", "CAM",photosynthetic_pathway))) %>% 
  select(-photosynthetic_pathway)

pairs(contTraits[,2:6])


# Read in dci diff

# dcidiff<-read.csv("C:/Users/megha/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/Species_DCiDiff_newtrts.csv")

dcidiff_models<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/Species_DCiDiff_formixedmodelsNov22.csv")

test<-dcidiff_models %>% 
   filter(trt_type2=="all mult") %>% 
   select(species_matched) %>% 
  unique()

alldat_cont<-dcidiff_models%>%
  right_join(contTraits)%>%
  gather(LDMC:seed_dry_mass, key="trait", value="value")%>%
  na.omit()


#Making a graph of all the models I am goign to run now, just to see if there are patterns.
ggplot(data=alldat_cont, aes(x=value, y=diff))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(trt_type2~trait, scales="free")
###
#making mixed models - run through all traits.

##-1 says don't give me an overall intercept. just changing the output. I want an intercept per level of each trt_type2
##fit a seperate slope for the value term for each of trt_type levels
## by using trt|species - come up with a value for each species, but do that for each trt type category
#(trt_type|species matches) = get a separate estimate for every species in each trt type
#(1|species matches) = each species can have a different ave regardless of trt type
#fixef(m1) #should give fixed effects of model

##SLA
mSLA<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="SLA"))
summary(mSLA)

plot.sla<-as.data.frame(summary(mSLA)$coefficients)

plot.sla<-plot.sla%>%
  mutate(fixedef=row.names(plot.sla))

toplot.SLA<-plot.sla%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="SLA")

#not using SRL b/c data is no good on TRY.
# ###SRL
# 
# mSRL<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="SRL"))
# summary(mSRL)
# 
# plot.srl<-as.data.frame(summary(mSRL)$coefficients)%>%
#   mutate(fixedef=row.names(plot.srl))
# 
# toplot.SRL<-plot.srl%>%
#   separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
#   filter(!is.na(interaction))%>%
#   separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
#   select(-interaction, -drop)%>%
#   mutate(trait="SRL")

###seed mass
#this model failed to converge what does that mean?
mSM<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="seed_dry_mass"))
summary(mSM)

plot.sm<-as.data.frame(summary(mSM)$coefficients)

plot.sm<-plot.sm %>% mutate(fixedef=row.names(plot.sm))

toplot.SM<-plot.sm%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="Seed Mass")


# ###seed number - not doing this the data was weird.
# 
# mSN<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="seed_number"))
# summary(mSN)
# 
# plot.sn<-as.data.frame(summary(mSN)$coefficients)%>%
#   mutate(fixedef=row.names(plot.sn))
# 
# toplot.SN<-plot.sn%>%
#   separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
#   filter(!is.na(interaction))%>%
#   separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
#   select(-interaction, -drop)%>%
#   mutate(trait="Seed Number")

###rooting depth

mRD<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="rooting_depth"))
summary(mRD)

plot.rd<-as.data.frame(summary(mRD)$coefficients)

plot.rd<-plot.rd%>%
  mutate(fixedef=row.names(plot.rd))

toplot.RD<-plot.rd%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="Rooting Depth")

###LDMC

mLDMC<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="LDMC"))
summary(mLDMC)

plot.ldmc<-as.data.frame(summary(mLDMC)$coefficients)

plot.ldmc<-plot.ldmc%>%
  mutate(fixedef=row.names(plot.ldmc))

toplot.LDMC<-plot.ldmc%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="LDMC")

###plot height vegetative

mhght<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="plant_height_vegetative"))
summary(mhght)

plot.hght<-as.data.frame(summary(mhght)$coefficients)

plot.hght<-plot.hght%>%
  mutate(fixedef=row.names(plot.hght))

toplot.hgt<-plot.hght%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="PlantHeight")


toplot<-toplot.SLA%>%
  bind_rows(toplot.SM, toplot.RD, toplot.LDMC, toplot.hgt)

colnames(toplot)[2] <- "SE"

toplot$trt_type2=factor(toplot$trt_type, levels=c("co2", "drought", "irrigation", "temp", "n", "p", "multnuts", "all mult"))

#gave up here trying to nicely label figures
trt.labels=c("CO2", "Drt", "Irg", "Temp", "N", "P", "Nutrients", "Interacting")

ggplot(data=toplot, aes(y=Estimate, x=1))+
  geom_point()+
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width=0.05)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank())+
  xlab("")+
  scale_x_continuous(limits=c(0, 2))+
  geom_hline(yintercept=0, linetype="dashed")+
  facet_grid(trt_type2~trait, scales="free", labeller = labeller(trt_type2=trt.labels))


####Categorical data

table(catTraits$photo_path)#drop possible hybird, parasitic
table(catTraits$growth_form)#drop cactus
table(catTraits$clonal)#drop uncertain
table(catTraits$n_fixation)
table(catTraits$mycorrhizal)#drop unk
table(catTraits$lifespan)#drop uncertain

alldat_cat<-dcidiff_models%>%
  right_join(catTraits)%>%
  gather(growth_form:photo_path, key="trait", value="value")%>%
  na.omit()

#photo path
mpp<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="photo_path"&value!="parasitic"&value!="possible hybrid"))
summary(mpp)

plot.mpp<-as.data.frame(emmeans(mpp, ~ value*trt_type2)) %>% 
  mutate(trait="photo_path")

#lifespan
ml<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="lifespan"&value!='uncertain'))
summary(ml)

plot.ml<-as.data.frame(emmeans(ml, ~ value*trt_type2)) %>% 
  mutate(trait="lifespan")

#growth form
mgf<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="growth_form"&value!="cactus"))
summary(mgf)

plot.mgf<-as.data.frame(emmeans(mgf, ~ value*trt_type2)) %>% 
  mutate(trait="growth_form")

###clonaltiy
mc<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="clonal"&value!='uncertain'))

plot.mc<-as.data.frame(emmeans(mc, ~ value*trt_type2))%>%
  mutate(trait="Clonal")

##nfixer
mnf<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="n_fixation"))

plot.nf<-as.data.frame(emmeans(mnf, ~ value*trt_type2))%>%
  mutate(trait="Nfix")


##mycorhizal
mmyc<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="mycorrhizal"&value!="unk"))

plot.myc<-as.data.frame(emmeans(mmyc, ~ value*trt_type2))%>%
  mutate(trait="Myc")


toplot.cat<-plot.mpp%>%
  bind_rows(plot.ml, plot.mgf, plot.mc, plot.nf, plot.myc)

toplotesacat<-toplot.cat

##for now need to subset out each cat trait one by one
ggplot(data=subset(toplotesacat, trait=="lifespan"), aes(y=emmean, x=1))+
  geom_point()+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.05)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank())+
  xlab("")+
  scale_x_continuous(limits=c(0, 2))+
  geom_hline(yintercept=0, linetype="dashed")+
  facet_grid(trt_type2~value, scales="free")

ggplot(data=ESA, aes(y=Estimate, x=1))+
  geom_point()+
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width=0.05)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank())+
  xlab("")+
  scale_x_continuous(limits=c(0, 2))+
  geom_hline(yintercept=0, linetype="dashed")+
  facet_grid(Trt~trait, scales="free")

####making figures of it all together
toplot.cat2<-toplot.cat %>% 
  filter(trait!="growth_form") %>% 
    mutate(value2=paste(trait, value, sep="_" )) %>% 
  select(-trait) %>% 
  rename(trait=value2, trt_type=trt_type2, Estimate=emmean) %>% 
  filter(trait!="growth_form"&trait!="Nfix_no"&trait!='Myc_yes'&trait!='lifespan_biennial'&trait!='photo_path_CAM') %>% 
  select(-asymp.LCL, -asymp.UCL, -df, -value)

toplot.cont<-toplot %>% 
  select(-"t value")
         
alltraitmodels<-toplot.cont %>% 
  bind_rows(toplot.cat2) %>% 
  mutate(Trait_name=ifelse(trait=='Clonal_no', "Non-Clonal", ifelse(trait=='Clonal_yes', 'Clonal', ifelse(trait=='lifespan_annual', 'Annual', ifelse(trait=='lifespan_perennial', 'Perennial', ifelse(trait=='Myc_facultative', 'Mycorrhizal', ifelse(trait=='Myc_no', 'Non-Mycorrhizal', ifelse(trait=='Nfix_yes', 'Nitrogen Fixer', ifelse(trait=='photo_path_C3', "C3", ifelse(trait=='photo_path_C4', 'C4', ifelse(trait=='PlantHeight', 'Plant Height', trait))))))))))) %>% 
  mutate(Traits2=factor(Trait_name,levels=c("Annual", 'Perennial', 'Clonal', 'Non-Clonal', 'C3', 'C4', 'Mycorrhizal', 'Non-Mycorrhizal', 'Nitrogen Fixer', 'LDMC', 'Plant Height', 'Rooting Depth', 'Seed Mass', 'SLA')))


theme_set(theme_bw(16))
ggplot(data=subset(alltraitmodels), aes(y=Estimate, x=1))+
  geom_point()+
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width=0.05)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank())+
  xlab("")+
  scale_x_continuous(limits=c(0, 2))+
  geom_hline(yintercept=0, linetype="dashed")+
  facet_grid(trt_type~Traits2, scales="free")


###trying to do a multiple regression

combined<-dcidiff_models %>% 
  right_join(contTraits) %>% 
  right_join(catTraits) %>% 
  drop_na() %>% 
  filter(photo_path!='parasitic'&photo_path!='possible hybrid') %>% 
  filter(mycorrhizal!='unk') %>% 
  filter(clonal!='uncertain') %>% 
  filter(lifespan!='uncertain') %>% 
  mutate(photo_path2=ifelse(photo_path=="CAM"|photo_path=='C4', 'CAM/C4', photo_path)) %>% 
  mutate(lifespan2=ifelse(lifespan=='annual'|lifespan=='biennial', 'ann/bien', lifespan)) %>% 
  mutate(mycorrhizal=as.factor(mycorrhizal),
         lifespan2=as.factor(lifespan2),
         clonal=as.factor(clonal),
         n_fixation=as.factor(n_fixation),
         photo_path2=as.factor(photo_path2))

sp<-dcidiff_models %>% 
  select(species_matched) %>% 
  unique()

str(combined) #check categorial variables are factors

##testing a multiple regression###
##just putting in each model what came out significant in previous analyses

##co2 model
mCO2<-lm(diff~SLA+plant_height_vegetative+rooting_depth+seed_dry_mass+photo_path2+n_fixation, data=subset(combined, trt_type2=="co2"))

summary(mCO2)
calc.relimp(mCO2)

##N model
mN<-lm(diff~LDMC+SLA+rooting_depth+lifespan+clonal+photo_path2+n_fixation+mycorrhizal, data=subset(combined, trt_type2=="n"))

summary(mN)
calc.relimp(mN)

##drought
mdrt<-lm(diff~LDMC+plant_height_vegetative+clonal+lifespan2+n_fixation+mycorrhizal, data=subset(combined, trt_type2=="drought"))

summary(mdrt)
calc.relimp(mdrt)

##irrigation
mirg<-lm(diff~seed_dry_mass+plant_height_vegetative+clonal+lifespan2+photo_path2+mycorrhizal, data=subset(combined, trt_type2=="irrigation"))

summary(mirg)
calc.relimp(mirg)


##temp
mtemp<-lm(diff~rooting_depth+plant_height_vegetative+SLA+lifespan2+clonal+photo_path2+n_fixation, data=subset(combined, trt_type2=="temp"))

summary(mtemp)
calc.relimp(mtemp)

##P
mP<-lm(diff~seed_dry_mass+SLA+clonal+lifespan2+photo_path2+mycorrhizal, data=subset(combined, trt_type2=="p"))

summary(mP)
calc.relimp(mP)


#all int
mAll<-lm(diff~SLA+LDMC+seed_dry_mass+clonal+lifespan2+photo_path2+mycorrhizal+n_fixation, data=subset(combined, trt_type2=="all mult"))

summary(mAll)
calc.relimp(mAll)



####old code here.
#ways to try to get the estimates easily from the categorical models.
anova.lme()

summary(m2_b)
anova(m2)
aov(m2)

contrasts(as.factor(alldat_cat$value))
emmeans(m2, ~ value*treat_type)

#  filter(trt_type2!="co2_other"&trt_type2!="dist_other"&trt_type2!="drt_other"&trt_type2!="herb_rem_other"&trt_type2!="irg_other"&trt_type2!="nuts_other"&trt_type2!="temp_other")%>%
  filter(trait %in% c("lifespan", "clonal", "photosynthetic_pathway", "mycorrhizal", "n_fixation", "growth_form"))

ggplot(data=barGraphStats(data=subset(alldat_cat, ave_diff>0), variable="ave_diff", byFactorNames=c("trait", "trt_type2", "value")), aes(x=value, y=mean))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), width=0.2)+
  facet_grid(trt_type2~trait, scales="free")+
  theme(axis.text.x = element_text(angle=90))



