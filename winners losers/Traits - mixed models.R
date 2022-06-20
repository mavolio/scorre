######
#####
##### code to study how species traits affect response to GCDs
##### code by M.Avolio with help from Adam Clark and Tamara Munkemuller
#### created summer 2021, updated with new trait dataset June 20, 2022
######

library(tidyverse)
library(lme4)
library(emmeans)

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

# ##code to drop outliers
# contLong<-contTraits %>% 
#   pivot_longer(LDMC:seed_number, names_to = "trait", values_to = "value")
# ggplot(data=contLong, aes(x=value))+
#   geom_histogram()+
#   facet_wrap(~trait, scales = "free")


catTraits <- read.csv('C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\CoRRE data\\trait data\\Final TRY Traits\\sCoRRE categorical trait data_final_20211209.csv') %>% 
  select(species_matched, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation) %>% 
  filter(growth_form!="moss", mycorrhizal_type!="uncertain", lifespan!="uncertain", clonal!="uncertain", species_matched!="") %>% 
  mutate(mycorrhizal=ifelse(mycorrhizal_type=="none", 'no', ifelse(mycorrhizal_type=="facultative_AM"|mycorrhizal_type=="facultative_AM_EcM", "facultative", "yes"))) %>% 
  select(-mycorrhizal_type) %>% 
  mutate(photo_path=ifelse(photosynthetic_pathway=="possible C4"|photosynthetic_pathway=="possible C4/CAM", "C4", ifelse(photosynthetic_pathway=="possible CAM", "CAM",photosynthetic_pathway))) %>% 
  select(-photosynthetic_pathway)


traitsScaled <- contTraits%>% ## only scales continuous traits
  mutate_at(vars(ssd:SRL), scale)


# Read in dci diff

# dcidiff<-read.csv("C:/Users/megha/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/Species_DCiDiff_newtrts.csv")

dcidiff_models<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/Species_DCiDiff_formixedmodels.csv")

alldat_cont<-dcidiff_models%>%
  right_join(traitsScaled)%>%
  select(-leaf_type, -leaf_compoundness, -growth_form, -photosynthetic_pathway, -lifespan, -stem_support, -clonal, -mycorrhizal, -mycorrhizal_type, -n_fixation, -rhizobial,-actinorhizal)%>%
  filter(seed_mass<4, seed_number<3) %>%
  gather(ssd:SRL, key="trait", value="value")%>%
  na.omit()%>%
  filter(trait %in% c("seed_mass", "seed_number", "rooting_depth", "SRL", "LDMC", "SLA"))

# ggplot(data=alldat_cont, aes(x=value, y=ave_diff))+
#   geom_point()+
#   geom_smooth(method="lm")+
#   facet_grid(trt_type2~trait, scales="free")
####
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

plot.sla<-as.data.frame(summary(mSLA)$coefficients)%>%
  mutate(fixedef=row.names(plot.sla))

toplot.SLA<-plot.sla%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="SLA")

###SRL

mSRL<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="SRL"))
summary(mSRL)

plot.srl<-as.data.frame(summary(mSRL)$coefficients)%>%
  mutate(fixedef=row.names(plot.srl))

toplot.SRL<-plot.srl%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="SRL")

###seed mass

mSM<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="seed_mass"))
summary(mSM)

plot.sm<-as.data.frame(summary(mSM)$coefficients)%>%
  mutate(fixedef=row.names(plot.sm))

toplot.SM<-plot.sm%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="Seed Mass")


###seed number

mSN<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="seed_number"))
summary(mSN)

plot.sn<-as.data.frame(summary(mSN)$coefficients)%>%
  mutate(fixedef=row.names(plot.sn))

toplot.SN<-plot.sn%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="Seed Number")

###rooting depth

mRD<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cont, trait=="rooting_depth"))
summary(mRD)

plot.rd<-as.data.frame(summary(mRD)$coefficients)%>%
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

plot.ldmc<-as.data.frame(summary(mLDMC)$coefficients)%>%
  mutate(fixedef=row.names(plot.ldmc))

toplot.LDMC<-plot.ldmc%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="LDMC")



toplot<-toplot.SLA%>%
  bind_rows(toplot.SRL, toplot.SN, toplot.SM, toplot.RD, toplot.LDMC)

colnames(toplot)[2] <- "SE"

toplotESA<-toplot%>%
  filter(trt_type=="n"|trt_type=="all mult",
         trait=="SLA"|trait=="Rooting Depth"|trait=="Seed Mass")%>%
  mutate(Trt=ifelse(trt_type=="all mult", "Multiple Trts.", "N"))

ggplot(data=toplotESA, aes(y=Estimate, x=1))+
  geom_point()+
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width=0.05)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank())+
  xlab("")+
  scale_x_continuous(limits=c(0, 2))+
  geom_hline(yintercept=0, linetype="dashed")+
  facet_grid(trait~Trt, scales="free")



####Categorical data
alldat_cat<-dcidiff_models%>%
  right_join(traitsScaled)%>%
  select(species_matched, trt_type2, site_code, diff, leaf_type, leaf_compoundness, growth_form, photosynthetic_pathway, lifespan, stem_support, clonal,mycorrhizal, mycorrhizal_type, n_fixation, rhizobial,actinorhizal)%>%
  gather(leaf_type:actinorhizal, key="trait", value="value")%>%
  na.omit()

#photo path
mpp<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="photosynthetic_pathway"))
summary(mpp)

plot.mpp<-as.data.frame(emmeans(mpp, ~ value*trt_type2))%>%
  filter(value!="CAM")

#lifespan
ml<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="lifespan"))
summary(ml)

plot.ml<-as.data.frame(emmeans(ml, ~ value*trt_type2))%>%
  filter(value!="biennial")

###clonaltiy
mc<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="clonal"))

plot.mc<-as.data.frame(emmeans(mc, ~ value*trt_type2))%>%
  mutate(trait="Clonal")%>%
  rename(v=value)%>%
  mutate(value=paste(trait, v, sep=" "))

##nfixer
mnf<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="n_fixation"))

plot.nf<-as.data.frame(emmeans(mnf, ~ value*trt_type2))%>%
  mutate(trait="Nfix")%>%
  rename(v=value)%>%
  mutate(value=paste(trait, v, sep=" "))


##mycorhizal
mmyc<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=subset(alldat_cat, trait=="mycorrhizal"))

plot.myc<-as.data.frame(emmeans(mmyc, ~ value*trt_type2))%>%
  mutate(trait="Myc")%>%
  rename(v=value)%>%
  mutate(value=paste(trait, v, sep=" "))%>%
  filter(v!="", v!="uncertain")


toplot.cat<-plot.mpp%>%
  bind_rows(plot.ml, plot.mc, plot.nf, plot.myc)

toplotesacat<-toplot.cat%>%
  filter(trt_type2=="n"|trt_type2=="all mult",
         value=="C3"|value=="C4"|value=="annual"|value=="Nfix yes")%>%
  mutate(Trt=ifelse(trt_type2=="all mult", "Multiple Trts.", "N"),
         trait=ifelse(value=="C3", "C3", ifelse(value=="C4", "C4", ifelse(value=="annual", "Annual", "N-fixer"))))%>%
  rename(Estimate=emmean)%>%
  select(trait, Trt, Estimate, SE)

ggplot(data=toplotesacat, aes(y=emmean, x=1))+
  geom_point()+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.05)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank())+
  xlab("")+
  scale_x_continuous(limits=c(0, 2))+
  geom_hline(yintercept=0, linetype="dashed")+
  facet_grid(trait~Trt, scales="free")


ESA<-toplotESA%>%
  select(trait, Trt, Estimate, SE)%>%
  bind_rows(toplotesacat)


ggplot(data=ESA, aes(y=Estimate, x=1))+
  geom_point()+
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width=0.05)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank())+
  xlab("")+
  scale_x_continuous(limits=c(0, 2))+
  geom_hline(yintercept=0, linetype="dashed")+
  facet_grid(Trt~trait, scales="free")


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



