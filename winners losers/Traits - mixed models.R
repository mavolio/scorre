library(tidyverse)
library(lme4)
library(emmeans)

theme_set(theme_bw(12))

### Trait data

#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  

#read in data
contTraits <- read.csv('C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/Trait Data/TRY Data/Gap_Filled/TRY_new.csv')%>%
  rename(species_matched=Species)%>%
  select(-X.1, -X, -Family, -Genus, -ObservationID)%>%
  group_by(species_matched)%>%
  summarise_all(funs(mean))%>%
  ungroup()

contTraitsSubset <- contTraits%>%
  rename(ssd=X4, rooting_depth=X6, SLA=X11, leaf_C_mass=X13, leaf_N_mass=X14, leaf_P_mass=X15, stem_diameter=X21, seed_mass=X26, seed_length=X27, leaf_thickness=X46, LDMC=X47, leaf_dry_mass=X55, germination_rate=X95, leaf_length=X144, leaf_width=X145, leaf_CN=X146, stem_conduit_density=X169, stem_conduit_diameter=X281, seed_number=X138, SRL=X1080)%>%
  select(-X18, -X50, -X78, -X163, -X223, -X224, -X237, -X282, -X289, -X3112, -X3113, -X3114, -X3120)

traits <- read.csv('C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/CoRRE data\\CoRRE data\\trait data\\sCoRRE categorical trait data - traits_complete_pre spot check_03102021.csv')%>%
  full_join(contTraitsSubset) %>%
  drop_na()%>%
  filter(leaf_P_mass<20, stem_diameter<0.5, seed_mass<50, seed_number<10000, leaf_width<40, stem_conduit_density<1000, stem_conduit_diameter<200)

traitsOutliersRemoved <- traits %>%
  filter(!leaf_type %in% c("microphyll","frond")) %>%
  filter(!species_matched %in% c("Centrolepis aristata", "Centrolepis strigosa", "Acorus calamus"))

traitsScaled <- traitsOutliersRemoved %>% ## only scales continuous traits
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
  filter(trait %in% c("seed_mass", "seed_number", "rooting_depth", "SRL", "LDMC", "SLA"))%>%
  filter(trait=="SRL")

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


mSLA<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=alldat_cont)
summary(mSLA)

plot.sla<-as.data.frame(summary(mSLA)$coefficients)%>%
  mutate(fixedef=row.names(plot))

toplot.SLA<-plot%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="SLA")

mSRL<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=alldat_cont)
summary(mSRL)

plot.srl<-as.data.frame(summary(mSRL)$coefficients)%>%
  mutate(fixedef=row.names(plot))

toplot.SRL<-plot%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="SRL")

toplot<-toplot.SLA%>%
  bind_rows(toplot.SRL)

colnames(toplot)[2] <- "SE"

ggplot(data=toplot, aes(y=Estimate, x=1))+
  geom_point()+
  geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width=0.05)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank())+
  xlab("")+
  scale_x_continuous(limits=c(0, 2))+
  geom_hline(yintercept=0, linetype="dashed")+
  facet_grid(trait~trt_type)



####Categorical data
alldat_cat<-dcidiff_models%>%
  right_join(traitsScaled)%>%
  select(species_matched, trt_type2, site_code, diff, leaf_type, leaf_compoundness, growth_form, photosynthetic_pathway, lifespan, stem_support, clonal,mycorrhizal, mycorrhizal_type, n_fixation, rhizobial,actinorhizal)%>%
  gather(leaf_type:actinorhizal, key="trait", value="value")%>%
  na.omit()%>%
  filter(trait=="photosynthetic_pathway")


m2<-lmer(diff ~ -1 + trt_type2 + value:trt_type2 + (1|species_matched) + (1|site_code), data=alldat_cat)

emmeans(m2, ~ value*trt_type2)

m2_b = lmer(diff ~ -1 + trt_type2 +  (trt_type2|species_matched) + (trt_type2|site_code), data=alldat_cat)


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



