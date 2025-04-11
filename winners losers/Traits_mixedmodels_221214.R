######
#####
##### code to study how species traits affect response to GCDs
##### code by M.Avolio with help from Adam Clark and Tamara Munkemuller
#### created summer 2021, updated with new trait dataset June 20, 2022
#### updated 14.12.2022 by A. Clark
######
rm(list=ls())

library(tidyverse)
library(lme4)
library(emmeans)
require(RColorBrewer)
#library(relaimpo)

scale_fun = function(x) {
  x.new = log(x)
  x.new = x-mean(x,na.rm=T)
  x.new = x/sd(x, na.rm=T)
  x.new
}

#open TraitBoxPlots code and run the function make_boxplots

setwd('C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/manuscript')

# Read in dci diff
dcidiff_models<-read.csv("C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/WinnersLosers paper/data/Species_DCiDiff_formixedmodelsMarch2024.csv") %>% 
  rename(species=species_matched)

length(unique(dcidiff_models$species))


# Continuous trait data ---------------------------------------------------


#read in data from EDI
inUrl2  <- "https://pasta.lternet.edu/package/data/eml/edi/1533/1/169fc12d10ac20b0e504f8d5ca0b8ee8" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl"))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")


conttraits1<-read.csv(infile2,header=F 
                      ,skip=1
                      ,sep=","  
                      ,quot='"' 
                      , col.names=c(
                        "family",     
                        "species",     
                        "trait",     
                        "trait_value",     
                        "error_risk_overall",     
                        "error_risk_family",     
                        "error_risk_genus",     
                        "source"    ), check.names=TRUE)

 unlink(infile2)

#error risk < 3 drops ~28 observations
#error risk < 2.5 drops ~99 observations
#error risk < 2 drop 328 obs
remove_risk<-conttraits1 %>% 
  filter(error_risk_overall<2|is.na(error_risk_overall))

contTraits<-remove_risk%>%
  filter(trait %in% c('LDMC', 'SLA', 'plant_height_vegetative', 'SRL', 'leaf_N', 'seed_dry_mass')) %>%   select(-error_risk_overall, -error_risk_family, -error_risk_genus) %>%
  group_by(species, trait)%>%
  summarise(value=mean(trait_value, na.rm=T))

###looking for outliers
ggplot(data=contTraits, aes(x=value))+
  geom_histogram()+
  facet_wrap(~trait, scales = "free")

##merge trait data with species responses
alldat_cont<-dcidiff_models%>%
  left_join(contTraits, by="species") %>% 
  drop_na()

###how many species included in the continuous analyses?
length(unique(alldat_cont$species))

#summarizing the data
summaryCont<-alldat_cont %>% 
  select(species, trait, value) %>% 
  unique() %>% 
  drop_na() %>% 
  group_by(trait) %>% 
  summarize(min=min(value), max=max(value), n=length(value))

#Making a graph of all the models I am goign to run now, just to see if there are patterns.
ggplot(data=alldat_cont, aes(x=scale_fun(value), y=diff))+
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

#getting number of species per trait
tmp = alldat_cont
tmp = tmp[!is.na(tmp$diff),]
cont_n_data = table(unique(tmp[,c("species", "trt_type2", "trait")]))
cont_n_data = apply(cont_n_data, 2:3, sum)
#write.csv(cont_n_data, "cont_n_data.csv")


# Continuous models -------------------------------------------------------

##SLA
mSLA<-lmer(diff ~ -1 + trt_type2 + scale_fun(value):trt_type2 + (1|species) + (1|site_code), data=subset(alldat_cont, trait=="SLA"))
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

###SRL

mSRL<-lmer(diff ~ -1 + trt_type2 + scale_fun(value):trt_type2 + (1|species) + (1|site_code), data=subset(alldat_cont, trait=="SRL"))
summary(mSRL)

plot.srl<-as.data.frame(summary(mSRL)$coefficients)

plot.srl<-plot.srl%>%
  mutate(fixedef=row.names(plot.srl))

toplot.SRL<-plot.srl%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="SRL")

###seed mass
mSM<-lmer(diff ~ -1 + trt_type2 + scale_fun(value):trt_type2 + (1|species) + (1|site_code), data=subset(alldat_cont, trait=="seed_dry_mass"))
summary(mSM)

plot.sm<-as.data.frame(summary(mSM)$coefficients)

plot.sm<-plot.sm %>% mutate(fixedef=row.names(plot.sm))

toplot.SM<-plot.sm%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="Seed Mass")

###LDMC

mLDMC<-lmer(diff ~ -1 + trt_type2 + scale_fun(value):trt_type2 + (1|species) + (1|site_code), data=subset(alldat_cont, trait=="LDMC"))
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

mhght<-lmer(diff ~ -1 + trt_type2 + scale_fun(value):trt_type2 + (1|species) + (1|site_code), data=subset(alldat_cont, trait=="plant_height_vegetative"))
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


##LeafN
mN<-lmer(diff ~ -1 + trt_type2 + scale_fun(value):trt_type2 + (1|species) + (1|site_code), data=subset(alldat_cont, trait=="leaf_N"))
summary(mN)

plot.leafN<-as.data.frame(summary(mN)$coefficients)

plot.leafN<-plot.leafN%>%
  mutate(fixedef=row.names(plot.leafN))

toplot.leafN<-plot.leafN%>%
  separate(fixedef, into=c("trt_type2", "interaction"), sep=":")%>%
  filter(!is.na(interaction))%>%
  separate(trt_type2, into=c("drop", "trt_type"), sep=9)%>%
  select(-interaction, -drop)%>%
  mutate(trait="LeafN")


toplot<-toplot.SLA%>%
  bind_rows(toplot.SM, toplot.SRL, toplot.LDMC, toplot.hgt, toplot.leafN) %>% 
  rename(SE="Std. Error") %>% 
  mutate(trt_type2=factor(trt_type, levels=c("co2", "drought", "irrigation", "temp", "n", "p", "multnuts", "all mult"))) %>% 
  mutate(min=Estimate-SE, max=Estimate+SE)

# Continuous plot ---------------------------------------------------------

trt.labels=c(co2="CO2", drought="Drt", irrigation="Irg.", temp="Temp.", n="N", p="P", multnuts="Mult. Nut.","all mult" ="Interact.")
trait.labels=c(LDMC="LDMC", LeafN="Leaf N",
               PlantHeight="Plant Height", SRL="SRL",
               "Seed Mass"="Seed Mass", SLA="SLA")


pdf("traits_by_treat_contMarch2024_2.pdf", width = 5.2, height=10)
make_boxplot(toplot_data = toplot,
                        trt.labels_data = trt.labels,
                        trait.labels_data=trait.labels,
                        groupbytrait = FALSE,
                        legend_line_length = 1.4,
                        legend_textwidth = 0.012,
                        legend_yadj = 0.175,
                        legend_xadj = -0.002,
                        p_alpha = 0.05,
                        lower_margin = 6.8,
                        sigadj = -0.03,
                        traitorder = rev(c(2,1,6,4,3,5)),#number refers to position in trt label vector, do in reverse
                        group_colors = adjustcolor(rev(c("darkgreen",
                                                         "darkgreen",
                                                         "darkgreen",
                                                         "blue",
                                                         "blue",
                                                         "orange")), alpha.f = 0.08))
dev.off()



# Categorical data --------------------------------------------------------

#read in data from EDI
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/1533/1/5ebbc389897a6a65dd0865094a8d0ffd" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

cattraits1 <-read.csv(infile1,header=F 
                      ,skip=1
                      ,sep=","  
                      ,quot='"' 
                      , col.names=c(
                        "family",     
                        "species",     
                        "trait",     
                        "trait_value",     
                        "source",     
                        "error_risk_overall"    ), check.names=TRUE)

unlink(infile1)


#simplify groupings and select categories
catTraits <- cattraits1 %>% 
  filter(trait %in% c('clonal', 'growth_form', 'lifespan', 'mycorrhizal_type', 'n_fixation_type', 'photosynthetic_pathway')) %>% 
  mutate(trait_value2=ifelse(trait=='photosynthetic_pathway' & trait_value %in% c('C4', 'CAM'), 'C4/CAM', 
                      ifelse(trait=='n_fixation_type' & trait_value %in% c('actinorhizal', 'rhizobial'), 'N-fixer',
                      ifelse(trait=='mycorrhizal_type' & trait_value %in% c('AM', 'EcM', 'ErM', 'OM', 'multiple'), 'yes',
                      ifelse(trait=='lifespan' & trait_value %in% c('perennial', 'biennial'), 'Perenn./Bienn.', trait_value))))) %>% 
  mutate(drop=ifelse(trait=='photosynthetic_pathway' & trait_value2 %in% c('C3', 'C4/CAM'), 0, 
              ifelse(trait=='n_fixation_type' & trait_value2 %in% c('N-fixer', 'none'), 0,
              ifelse(trait=='mycorrhizal_type' & trait_value2 %in% c('yes', 'none'), 0, 
              ifelse(trait=='clonal' & trait_value2 %in% c('yes', 'no'), 0,
              ifelse(trait=='growth_form' & trait_value2 %in% c('forb', 'graminoid', 'woody'), 0,
              ifelse(trait=='lifespan' & trait_value2 %in% c('Perenn./Bienn.', 'annual'), 0, 1))))))) %>% 
  filter(drop==0)



alldat_cat<-dcidiff_models%>%
  left_join(catTraits) %>% 
  select(-trait_value, -source, -error_risk_overall) %>% 
  drop_na() #this drops one species that is a hybrid and had no trait data

length(unique(alldat_cat$species))

cattrait_values<-alldat_cat %>% 
  select(species, trait, trait_value2) %>% 
  unique() %>% 
  group_by(trait, trait_value2) %>% 
  summarize(n=length(species))

###This give the number of obs for the figure
tmp = alldat_cat
tmp = tmp[!is.na(tmp$diff),]
tmp$trait_value = paste(tmp$trait, tmp$trait_value2, sep = "_")
cat_n_data = table(unique(tmp[,c("species", "trt_type2", "trait_value")]))
cat_n_data = apply(cat_n_data, 2:3, sum)

value2 = colnames(cat_n_data)
tmp=ifelse(value2=='clonal_no', "Non-Clonal",
  ifelse(value2=='clonal_yes', 'Clonal',
  ifelse(value2=='lifespan_annual', 'Annual',
  ifelse(value2=='lifespan_Perenn./Bienn.', 'Perenn./Bienn.',
  ifelse(value2=='mycorrhizal_type_yes', 'Mycorr.',
  ifelse(value2=='mycorrhizal_type_none', 'Non-Mycorr.',
  ifelse(value2=='n_fixation_type_N-fixer', 'N-Fixer',
  ifelse(value2=='n_fixation_type_none', 'Non-N-Fixer',
  ifelse(value2=='photosynthetic_pathway_C3', "C3",
  ifelse(value2=='photosynthetic_pathway_C4/CAM', 'C4/CAM', value2))))))))))
tmp = gsub("growth_form_", "", tmp, fixed = TRUE)
#tmp = paste(toupper(substr(tmp,1,1)), substr(tmp,2,99), sep = "")
colnames(cat_n_data) = tmp
#write.csv(cat_n_data, "cat_n_data.csv")


# Categorical trait models ------------------------------------------------

##overall response of species
overall<-lmer(diff ~ -1 + trt_type2 + (1|species) + (1|site_code), data=dcidiff_models)
summary(overall)

plot.overall<-as.data.frame(emmeans(overall, ~ trt_type2)) %>% 
  mutate(trait='overall')


#photo path
mpp<-lmer(diff ~ -1 + trt_type2 + trait_value2:trt_type2 + (1|species) + (1|site_code), data=subset(alldat_cat, trait=="photosynthetic_pathway"))
summary(mpp)

plot.mpp<-as.data.frame(emmeans(mpp, ~ trait_value2*trt_type2)) %>% 
  mutate(trait="photosynthetic_pathway")

#lifespan
ml<-lmer(diff ~ -1 + trt_type2 + trait_value2:trt_type2 + (1|species) + (1|site_code), data=subset(alldat_cat, trait=="lifespan"))
summary(ml)

plot.ml<-as.data.frame(emmeans(ml, ~ trait_value2*trt_type2)) %>% 
  mutate(trait="lifespan")

###clonaltiy
mc<-lmer(diff ~ -1 + trt_type2 + trait_value2:trt_type2 + (1|species) + (1|site_code), data=subset(alldat_cat, trait=="clonal"))

plot.mc<-as.data.frame(emmeans(mc, ~ trait_value2*trt_type2))%>%
  mutate(trait="Clonal")

##nfixer
mnf<-lmer(diff ~ -1 + trt_type2 + trait_value2:trt_type2 + (1|species) + (1|site_code), data=subset(alldat_cat, trait=="n_fixation_type"))

plot.nf<-as.data.frame(emmeans(mnf, ~ trait_value2*trt_type2))%>%
  mutate(trait="n_fixation_type")


##mycorhizal
mmyc<-lmer(diff ~ -1 + trt_type2 + trait_value2:trt_type2 + (1|species) + (1|site_code), data=subset(alldat_cat, trait=="mycorrhizal_type"))

plot.myc<-as.data.frame(emmeans(mmyc, ~ trait_value2*trt_type2))%>%
  mutate(trait="mycorrhizal_type")

##growthform
mgf<-lmer(diff ~ -1 + trt_type2 + trait_value2:trt_type2 + (1|species) + (1|site_code), data=subset(alldat_cat, trait=="growth_form"))
plot.mgf<-as.data.frame(emmeans(mgf, ~ trait_value2*trt_type2))%>%
  mutate(trait="GF")

# Graphing categorical traits ---------------------------------------------



toplot.cat<-plot.mpp%>%
  bind_rows(plot.ml, plot.mc, plot.nf, plot.myc, plot.mgf, plot.overall)


toplotesacat<-toplot.cat %>% 
  mutate(trt_type3=factor(trt_type2, levels=c("co2", "drought", "irrigation", "temp", "n", "p", "multnuts", "all mult"))) %>% 
  mutate(value2=paste(trait, trait_value2, sep="_" )) %>% 
  mutate(min=emmean-SE, max=emmean+SE) %>% 
  mutate(sig=ifelse(min>0&max>0, "*", ifelse(min<0&max<0, "*", ""))) %>% 
  mutate(Trait_name=ifelse(value2=='Clonal_no', "Non-Clonal",
                    ifelse(value2=='Clonal_yes', 'Clonal',
                    ifelse(value2=='lifespan_annual', 'Annual',
                    ifelse(value2=='lifespan_Perenn./Bienn.', 'Perenn./Bienn.',
                    ifelse(value2=='mycorrhizal_type_yes', 'Mycorr.',
                    ifelse(value2=='mycorrhizal_type_none', 'Non-Mycorr.',
                    ifelse(value2=='n_fixation_type_N-fixer', 'N-Fixer',
                    ifelse(value2=='n_fixation_type_none', 'Non-N-Fixer',
                    ifelse(value2=='photosynthetic_pathway_C3', "C3",
                    ifelse(value2=='photosynthetic_pathway_C4/CAM', 'C4/CAM', 
                    ifelse(value2=='overall_NA', 'Overall', "GF"))))))))))))# %>%
  #mutate(Traits2=factor(Trait_name,levels=c("Annual", 'Perennial', 'Clonal', 'Non-Clonal', 'C3', 'C4', 'Mycorr.', 'Non-Mycorr.', 'N-Fixer', 'LDMC', 'Plant Height', 'Rooting Depth', 'Seed Mass', 'SLA'))) #%>% 
  #na.omit()
head(toplotesacat)
toplotesacat$value = as.character(toplotesacat$trait_value2)
toplotesacat$Trait_name = as.character(toplotesacat$Trait_name)
toplotesacat$value[toplotesacat$Trait_name!="GF"] = toplotesacat$Trait_name[toplotesacat$Trait_name!="GF"]
toplotesacat$value = factor(toplotesacat$value)

unique(data.frame(toplotesacat$value, toplotesacat$Trait_name))


trt.labels=c(co2="CO2", drought="Drt", irrigation="Irg", temp="Temp", n="N", p="P", multnuts="Mult. Nut.","all mult" ="Interact.")
trait.labels = sort(unique(as.character(toplotesacat$value)))
tmp = trait.labels
trait.labels = paste(toupper(substr(trait.labels,1,1)),
      substr(trait.labels,2,99), sep = "")
names(trait.labels) = tmp
trait.categories = c("Lifespan", "Photo. Pathway", "Photo. Pathway",
                     "Clonality", "Growth Form",
                     "Growth Form", "Mycorr. Assoc.",
                     "N. Fixation", "Clonality", "Mycorr. Assoc.",
                     "N. Fixation", "Overall", "Lifespan", "Growth Form")
trait.super.categories = c("Reproduction", "Leaf Traits", "Leaf Traits",
                     "Reproduction", "Growth Form",
                     "Growth Form", "Symbiosis",
                     "Symbiosis", "Reproduction", "Symbiosis",
                     "Symbiosis", "Overall", "Reproduction", "Growth Form")
cbind(trait.labels, trait.categories, trait.super.categories)
toplotesacat$Estimate = toplotesacat$emmean
toplotesacat$trait = toplotesacat$value

##all cat in one fig

#ggplot(data=toplotesacat, aes(y=emmean, x=1, label=sig))+
#  geom_point()+
#  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.05)+
#  scale_y_continuous(breaks = c(-0.05, 0, 0.05))+
#  coord_flip()+
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank())+
#  xlab("")+
#  ylab("Estimate")+
#  scale_x_continuous(limits=c(0, 2))+
#  geom_hline(yintercept=0, linetype="dashed")+
#  geom_text(nudge_x = 0.2, nudge_y = 0.01, size=10, color="red")+
#  facet_grid(trt_type3~Traits2, labeller = labeller(trt_type3=trt.labels))

##############
# new plot
##############
gcol = adjustcolor(rev(c(rep("gray",1),
                         rep("red", 3),
                         rep("darkgreen",2),
                         rep("purple",4),
                         rep("orange",4))), alpha.f = 0.08)
# tord = rev(c(2,3,
#              7,10,
#              8,11,
#              1,12,
#              4,9,
#              5,6,
#              13,14))

tord = rev(c(12, 6,
             5, 14,
             2,3,
             7,10,
             8,11,
             4,9,
             13,1))

if(FALSE) {
pdf("traits_by_treat_catMarch2024.pdf", width = 5.2, height=10)
make_boxplot(toplot_data = toplotesacat,
                        trt.labels_data = trt.labels,
                        trait.labels_data=trait.labels,
                        groupbytrait = FALSE,
                        legend_line_length = 1.4,
                        legend_textwidth = 0.037,
                        legend_cex = 0.9,
                        legend_yadj = 0.1,
                        legend_xadj = -0.01,
                        p_alpha = 0.05,
                        lower_margin = 6.8,
                        sigadj = -0.09,
                        traitorder = tord,
                        group_colors = gcol,
             xtit = "Effect Size, Mean DCi diff. by. Trait",
             n_table = cat_n_data)
dev.off()
}

# tord1 = rev(c(2,11,7,12,4,5,14))
# tord2 = rev(c(3,8,10,1,9,6,13))
tord1 = rev(c(12,5,2,11,7,4,13))
tord2 = rev(c(6,14,3,8,10,9,1))
gcol_split1 = adjustcolor(rev(c(rep("gray", 1),
                               rep("red", 1),
                               rep("darkgreen",1),
                         rep("purple",2),
                         rep("orange",2))), alpha.f = 0.08)

gcol_split2 = adjustcolor(rev(c(rep("red", 2),
                                rep("darkgreen",1),
                                rep("purple",2),
                                rep("orange",2))), alpha.f = 0.08)


pdf("traits_by_treat_cat_2colMarch2024.pdf", width = 10.4, height=10)
par(mar=c(2,6.8,3.5,0.2), oma =c(3,1,0,0), mfrow=c(1,2))#controlling margins of plots
make_boxplot(toplot_data = toplotesacat,
             trt.labels_data = trt.labels,
             trait.labels_data=trait.labels,
             groupbytrait = FALSE,
             legend_line_length = 1.4,
             legend_textwidth = 0.037,
             legend_cex = 1.2,
             legend_yadj = 0.185,
             legend_xadj = 0.14,
             p_alpha = 0.05,
             sigadj = -0.04,
             xlm = c(-0.12,0.12),
             traitorder = tord1,
             group_colors = gcol_split1, 
             xtit = "", autopar = FALSE, n_table = cat_n_data)

par(mar=c(2,0.2,3.5,6.8))
make_boxplot(toplot_data = toplotesacat,
             trt.labels_data = trt.labels,
             trait.labels_data=trait.labels,
             groupbytrait = FALSE,
             legend_line_length = 1.4,
             legend_textwidth = 0.037,
             legend_cex = 0.9,
             legend_yadj = 10,
             legend_xadj = -0.01,
             p_alpha = 0.05,
             sigadj = -0.04,
             xlm = c(-0.12,0.12),
             traitorder = tord2,
             group_colors = gcol_split2,
             xtit = "", axisside = 4, autopar = FALSE, n_table = cat_n_data)
mtext("Effect Size, Mean DCi diff. by. Trait", side = 1, outer = TRUE, line = 1.2, cex = 1.7)
dev.off()
dev.off()



# TraitSyndromes ----------------------------------------------------------



###thinking about trait categories
catTraits_full <- cattraits1 %>% 
  filter(trait %in% c('clonal', 'growth_form', 'lifespan', 'mycorrhizal_type', 'n_fixation_type', 'photosynthetic_pathway')) %>% 
  mutate(trait_value2=ifelse(trait=='photosynthetic_pathway' & trait_value %in% c('C4', 'CAM'), 'C4/CAM',    ifelse(trait=='n_fixation_type' & trait_value %in% c('actinorhizal', 'rhizobial'), 'N-fixer',ifelse(trait=='mycorrhizal_type' & trait_value %in% c('AM', 'EcM', 'ErM', 'OM', 'multiple'), 'yes',ifelse(trait=='lifespan' & trait_value %in% c('perennial', 'biennial'), 'Perenn./Bienn.', trait_value))))) 

alldat_cat_full<-dcidiff_models%>%
  left_join(catTraits_full) %>% 
  select(-trait_value, -source, -error_risk_overall) %>% 
  drop_na()

traitcats_full<-alldat_cat_full %>%
  pivot_wider(names_from=trait, values_from = trait_value2, values_fill = NA) %>% 
  drop_na() %>% 
  select(species, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation_type) %>% 
  unique() %>% 
  group_by(growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation_type) %>% 
  summarize(n=length(species))

sum(traitcats_full$n)

TraitSyndrome<-alldat_cat_full %>%
  pivot_wider(names_from=trait, values_from = trait_value2, values_fill = NA) %>% 
  select(species, growth_form, photosynthetic_pathway, lifespan, clonal, mycorrhizal_type, n_fixation_type) %>% 
  unique() %>% 
 mutate(syndrome=
          ifelse(photosynthetic_pathway=="C4/CAM", 'C4/CAM', 
          ifelse(growth_form=='forb'&lifespan=='annual', 'Annual forb',
          ifelse(growth_form=='graminoid'&lifespan=='annual', 'Annual gram.', 
          ifelse(growth_form=='forb'&lifespan=='Perenn./Bienn.'&mycorrhizal_type=='none'&n_fixation_type=='none', 'Non-mutualistic peren. forb', 
          ifelse(growth_form=='forb'&lifespan=='Perenn./Bienn.'&n_fixation_type=='N-fixer', 'N-fixing peren. forb', 
          ifelse(growth_form=='forb'&lifespan=='Perenn./Bienn.'&mycorrhizal_type=='yes','Mycorrhizal forb', 
          #ifelse(growth_form=='forb'&lifespan=='Perenn./Bienn.'&clonal=='no'&mycorrhizal_type=='yes', 'ForbPernNoClonMyc',
          ifelse(growth_form=='woody', 'Woody', 
          ifelse(growth_form=='graminoid'&lifespan=='Perenn./Bienn.', 'Peren. Gram.', 'todo')))))))))

SyndromeN<-TraitSyndrome %>% 
  group_by(syndrome) %>% 
  summarise(n=length(species))

SpSyndrome<-TraitSyndrome %>% 
  select(species, syndrome) %>% 
  filter(syndrome!='todo')

syndrome_diff_all <-dcidiff_models %>% 
  left_join(SpSyndrome)


###This give the number of obs for the figure
tmp = syndrome_diff_all
tmp = tmp[!is.na(tmp$diff),]
synd_n_data = table(unique(tmp[,c("species", "trt_type2", "syndrome")]))
synd_n_data = apply(synd_n_data, 2:3, sum)

value2 = colnames(synd_n_data)

#syndromes
msynd<-lmer(diff ~ -1 + trt_type2 + syndrome:trt_type2 + (1|species) + (1|site_code), data=syndrome_diff_all)

plot.msynd<-as.data.frame(emmeans(msynd, ~ syndrome*trt_type2))%>%
  mutate(trait="Synd") %>% 
  mutate(sig=ifelse(asymp.LCL<0&asymp.UCL<0|asymp.LCL>0&asymp.UCL>0, 1, 0))

ggplot(data=plot.msynd, aes(x=emmean, y=trt_type2, color=as.factor(sig)))+
  geom_point()+
  geom_errorbar(aes(xmin=emmean-1.96*SE, xmax=emmean+1.96*SE), width=0.05)+
  geom_vline(xintercept=0)+
  facet_wrap(~syndrome)

###modifying text for trait syndromes
toplot.synd<-plot.msynd

toplotsynd<-toplot.synd %>% 
  mutate(trt_type3=factor(trt_type2, levels=c("co2", "drought", "irrigation", "temp", "n", "p", "multnuts", "all mult"))) %>% 
  mutate(Trait_name=syndrome)# %>% 
#mutate(Traits2=factor(Trait_name,levels=c("Annual", 'Perennial', 'Clonal', 'Non-Clonal', 'C3', 'C4', 'Mycorr.', 'Non-Mycorr.', 'N-Fixer', 'LDMC', 'Plant Height', 'Rooting Depth', 'Seed Mass', 'SLA'))) #%>% 
#na.omit()
head(toplotsynd)
toplotsynd$value = as.character(toplotsynd$syndrome)
toplotsynd$Trait_name = as.character(toplotsynd$Trait_name)
toplotsynd$value = factor(toplotsynd$value)

unique(data.frame(toplotsynd$value, toplotsynd$Trait_name))


trt.labels=c(co2="CO2", drought="Drt", irrigation="Irg", temp="Temp", n="N", p="P", multnuts="Mult. Nut.","all mult" ="Interact.")
trait.labels = sort(unique(as.character(toplotsynd$value)))
tmp = trait.labels
trait.labels = paste(toupper(substr(trait.labels,1,1)),
                     substr(trait.labels,2,99), sep = "")
names(trait.labels) = tmp

toplotsynd$Estimate = toplotsynd$emmean
toplotsynd$trait = toplotsynd$value


pdf("syndromes_by_treat_contMarch2024.pdf", width = 6, height=10)
make_boxplot(toplot_data = toplotsynd,
             trt.labels_data = trt.labels,
             trait.labels_data=trait.labels,
             groupbytrait = FALSE,
             legend_line_length = 1.4,
             legend_textwidth = 0.02,
             legend_yadj = 0.175,
             legend_xadj = -0.002,
             p_alpha = 0.05,
             lower_margin = 6.8,
             sigadj = -0.03,
             traitorder = rev(c(1,2,3,4,5, 6, 7, 8)),#number refers to position in trt label vector, do in reverse
             #group_colors = adjustcolor(rev(c('yellow', 'pink', 'pink', 'pink','pink', 'green', 'green', 'brown')), alpha.f = 0.08),
             n_table = synd_n_data)
dev.off()

                 