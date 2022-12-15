
#investigate winners losers

library(tidyverse)
library(gridExtra)

###read in data

my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "E:/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Working groups\\sDiv\\"
  
#read in the data

#raw abundance data and drop pretreatment years
dat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Dec2021.csv",sep="")) %>% 
  filter(treatment_year!=0)

sp <-read.csv(paste(my.wd,"CoRRE data/trait data/corre2trykey_2021.csv", sep=""))%>%
  ungroup() %>% 
  select(genus_species, species_matched)%>%
  unique()

#clean species names and sum over them, drop non-vascular plants, drop plants not IDd to species

dat_cleansp<-dat%>%
  left_join(sp) %>% 
  na.omit()%>% ###get's rid of plants no identified to species
  mutate(drop=ifelse(species_matched %in% c("Andreaea obovata", "Anthelia juratzkana" , "Aulacomnium turgidum", "Barbilophozia hatcheri", "Barbilophozia kunzeana" , "Blepharostoma trichophyllum", "Brachythecium albicans", "Bryum arcticum", "Bryum pseudotriquetrum",  "Campylium stellatum",         "Cyrtomnium hymenophyllum" ,   "Dicranoweisia crispula",  "Dicranum brevifolium", "Dicranum elongatum",  "Dicranum fuscescens", "Dicranum groenlandicum",  "Dicranum scoparium" , "Distichium capillaceum",  "Ditrichum flexicaule",        "Gymnomitrion concinnatum" ,   "Hamatocaulis vernicosus",   "Homalothecium pinnatifidum",  "Hylocomium splendens",        "Hypnum cupressiforme", "Hypnum hamulosum", "Isopterygiopsis pulchella",   "Kiaeria starkei", "Leiocolea heterocolpos",      "Marchantia polymorpha",       "Marsupella brevissima",  "Meesia uliginosa", "Myurella tenerrima", "Oncophorus virens",  "Oncophorus wahlenbergii", "Pleurozium schreberi","Pogonatum urnigerum" ,  "Pohlia cruda" , "Pohlia nutans","Polytrichastrum alpinum", "Polytrichum juniperinum",     "Polytrichum piliferum",  "Polytrichum strictum", "Preissia quadrata", "Ptilidium ciliare",           "Racomitrium lanuginosum", "Rhytidium rugosum", "Saelania glaucescens", "Sanionia uncinata",   "Schistidium apocarpum", "Syntrichia ruralis","Tomentypnum nitens", "Tortella tortuosa",           "Tritomaria quinquedentata",   "Nephroma arcticum" , "Unknown NA", "Campylopus flexuosus",        "Hypnum jutlandicum","Plagiothecium undulatum",   "Polytrichum commune","Pseudoscleropodium purum",    "Rhytidiadelphus loreus",  "Rhytidiadelphus triquetrus",  "Thuidium tamariscinum"), 1, 0)) %>% 
  filter(drop==0)%>% #drops mosses
  group_by(site_code, project_name, community_type, calendar_year, treatment_year, treatment, plot_id, species_matched)%>% #sums plants that were counted twice in the same plot but with different names that once cleaned were combined
  summarize(relcov=sum(relcov))


#info on treatments
trts<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_ExperimentInfo_Dec2021.csv", sep=""))%>%
  select(site_code, project_name, community_type, treatment, trt_type, pulse, plot_mani,resource_mani)%>%
  unique()

#getting a list of all the treatments I want to focus on
trt_analysis<-trts%>%
  mutate(alltrts=ifelse(trt_type %in% c("CO2","CO2*temp", "drought","drought*CO2*temp","irr*CO2","irr*CO2*temp","N*CO2*temp","N*irr*CO2", "mult_nutrient*irr","N*irr*CO2*temp", "N","mult_nutrient","N*P","P","N*CO2","N*drought","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","temp","drought*temp","irr*temp","irr"),1,0))%>%
  filter(alltrts==1)%>%
  mutate(CO2=ifelse(trt_type %in% c("CO2"), 1, 0),
         drought=ifelse(trt_type %in% c("drought"), 1, 0),
         irg=ifelse(trt_type %in% c("irr"), 1, 0),
         temp=ifelse(trt_type %in% c("temp"), 1, 0),
         n=ifelse(trt_type=="N", 1, 0),
         p=ifelse(trt_type=="P", 1, 0),
         multnuts=ifelse(trt_type %in% c("mult_nutrient","N*P"), 1, 0),
         multtrts=ifelse(trt_type %in% c("CO2*temp", "drought*CO2*temp","irr*CO2","irr*CO2*temp","N*CO2*temp","N*irr*CO2", "mult_nutrient*irr","N*irr*CO2*temp", "N*CO2","N*drought","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","drought*temp","irr*temp"),1,0))


##Getting DCI

#combine relative abundance data with treatment na.omit removes unidentified species 
alldat<-dat_cleansp%>%
  left_join(trts)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

## Need to add in zeros

spc<-unique(alldat$site_project_comm)

datfilled<-data.frame()

for (i in 1:length(spc)){
  
  subset<-alldat%>%
    filter(site_project_comm==spc[i])
  
  addzero<-subset %>% 
    select(site_code, project_name, community_type, calendar_year, treatment_year, treatment, plot_id, trt_type, plot_mani, site_project_comm, species_matched, relcov)%>%
    pivot_wider(names_from = "species_matched", values_from = "relcov", values_fill=0)
  
  datfill<-addzero %>% 
    pivot_longer(11:ncol(addzero), names_to = "species_matched", values_to = "relcov")
    
  datfilled<-datfilled%>%
    bind_rows(datfill)
}

###get relative cover averaged each species in a treatment, over all years of the experiment over all plots.

relave<-datfilled%>%
  group_by(site_code, project_name, community_type, species_matched, plot_mani, treatment, trt_type) %>% 
  summarize(mean=mean(relcov))

#subset out control plots, plot_mani==0
Crelave<-relave%>%
  filter(plot_mani==0)%>%
  ungroup()%>%
  select(-treatment, -plot_mani)

#subset out treated plots, plot_mani!=0
Trelave<-relave%>%
  filter(plot_mani!=0)%>%
  ungroup()%>%
  select(-plot_mani)

#to get relative frequency, determine number of control plots
controlplots<-alldat%>%
  filter(plot_mani==0)%>%
  ungroup()%>%
  select(site_code, project_name, community_type, plot_id)%>%
  unique()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(ncplots=length(plot_id))

#to get relative frequency, determine number of control plots a species is found in, merge in total number of plots and calculate relative frequency. By taking this from all dat there are no options for zero and thus all using n to calcualte length will work.
control_freq<-alldat%>%
  filter(plot_mani==0)%>%
  ungroup() %>% 
  select(site_code, project_name, community_type, species_matched, plot_id)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, species_matched)%>%
  summarize(nplots=length(plot_id))%>%
  left_join(controlplots)%>%
  mutate(freq=nplots/ncplots)

#calculate DCi
control_dom<-Crelave%>%
  left_join(control_freq)%>%
  mutate(freq=replace_na(freq, 0)) %>% 
  mutate(DCi=(mean+freq)/2)%>%
  select(site_code, project_name, community_type, species_matched, mean, freq, DCi)%>%
  rename(cmean=mean, cfreq=freq)

#getting frequency of treated plots, same code as above but for treated plots
treatplots<-alldat%>%
  filter(plot_mani!=0)%>%
  ungroup() %>% 
  select(site_code, project_name, community_type, treatment, plot_id)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarize(ntplots=length(plot_id))

treat_freq<-alldat%>%
  filter(plot_mani!=0)%>%
  ungroup()%>%
  select(site_code, project_name, community_type, treatment, species_matched, plot_id)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment, species_matched)%>%
  summarize(nplots=length(plot_id))%>%
  left_join(treatplots)%>%
  mutate(freq=nplots/ntplots)

treat_dom<-Trelave%>%
  left_join(treat_freq)%>%
  mutate(freq=replace_na(freq, 0)) %>% 
  mutate(treatDCi=(mean+freq)/2)%>%
  select(site_code, project_name, community_type, treatment, species_matched, mean, freq, treatDCi)


# Calculating changes in DCi for each species comparing control to treatments. Dropping extra treatments in CDR, and pulse experiments, and dropping species that are not found in either control or treatments.
CT_diff<-treat_dom%>%
  full_join(control_dom)%>%
  mutate(diff=treatDCi-DCi)%>%
  filter(!is.na(treatment))%>%
  right_join(trt_analysis) %>% 
  mutate(drop=ifelse(site_code=="CDR"&treatment==2|site_code=="CDR"&treatment==3|site_code=="CDR"&treatment==4|site_code=="CDR"&treatment==5|site_code=="CDR"&treatment==7, 1, ifelse(pulse==1, 1, ifelse(treatDCi==0&DCi==0, 1, 0))))%>%
  filter(drop==0)

theme_set(theme_bw(20))
ggplot(data=CT_diff, aes(x=diff))+
  geom_histogram(binwidth = .08)+
  xlab("DCi Difference")+
  ylab("Count")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0)

#dataset of treatment responses, ave, se, and how often species is found for phylogenetic analyses. and calculating number of sites and experiments at treatment is at.

##Temperature
temp_sites<-CT_diff%>%
  filter(temp==1)%>%
  select(site_code)%>%
  unique()

temp_experimenets<-CT_diff%>%
  filter(temp==1)%>%
  select(site_code, project_name, community_type)%>%
  unique()

temp_mean<-CT_diff%>%
  filter(temp==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff)) %>% 
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="temp")%>%
  select(-sd)

##CO2
co2_sites<-CT_diff%>%
  filter(CO2==1)%>%
  select(site_code)%>%
  unique()

co2_experimenets<-CT_diff%>%
  filter(CO2==1)%>%
  select(site_code, project_name, community_type)%>%
  unique()

co2_mean<-CT_diff%>%
  filter(CO2==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff)) %>% 
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="CO2")%>%
  select(-sd)


##Irrigation
irg_sites<-CT_diff%>%
  filter(irg==1)%>%
  select(site_code)%>%
  unique()

irg_experimenets<-CT_diff%>%
  filter(irg==1)%>%
  select(site_code, project_name, community_type)%>%
  unique()

irg_mean<-CT_diff%>%
  filter(irg==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff)) %>% 
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="irg")%>%
  select(-sd)


##Drought
drt_sites<-CT_diff%>%
  filter(drought==1)%>%
  select(site_code)%>%
  unique()

drt_experimenets<-CT_diff%>%
  filter(drought==1)%>%
  select(site_code, project_name, community_type)%>%
  unique()

drt_mean<-CT_diff%>%
  filter(drought==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff)) %>% 
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="drt")%>%
  select(-sd)


##N
n_sites<-CT_diff%>%
  filter(n==1)%>%
  select(site_code)%>%
  unique()

n_experimenets<-CT_diff%>%
  filter(n==1)%>%
  select(site_code, project_name, community_type)%>%
  unique()

n_mean<-CT_diff%>%
  filter(n==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff)) %>% 
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="n")%>%
  select(-sd)

##P
p_sites<-CT_diff%>%
  filter(p==1)%>%
  select(site_code)%>%
  unique()

p_experimenets<-CT_diff%>%
  filter(p==1)%>%
  select(site_code, project_name, community_type)%>%
  unique()

p_mean<-CT_diff%>%
  filter(p==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff)) %>% 
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="p")%>%
  select(-sd)


# ##multiple treatments only use species that are found in 3 or more experiments
## deciding to no longer do this b/c we are changing the mult trt
# 
# allmult_subset<-CT_diff%>%
#   filter(multtrts==1)%>%
#   ungroup()%>%
#   select(species_matched, trt_type)%>%
#   unique()%>%
#   group_by(species_matched)%>%
#   summarize(n=length(trt_type))%>%
#   filter(n>2)%>%
#   select(-n)

allmult_sites<-CT_diff%>%
  filter(multtrts==1)%>%
  select(site_code)%>%
  unique()

allmult_exmt<-CT_diff%>%
  filter(multtrts==1)%>%
  select(site_code, project_name, community_type)%>%
  unique()%>%
  mutate(presmult=1)


allmult_mean<-CT_diff%>%
  filter(multtrts==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="all mult")%>%
  select(-sd)
 
##what are the treatments here, what % have N?
allmult_treats<-CT_diff%>%
  filter(multtrts==1)%>%
  right_join(allmult_subset) %>% 
  select(site_code, project_name, community_type, trt_type) %>% 
  unique() %>% 
  group_by(trt_type)%>%
  summarize(nobs=length(community_type))

sum(allmult_treats$nobs)

###are certain families only in some treatments?
fam<-read.csv("C:\\Users\\mavolio2\\Dropbox\\sDiv_sCoRRE_shared\\CoRRE data\\trait data\\species_families_trees_2021.csv")

allmult_sp<-allmult_mean %>% 
  left_join(fam)


##multiple nutrients
allnuts_sites<-CT_diff%>%
  filter(multnuts==1)%>%
  select(site_code)%>%
  unique()

allnuts_exmt<-CT_diff%>%
  filter(multnuts==1)%>%
  select(site_code, project_name, community_type, trt_type)%>%
  unique()%>%
  mutate(presmult=1)


allnut_mean<-CT_diff%>%
  filter(multnuts==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="all nuts")%>%
  select(-sd)


Fulldataset<-allmult_mean%>%
  bind_rows(allnut_mean, co2_mean, drt_mean, irg_mean, n_mean, p_mean, temp_mean)

write.csv(Fulldataset, paste(my.wd, "WinnersLosers paper/data/Species_DCiDiff_Nov2022.csv", sep=""), row.names=F)

###checking normality of DCI values
ggplot(data=Fulldataset, aes(x=ave_diff))+
  geom_histogram()+
  facet_wrap(~trt_type2)

#figure
aves<-Fulldataset%>%
  group_by(trt_type2)%>%
  summarize(mdiff=mean(ave_diff),
            nobs=length(ave_diff),
            sd=sd(ave_diff))%>%
  mutate(se=sd/sqrt(nobs))

theme_set(theme_bw(14))
ggplot(data=aves, aes(x=trt_type2, y=mdiff, label=nobs))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Treatment")+
  ylab("Average Trt-Ctl Diff")+
  geom_errorbar(aes(ymin=mdiff-se, ymax=mdiff+se), position = position_dodge(0.9), width=0.5)+
  theme(axis.text.x = element_text(angle = 90))+
  geom_text(y=0.01)


# #For this analysis I am using only species that are found in treatments and controls, using differneces btwn trt controls, and the more stringent approach to rare/dominant, 
# categories<-CT_diff%>%
#   na.omit%>%
#   filter(treatDCi!=0|DCi!=0)%>%
#   mutate(cat=ifelse(DCi>0.63&treatDCi<0.08, "Sup. Lose",
#                     ifelse(DCi<0.08&treatDCi<0.08, "Rare Tol.", 
#                            ifelse(DCi>0.63&treatDCi>0.63, "Dom. Tol.", 
#                                   ifelse(DCi>0.63&treatDCi<0.63, "Loose Dom.",
#                                          ifelse(DCi<0.08&treatDCi>0.08, "No Long. Rare",
#                                                 ifelse(diff>0.1, "Winner", 
#                                                        ifelse(absdiff>0.1&diff<0, "Loser", 
#                                                               ifelse(absdiff<0.1, "Tol.", "none")))))))),
#          DR=ifelse(DCi>0.63, "Dom", 
#                    ifelse(DCi<0.08, "Rare", "none")),
#          Pres=ifelse(DCi>0.08&DCi<0.63, 1, 0))

#old code
# ifelse(DCi!=0&DCi<0.08&treatDCi>0.6, "Sup. Win.", 
#        ifelse(DCi==0&treatDCi>0.6, "Sup. Dup. Win", 
# ifelse(DCi==0&treatDCi<0.08, "Appear Rare",
#        ifelse(DCi==0&treatDCi>0.08, "Appear", 
# AppDis=ifelse(DCi==0&treatDCi>0, "Appear",
#               ifelse(DCi>0&treatDCi==0, "Disapp.", "none")),

# #getting a table of responses
# cat_sum<-categories%>%
#   group_by(cat)%>%
#   summarize(n=length(cat))%>%
#   mutate(prop=n/26497)
# 
# # dr_sum<-categories%>%
#   group_by(DR)%>%
#   summarize(n=length(DR))
# 
# ad_sum<-categories%>%
#   group_by(AppDis)%>%
#   summarize(n=length(AppDis))
# 
# sum(categories$Pres)


# ###contrasting different treatments
trts2<-trts%>%
  select(site_code, project_name, community_type, treatment, trt_type)%>%
  unique
# 
# cat_trts<-categories%>% 
#   left_join(trts2)%>%
#   mutate(cat2=ifelse(cat=="Appear"|cat=="Winner"|cat=="No Long. Rare"|cat=="Appear Rare"|cat=="Sup. Dup. Win", "Winner", ifelse(cat=="Loser"|cat=="Loose Dom.", "Loser", ifelse(cat=="Tol."|cat=="Dom. Tol."|cat=="Rare Tol.", "Tol.", "None"))))%>%
#   group_by(trt_type)%>%
#   summarize(avediff=mean(diff))
# 
# ntrt<-cat_trts%>%
#   filter(trt_type=="N")
# ntrt_sum<-ntrt%>%
#   group_by(trt_type, cat2)%>%
#   summarize(n=length(cat2))%>%
#   mutate(prop=n/2249)
# 
# drttrt<-cat_trts%>%
#   filter(trt_type=="drought")
# drttrt_sum<-drttrt%>%
#   group_by(trt_type,cat2)%>%
#   summarize(n=length(cat2))%>%
#   mutate(prop=n/939)
# 
# irrtrt<-cat_trts%>%
#   filter(trt_type=="irr")
# irrtrt_sum<-irrtrt%>%
#   group_by(trt_type,cat2)%>%
#   summarize(n=length(cat2))%>%
#   mutate(prop=n/798)
# 
# multnuttrt<-cat_trts%>%
#   filter(trt_type=="mult_nutrient")
# multnuttrt_sum<-multnuttrt%>%
#   group_by(trt_type,cat2)%>%
#   summarize(n=length(cat2))%>%
#   mutate(prop=n/3548)
# 
# trt_contrast<-ntrt_sum%>%
#   bind_rows(drttrt_sum, irrtrt_sum, multnuttrt_sum)
# 
# theme_set(theme_bw(16))
# 
# ggplot(data=trt_contrast, aes(x=trt_type, y=prop, fill=cat2))+
#   geom_bar(stat="identity", position = position_dodge())+
#   xlab("Treatment")+
#   ylab("Proportion of Sp. Responses")+
#   scale_fill_manual(name="Response", values=c("red", "gray",'blue'))
# 
# ggplot(data=filter(cat_trts, trt_type=="N"|trt_type=="mult_nutrient"|trt_type=="irr"|trt_type=="drought"), aes(x=trt_type, y=avediff))+
#   geom_bar(stat="identity", position = position_dodge())+
#   xlab("Treatment")+
#   ylab("Average Control-Trt Diff")
# 
# #exploring species specfic responses
# ###example with A. gerardii
# ange<-categories%>%
#   filter(species_matched=="Andropogon gerardii")%>%
#   left_join(trts2)
# 
# 
# 
# angecat_sum<-ange%>%
#   group_by(cat)%>%
#   summarize(n=length(cat))%>%
#   mutate(prop=n/130)
# 
# ange_trt<-ange%>%
#   filter(trt_type=="N"|trt_type=="mult_nutrient"|trt_type=="irr"|trt_type=="drought")%>%
#   group_by(trt_type)%>%
#   summarize(ave_diff=mean(diff),
#             n=length(diff),
#             sd=sd(diff))%>%
#   mutate(se=sd/sqrt(n))
# 
# ggplot(data=ange_trt, aes(x=trt_type, y=ave_diff, label=n))+
#   geom_bar(stat="identity", position = position_dodge())+
#   xlab("Treatment")+
#   ylab("Average Control-Trt Diff")+
#   geom_label(hjust=1)+
#   geom_errorbar(aes(ymin=ave_diff-se, ymax=ave_diff+se), position = position_dodge(0.9), width=0.5)
# 
# 
# 
#common species - get list of common species to explore
common<-categories%>%
  group_by(species_matched)%>%
  summarize(n=length(species_matched))

##erigeron canadensis
erca<-Fulldataset%>%
  filter(species_matched=="Erigeron canadensis")%>%
  left_join(trts2)

ercacat_sum<-erca%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/209)

erca_trt<-erca%>%
  filter(trt_type=="N"|trt_type=="mult_nutrient"|trt_type=="irr"|trt_type=="drought")%>%
  group_by(trt_type)%>%
  summarize(ave_diff=mean(diff),
            n=length(diff),
            sd=sd(diff))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=erca_trt, aes(x=trt_type, y=ave_diff, label=n))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Treatment")+
  ylab("Average Control-Trt Diff")+
  geom_label(hjust=1)+
  geom_errorbar(aes(ymin=ave_diff-se, ymax=ave_diff+se), position = position_dodge(0.9), width=0.5)

#possibilties
#plantago lanceolata
#anthoxanthus odoratum
#Achillea millefolium 155 examples
#briza minor
#Bromus hordeaceus - 58 obs varied responses
#Chondrosum gracile - 90 obs nice split and varied responses

#looking at a few examples
test<-Fulldataset_mixedmodels%>%
  filter(species_matched=="Achillea millefolium"|species_matched=="Bromus hordeaceus"|species_matched=='Chondrosum gracile')

test1<-test%>%
  filter(trt_type2=="n"|trt_type2=="all mult"|trt_type2=="co2"|trt_type2=="drought"|trt_type2=="temp")%>%
  group_by(species_matched,trt_type2)%>%
  summarize(ave_diff=mean(diff),
            n=length(diff),
            sd=sd(diff))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=test1, aes(x=trt_type2, y=ave_diff, fill=species_matched, group=species_matched))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Treatment")+
  ylab("Trt-Cntrl Difference")+
  #geom_text(aes(label=n), position=position_dodge(.6), vjust=1)+
  geom_errorbar(aes(ymin=ave_diff-se, ymax=ave_diff+se), position = position_dodge(0.9), width=0.2)+
  scale_x_discrete(limits=c('n', 'co2', 'drought', 'temp', 'all mult'), labels=c('N', 'CO2', 'Drought', 'Temp.', 'Mult Trts'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(name="Species", values=c('darkblue', 'darkgreen', 'green1'))


###species repsonses
common_site<-categories%>%
  ungroup()%>%
  select(species_matched, site_code)%>%
  unique()%>%
  group_by(species_matched)%>%
  summarize(n=length(species_matched))%>%
  filter(n>6)%>%
  select(-n)

#overall response
common_trt_sp_overall<-categories%>%
  left_join(trts2)%>%
  right_join(common_site)%>%
  filter(site_code!="Sil")%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),n=length(diff),
            sd=sd(diff))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(trt_type="overall")

##sp subset analysis
common_trt_sp<-categories%>%
  left_join(trts2)%>%
  right_join(common_site)%>%
  filter(trt_type=="mult_nutrient"|trt_type=="N"|trt_type=="irr"|trt_type=="drought")%>%
  group_by(species_matched, trt_type)%>%
  summarize(ave_diff=mean(diff),
            n=length(diff),
            sd=sd(diff))%>%
  mutate(se=sd/sqrt(n))%>%
  bind_rows(common_trt_sp_overall)

##big figure to share
ggplot(data=common_trt_sp, aes(x=trt_type, y=ave_diff, label=n))+
  geom_bar(stat="identity", position = position_dodge())+
  geom_text(y=-0.3)+
  xlab("Treatment")+
  ylab("Average Control-Trt Diff")+
  geom_errorbar(aes(ymin=ave_diff-se, ymax=ave_diff+se), position = position_dodge(0.9), width=0.5)+
  theme(axis.text.x = element_text(angle = 90))+
  geom_vline(aes(xintercept = 4.5))+
  facet_wrap(~species_matched)

# 
# andro<-categories%>%
#   left_join(trts2)%>%
#   filter(species_matched=="Andropogon gerardii")%>%
#   filter(trt_type=="N*P"|trt_type=="N"|trt_type=='drought'|trt_type=="irr"|trt_type=="mult_nutrient")
# 
# ggplot(data=andro, aes(x=trt_type, y=diff, fill=project_name))+
#   geom_bar(stat="identity", position = position_dodge())+
#   xlab("Treatment")+
#   ylab("Average Control-Trt Diff")+
#   theme(axis.text.x = element_text(angle = 90))+
#   facet_wrap(~site_code, scale="free_x")
# 
# ###now doing this for every year of an experiment.
# ####I did this for our Sept 2020 meeting, there were no major differences and I think we should not explore this further
# 
# relave_yr<-allreldat%>%
#   group_by(site_code, project_name, community_type, treatment, plot_mani, treatment_year, species_matched)%>%
#   summarize(mean=mean(relcov))
# 
# Crelave_yr<-relave_yr%>%
#   filter(plot_mani==0)%>%
#   ungroup()%>%
#   select(-treatment, -plot_mani)%>%
#   rename(cmean=mean)
# 
# Trelave_yr<-relave_yr%>%
#   filter(plot_mani!=0)%>%
#   ungroup()%>%
#   select(-plot_mani)
# 
# controlplots_yr<-allreldat%>%
#   filter(plot_mani==0)%>%
#   select(site_code, project_name, community_type, plot_id, treatment_year)%>%
#   unique()%>%
#   group_by(site_code, project_name, community_type, treatment_year)%>%
#   summarize(ncplots=length(plot_id))
# 
# control_freq_yr<-allreldat%>%
#   filter(plot_mani==0)%>%
#   select(site_code, project_name, community_type, species_matched, plot_id, treatment_year)%>%
#   unique()%>%
#   group_by(site_code, project_name, community_type, species_matched, treatment_year)%>%
#   summarize(nplots=length(plot_id))%>%
#   left_join(controlplots)%>%
#   mutate(freq=nplots/ncplots)
# 
# control_dom_yr<-control_freq_yr%>%
#   left_join(Crelave_yr)%>%
#   mutate(DCi=(cmean+freq)/2)%>%
#   select(site_code, project_name, community_type, species_matched, treatment_year, DCi)
# 
# treatplots_yr<-allreldat%>%
#   filter(plot_mani!=0)%>%
#   select(site_code, project_name, community_type, treatment, plot_id, treatment_year)%>%
#   unique()%>%
#   group_by(site_code, project_name, community_type, treatment, treatment_year)%>%
#   summarize(ntplots=length(plot_id))
# 
# treat_freq_yr<-allreldat%>%
#   filter(plot_mani!=0)%>%
#   select(site_code, project_name, community_type, treatment, species_matched, plot_id, treatment_year)%>%
#   unique()%>%
#   group_by(site_code, project_name, community_type, treatment, species_matched, treatment_year)%>%
#   summarize(nplots=length(plot_id))%>%
#   left_join(treatplots)%>%
#   mutate(freq=nplots/ntplots)
# 
# treat_dom_yr<-treat_freq_yr%>%
#   left_join(Trelave_yr)%>%
#   mutate(treatDCi=(mean+freq)/2)%>%
#   select(site_code, project_name, community_type, treatment, species_matched, treatment_year, treatDCi)
# 
# CT_yr<-treat_dom_yr%>%
#   full_join(control_dom_yr)#this introduces NA into the treatment column when it is present in the controls but not the treated plots.
# 
# CT_yr$DCi[is.na(CT_yr$DCi)] <- 0
# CT_yr$treatDCi[is.na(CT_yr$treatDCi)] <- 0
# 
# 
# CT_diff_yr<-CT_yr%>%
#   mutate(diff=treatDCi-DCi,
#          absdiff=abs(diff))
# 
# categories_yr<-CT_diff_yr%>%
#   na.omit%>%
#   mutate(cat=ifelse(DCi!=0&DCi<0.08&treatDCi>0.6, "Sup. Win.",
#              ifelse(DCi==0&treatDCi>0.6, "Sup. Dup. Win",
#              ifelse(DCi>0.6&treatDCi!=0&treatDCi<0.1, "Sup. Lose",
#              ifelse(DCi!=0&DCi<0.08&treatDCi!=0&treatDCi<0.08, "Rare Tol.",
#              ifelse(DCi>0.6&treatDCi>0.6, "Dom. Tol.",
#              ifelse(DCi==0&treatDCi<0.08, "Appear Rare",
#              ifelse(DCi==0&treatDCi>0.08, "Appear",
#              ifelse(DCi>0.6&treatDCi<0.6, "Loose Dom.",
#              ifelse(DCi<0.08&treatDCi>0.08, "No Long. Rare",
#              ifelse(diff>0.1, "Winner",
#              ifelse(absdiff>0.1&diff<0, "Loser",
#              ifelse(absdiff<0.1, "Tol.", "none")))))))))))),
#          DR=ifelse(DCi>0.6, "Dom",
#                    ifelse(DCi!=0&DCi<0.08, "Rare", "none")),
#          AppDis=ifelse(DCi==0&treatDCi>0, "Appear",
#                        ifelse(DCi>0&treatDCi==0, "Disapp.", "none")),
#          Pres=ifelse(DCi!=0&DCi>0.08&DCi<0.6, 1, 0),
#          cat2=ifelse(cat=="Appear"|cat=="Winner"|cat=="No Long. Rare"|cat=="Appear Rare"|cat=="Sup. Dup. Win", "Winner", ifelse(cat=="Loser"|cat=="Loose Dom."|cat=="Sup. Lose", "Loser", ifelse(cat=="Tol."|cat=="Dom. Tol."|cat=="Rare Tol.", "Tol.", "None"))))
# 
# cat_sum_yr<-categories_yr%>%
#   group_by(cat)%>%
#   summarize(n=length(cat))%>%
#   mutate(prop=n/72221)
# 
# dr_sum_yr<-categories_yr%>%
#   group_by(DR)%>%
#   summarize(n=length(DR))
# 
# ad_sum_yr<-categories_yr%>%
#   group_by(AppDis)%>%
#   summarize(n=length(AppDis))
# 
# cat2_sum_yr<-categories_yr%>%
#   group_by(cat2)%>%
#   summarise(n=length(cat2))%>%
#   mutate(prop=n/72221)
# 
# ####looking into trends through time for andro
# ange<-CT_diff_yr%>%
#   filter(species_matched=="Andropogon gerardii")%>%
#   left_join(trts2)
# 
# ggplot(data=ange, aes(x=treatment_year, y=diff, color=treatment))+
#   geom_point()+
#   geom_smooth(method="lm")+
#   facet_wrap(~trt_type, scales="free")+
#   theme(legend.position = "none")
# 
# koma<-CT_diff_yr%>%
#   filter(species_matched=="Koeleria macrantha")%>%
#   left_join(trts2)
# 
# ggplot(data=koma, aes(x=treatment_year, y=diff, color=treatment))+
#   geom_point()+
#   geom_smooth(method="lm")+
#   facet_wrap(~trt_type, scales="free")+
#   theme(legend.position = "none")


####
###dataset of all responses for trait mixed models
CT_Sp_temp<-CT_diff%>%
  filter(temp==1)%>%
  mutate(trt_type2="temp")

CT_Sp_co2<-CT_diff%>%
  filter(CO2==1)%>%
  mutate(trt_type2="co2")

CT_Sp_irg<-CT_diff%>%
  filter(irg==1)%>%
  mutate(trt_type2="irrigation")

CT_Sp_drt<-CT_diff%>%
  filter(drought==1)%>%
  mutate(trt_type2="drought")

CT_Sp_N<-CT_diff%>%
  filter(n==1)%>%
  mutate(trt_type2="n")

CT_Sp_P<-CT_diff%>%
  filter(p==1)%>%
  mutate(trt_type2="p")

CT_Sp_multnuts<-CT_diff%>%
  filter(multnuts==1)%>%
  mutate(trt_type2="multnuts")

CT_Sp_allint<-CT_diff%>%
  filter(multtrts==1)%>%
  mutate(trt_type2="all mult")

Fulldataset_mixedmodels<-CT_Sp_allint%>%
  bind_rows(CT_Sp_co2, CT_Sp_drt, CT_Sp_irg, CT_Sp_N, CT_Sp_P, CT_Sp_temp, CT_Sp_multnuts)%>%
  select(site_code, project_name, community_type, species_matched, trt_type2, diff)

write.csv(Fulldataset_mixedmodels, paste(my.wd, "WinnersLosers paper/data/Species_DCiDiff_formixedmodelsNov22.csv", sep=""), row.names=F)


