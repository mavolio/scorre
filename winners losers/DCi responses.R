
#investigate winners losers

library(tidyverse)
library(gridExtra)

###read in data

my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "E:/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"

#read in the data

#raw abundance data
dat<-read.csv(paste(my.wd, "CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Dec2021.csv",sep=""))

sp <-read.csv(paste(my.wd,"CoRRE data/CoRRE data/trait data/CoRRE2trykey_2021.csv", sep=""))%>%
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
  mutate(alltrts=ifelse(trt_type %in% c("CO2","CO2*temp", "mow_clip","burn","burn*graze","disturbance","burn*mow_clip","drought","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N","mult_nutrient","N*P","P","N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp","temp*mow_clip","drought*temp","irr*temp","irr"),1,0))%>%
  filter(alltrts==1)%>%
  mutate(dist=ifelse(trt_type %in% c("mow_clip","burn","burn*graze","disturbance","burn*mow_clip"), 1, 0), 
         dist_other=ifelse(trt_type %in% c("drought*mow_clip
","drought*temp*mow_clip", "irr*mow_clip","irr*temp*mow_clip","N*irr*mow_clip","N*P*burn*graze","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip", "mult_nutrient*mow_clip","N*burn*mow_clip", "N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp*mow_clip"), 1, 0),
         CO2=ifelse(trt_type %in% c("CO2"), 1, 0),
         CO2_other=ifelse(trt_type %in% c("CO2*temp","drought*CO2*temp","irr*CO2","irr*CO2*temp","N*CO2*temp","N*irr*CO2","N*irr*CO2*temp","N*CO2"), 1, 0),
         drought=ifelse(trt_type %in% c("drought"), 1, 0),
         drought_other=ifelse(trt_type %in% c("drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","N*drought","drought*temp"), 1, 0),
         herb_removal=ifelse(trt_type %in% c("herb_removal"), 1, 0),
         herb_removal_other=ifelse(trt_type %in% c("herb_removal*mow_clip","irr*herb_removal","N*herb_removal","P*herb_removal","mult_nutrient*herb_removal*mow_clip","mult_nutrient*herb_removal"), 1, 0),
         irg=ifelse(trt_type %in% c("irr"), 1, 0),
         irg_other=ifelse(trt_type %in% c("irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*irr*CO2","N*irr*mow_clip","mult_nutrient*irr","N*irr*CO2*temp","N*irr","N*irr*temp","irr*temp"), 1, 0),
         temp=ifelse(trt_type %in% c("temp"), 1, 0),
         temp_other=ifelse(trt_type %in% c("CO2*temp","drought*CO2*temp","drought*temp*mow_clip","irr*CO2*temp","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2*temp","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","temp*mow_clip","drought*temp","irr*temp"), 1, 0),
        # nuts_other=ifelse(trt_type %in% c("N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze","mult_nutrient*irr","N*irr*CO2*tempN*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip"), 1, 0),
         n=ifelse(trt_type=="N", 1, 0),
         n_other=ifelse(trt_type %in% c("N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze","N*irr*CO2*tempN*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","N*drought","N*herb_removal","N*irr","N*irr*temp","N*temp","N*P*temp","N*burn*mow_clip","N*P*burn","N*P*mow_clip"), 1, 0),
         p=ifelse(trt_type=="P", 1, 0),
         p_other=ifelse(trt_type %in% c("N*P*burn*graze","mult_nutrient*irr","P*burn*graze","P*burn*mow_clip","P*herb_removal","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal*mow_clip"), 1, 0),
         multtrts=ifelse(trt_type %in% c("mult_nutrient","N*P","CO2*temp", "burn*graze","burn*mow_clip","drought*CO2*temp","drought*mow_clip","drought*temp*mow_clip","herb_removal*mow_clip","irr*CO2","irr*CO2*temp","irr*mow_clip","irr*herb_removal","irr*temp*mow_clip","N*CO2*temp","N*irr*CO2","N*irr*mow_clip","N*P*burn*graze", "mult_nutrient*irr","N*irr*CO2*temp", "N*CO2","N*mow_clip","N*burn","N*burn*graze","N*disturbance","P*burn*graze","P*burn*mow_clip","N*drought","N*herb_removal","P*herb_removal","N*irr","N*irr*temp","N*temp","mult_nutrient*temp","N*P*temp","mult_nutrient*mow_clip","N*burn*mow_clip","N*P*burn","N*P*mow_clip","P*burn","P*mow_clip","mult_nutrient*herb_removal","mult_nutrient*herb_removal*mow_clip","temp*mow_clip","drought*temp","irr*temp"),1,0))


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
Crelave<-alldat%>%
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

#to get relative frequency, determine number of control plots a species is found in, merge in total number of plots and calculate relative frequency  
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

treat_dom<-treat_freq%>%
  left_join(Trelave)%>%
  mutate(treatDCi=(mean+freq)/2)%>%
  select(site_code, project_name, community_type, treatment, species_matched, mean, freq, treatDCi)


# Calculating changes in DCi for each species comparing control to treatments.
CT<-treat_dom%>%
  full_join(control_dom)#this introduces NA into the treatment column when it is present in the controls but not the treated plots.

#this is to replace the NAs with 0
CT$DCi[is.na(CT$DCi)] <- 0
CT$treatDCi[is.na(CT$treatDCi)] <- 0

#add treatments to plots were was only in control
CT_diff_control<-CT%>%
  mutate(diff=treatDCi-DCi)%>%
  filter(is.na(treatment))%>%#dropping control only plots
  ungroup()%>%
  select(-treatment)%>%
  left_join(trt_analysis)%>%
  filter(!is.na(treatment))

CT_diff<-CT%>%
  mutate(diff=treatDCi-DCi)%>%
  filter(!is.na(treatment))%>%
  right_join(trt_analysis)%>%
  filter(!is.na(diff))%>%# somthing about adding controls adds NA for treatemtns were there were control plots kept. investigate this. 3 datasets are missing, GVN_Face, KNZ_GFP, RIO_interaction and SCL_Lucero
  bind_rows(CT_diff_control)%>%
  mutate(drop=ifelse(site_code=="Sil"&resource_mani==0, 1, ifelse(site_code=="CDR"&treatment==2|site_code=="CDR"&treatment==3|site_code=="CDR"&treatment==4|site_code=="CDR"&treatment==5|site_code=="CDR"&treatment==7, 1, ifelse(pulse==1, 1, 0))))%>%
  filter(drop==0)

#dataset of treatment responses, ave, se, min, max, and how often species is found for phylogenetic analyses.

##
CT_Sp_herb<-CT_diff%>%
  filter(herb_removal==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="herb_removal")

CT_Sp_herb_other<-CT_diff%>%
  filter(herb_removal_other==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="herb_rem_other")

CT_Sp_temp<-CT_diff%>%
  filter(temp==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="temp")

CT_Sp_temp_other<-CT_diff%>%
  filter(temp_other==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="temp_other")

CT_Sp_co2<-CT_diff%>%
  filter(CO2==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="co2")
 

CT_Sp_co2_other<-CT_diff%>%
  filter(CO2_other==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="co2_other")

CT_Sp_dist<-CT_diff%>%
  filter(dist==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="disturbance")

CT_Sp_dist_other<-CT_diff%>%
  filter(dist_other==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="dist_other")

CT_Sp_irg<-CT_diff%>%
  filter(irg==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="irrigation")

CT_Sp_irg_other<-CT_diff%>%
  filter(irg_other==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="irg_other")

CT_Sp_drt<-CT_diff%>%
  filter(drought==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="drought")

CT_Sp_drt_other<-CT_diff%>%
  filter(drought_other==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="drt_other")

CT_Sp_N<-CT_diff%>%
  filter(n==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="n")

CT_Sp_P<-CT_diff%>%
  filter(p==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="p")

CT_Sp_P_other<-CT_diff%>%
  filter(p_other==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="p_other")

# CT_Sp_nuts_other<-CT_diff%>%
#   filter(nuts_other==1)%>%
#   group_by(species_matched)%>%
#   summarize(ave_diff=mean(diff),
#             nobs=length(diff),
#             sd=sd(diff),
#             min=min(diff),
#             max=max(diff))%>%
#   mutate(se=sd/sqrt(nobs))%>%
#   mutate(trt_type2="nuts_other")

CT_Sp_n_other<-CT_diff%>%
  filter(n_other==1)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="n_other")

allmul_subset<-CT_diff%>%
  filter(multtrts==1)%>%
  ungroup()%>%
  select(species_matched, trt_type)%>%
  unique()%>%
  group_by(species_matched)%>%
  summarize(n=length(trt_type))%>%
  filter(n>2)%>%
  select(-n)

CT_Sp_allint<-CT_diff%>%
  filter(multtrts==1)%>%
  right_join(allmul_subset)%>%
  group_by(species_matched)%>%
  summarize(ave_diff=mean(diff),
            nobs=length(diff),
            sd=sd(diff),
            min=min(diff),
            max=max(diff))%>%
  mutate(se=sd/sqrt(nobs))%>%
  mutate(trt_type2="all mult")


Fulldataset<-CT_Sp_allint%>%
  bind_rows(CT_Sp_co2, CT_Sp_co2_other, CT_Sp_dist, CT_Sp_dist_other, CT_Sp_drt,CT_Sp_drt_other,  CT_Sp_herb, CT_Sp_herb_other, CT_Sp_irg, CT_Sp_irg_other, CT_Sp_N, CT_Sp_n_other, CT_Sp_P, CT_Sp_P_other, CT_Sp_temp, CT_Sp_temp_other)%>%
  select(species_matched, trt_type2, nobs, ave_diff, min, max, se)

write.csv(Fulldataset, paste(my.wd, "WinnersLosers paper/data/Species_DCiDiff_Dec2021.csv", sep=""), row.names=F)

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


###contrasting different treatments
trts2<-trts%>%
  select(site_code, project_name, community_type, treatment, trt_type)%>%
  unique

cat_trts<-categories%>% 
  left_join(trts2)%>%
  mutate(cat2=ifelse(cat=="Appear"|cat=="Winner"|cat=="No Long. Rare"|cat=="Appear Rare"|cat=="Sup. Dup. Win", "Winner", ifelse(cat=="Loser"|cat=="Loose Dom.", "Loser", ifelse(cat=="Tol."|cat=="Dom. Tol."|cat=="Rare Tol.", "Tol.", "None"))))%>%
  group_by(trt_type)%>%
  summarize(avediff=mean(diff))

ntrt<-cat_trts%>%
  filter(trt_type=="N")
ntrt_sum<-ntrt%>%
  group_by(trt_type, cat2)%>%
  summarize(n=length(cat2))%>%
  mutate(prop=n/2249)

drttrt<-cat_trts%>%
  filter(trt_type=="drought")
drttrt_sum<-drttrt%>%
  group_by(trt_type,cat2)%>%
  summarize(n=length(cat2))%>%
  mutate(prop=n/939)

irrtrt<-cat_trts%>%
  filter(trt_type=="irr")
irrtrt_sum<-irrtrt%>%
  group_by(trt_type,cat2)%>%
  summarize(n=length(cat2))%>%
  mutate(prop=n/798)

multnuttrt<-cat_trts%>%
  filter(trt_type=="mult_nutrient")
multnuttrt_sum<-multnuttrt%>%
  group_by(trt_type,cat2)%>%
  summarize(n=length(cat2))%>%
  mutate(prop=n/3548)

trt_contrast<-ntrt_sum%>%
  bind_rows(drttrt_sum, irrtrt_sum, multnuttrt_sum)

theme_set(theme_bw(16))

ggplot(data=trt_contrast, aes(x=trt_type, y=prop, fill=cat2))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Treatment")+
  ylab("Proportion of Sp. Responses")+
  scale_fill_manual(name="Response", values=c("red", "gray",'blue'))

ggplot(data=filter(cat_trts, trt_type=="N"|trt_type=="mult_nutrient"|trt_type=="irr"|trt_type=="drought"), aes(x=trt_type, y=avediff))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Treatment")+
  ylab("Average Control-Trt Diff")

#exploring species specfic responses
###example with A. gerardii
ange<-categories%>%
  filter(species_matched=="Andropogon gerardii")%>%
  left_join(trts2)



angecat_sum<-ange%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/130)

ange_trt<-ange%>%
  filter(trt_type=="N"|trt_type=="mult_nutrient"|trt_type=="irr"|trt_type=="drought")%>%
  group_by(trt_type)%>%
  summarize(ave_diff=mean(diff),
            n=length(diff),
            sd=sd(diff))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=ange_trt, aes(x=trt_type, y=ave_diff, label=n))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Treatment")+
  ylab("Average Control-Trt Diff")+
  geom_label(hjust=1)+
  geom_errorbar(aes(ymin=ave_diff-se, ymax=ave_diff+se), position = position_dodge(0.9), width=0.5)



#common species - get list of common species to explore
common<-categories%>%
  group_by(species_matched)%>%
  summarize(n=length(species_matched))

##erigeron canadensis
erca<-categories%>%
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


#for poa pratensis
popr<-categories%>%
  filter(species_matched=="Poa pratensis")%>%
  left_join(trts2)

poprcat_sum<-popr%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/233)

popr_trt<-popr%>%
  filter(trt_type=="N"|trt_type=="mult_nutrient"|trt_type=="irr"|trt_type=="drought")%>%
  group_by(trt_type)%>%
  summarize(ave_diff=mean(diff),
            n=length(diff),
            sd=sd(diff))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=popr_trt, aes(x=trt_type, y=ave_diff, label=n))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Treatment")+
  ylab("Average Control-Trt Diff")+
  geom_label(hjust=1)+
  geom_errorbar(aes(ymin=ave_diff-se, ymax=ave_diff+se), position = position_dodge(0.9), width=0.5)


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


#investigating sil nash
sil<-categories%>%
  filter(site_code=="Sil")%>%
  left_join(trts2)%>%
  select(treatment, trt_type)%>%
  unique()

andro<-categories%>%
  left_join(trts2)%>%
  filter(species_matched=="Andropogon gerardii")%>%
  filter(trt_type=="N*P"|trt_type=="N"|trt_type=='drought'|trt_type=="irr"|trt_type=="mult_nutrient")

ggplot(data=andro, aes(x=trt_type, y=diff, fill=project_name))+
  geom_bar(stat="identity", position = position_dodge())+
  xlab("Treatment")+
  ylab("Average Control-Trt Diff")+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~site_code, scale="free_x")

###now doing this for every year of an experiment.
####I did this for our Sept 2020 meeting, there were no major differences and I think we shoudl not explore this further

relave_yr<-allreldat%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani, treatment_year, species_matched)%>%
  summarize(mean=mean(relcov))

Crelave_yr<-relave_yr%>%
  filter(plot_mani==0)%>%
  ungroup()%>%
  select(-treatment, -plot_mani)%>%
  rename(cmean=mean)

Trelave_yr<-relave_yr%>%
  filter(plot_mani!=0)%>%
  ungroup()%>%
  select(-plot_mani)

controlplots_yr<-allreldat%>%
  filter(plot_mani==0)%>%
  select(site_code, project_name, community_type, plot_id, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment_year)%>%
  summarize(ncplots=length(plot_id))

control_freq_yr<-allreldat%>%
  filter(plot_mani==0)%>%
  select(site_code, project_name, community_type, species_matched, plot_id, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, species_matched, treatment_year)%>%
  summarize(nplots=length(plot_id))%>%
  left_join(controlplots)%>%
  mutate(freq=nplots/ncplots)

control_dom_yr<-control_freq_yr%>%
  left_join(Crelave_yr)%>%
  mutate(DCi=(cmean+freq)/2)%>%
  select(site_code, project_name, community_type, species_matched, treatment_year, DCi)

treatplots_yr<-allreldat%>%
  filter(plot_mani!=0)%>%
  select(site_code, project_name, community_type, treatment, plot_id, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment, treatment_year)%>%
  summarize(ntplots=length(plot_id))

treat_freq_yr<-allreldat%>%
  filter(plot_mani!=0)%>%
  select(site_code, project_name, community_type, treatment, species_matched, plot_id, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment, species_matched, treatment_year)%>%
  summarize(nplots=length(plot_id))%>%
  left_join(treatplots)%>%
  mutate(freq=nplots/ntplots)

treat_dom_yr<-treat_freq_yr%>%
  left_join(Trelave_yr)%>%
  mutate(treatDCi=(mean+freq)/2)%>%
  select(site_code, project_name, community_type, treatment, species_matched, treatment_year, treatDCi)

CT_yr<-treat_dom_yr%>%
  full_join(control_dom_yr)#this introduces NA into the treatment column when it is present in the controls but not the treated plots.

CT_yr$DCi[is.na(CT_yr$DCi)] <- 0
CT_yr$treatDCi[is.na(CT_yr$treatDCi)] <- 0


CT_diff_yr<-CT_yr%>%
  mutate(diff=treatDCi-DCi,
         absdiff=abs(diff))

categories_yr<-CT_diff_yr%>%
  na.omit%>%
  mutate(cat=ifelse(DCi!=0&DCi<0.08&treatDCi>0.6, "Sup. Win.",
             ifelse(DCi==0&treatDCi>0.6, "Sup. Dup. Win",
             ifelse(DCi>0.6&treatDCi!=0&treatDCi<0.1, "Sup. Lose",
             ifelse(DCi!=0&DCi<0.08&treatDCi!=0&treatDCi<0.08, "Rare Tol.",
             ifelse(DCi>0.6&treatDCi>0.6, "Dom. Tol.",
             ifelse(DCi==0&treatDCi<0.08, "Appear Rare",
             ifelse(DCi==0&treatDCi>0.08, "Appear",
             ifelse(DCi>0.6&treatDCi<0.6, "Loose Dom.",
             ifelse(DCi<0.08&treatDCi>0.08, "No Long. Rare",
             ifelse(diff>0.1, "Winner",
             ifelse(absdiff>0.1&diff<0, "Loser",
             ifelse(absdiff<0.1, "Tol.", "none")))))))))))),
         DR=ifelse(DCi>0.6, "Dom",
                   ifelse(DCi!=0&DCi<0.08, "Rare", "none")),
         AppDis=ifelse(DCi==0&treatDCi>0, "Appear",
                       ifelse(DCi>0&treatDCi==0, "Disapp.", "none")),
         Pres=ifelse(DCi!=0&DCi>0.08&DCi<0.6, 1, 0),
         cat2=ifelse(cat=="Appear"|cat=="Winner"|cat=="No Long. Rare"|cat=="Appear Rare"|cat=="Sup. Dup. Win", "Winner", ifelse(cat=="Loser"|cat=="Loose Dom."|cat=="Sup. Lose", "Loser", ifelse(cat=="Tol."|cat=="Dom. Tol."|cat=="Rare Tol.", "Tol.", "None"))))

cat_sum_yr<-categories_yr%>%
  group_by(cat)%>%
  summarize(n=length(cat))%>%
  mutate(prop=n/72221)

dr_sum_yr<-categories_yr%>%
  group_by(DR)%>%
  summarize(n=length(DR))

ad_sum_yr<-categories_yr%>%
  group_by(AppDis)%>%
  summarize(n=length(AppDis))

cat2_sum_yr<-categories_yr%>%
  group_by(cat2)%>%
  summarise(n=length(cat2))%>%
  mutate(prop=n/72221)

####looking into trends through time for andro
ange<-CT_diff_yr%>%
  filter(species_matched=="Andropogon gerardii")%>%
  left_join(trts2)

ggplot(data=ange, aes(x=treatment_year, y=diff, color=treatment))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~trt_type, scales="free")+
  theme(legend.position = "none")

koma<-CT_diff_yr%>%
  filter(species_matched=="Koeleria macrantha")%>%
  left_join(trts2)

ggplot(data=koma, aes(x=treatment_year, y=diff, color=treatment))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~trt_type, scales="free")+
  theme(legend.position = "none")


####
###dataset of all responses for trait mixed models
CT_Sp_herb<-CT_diff%>%
  filter(herb_removal==1)%>%
  mutate(trt_type2="herb_removal")

CT_Sp_herb_other<-CT_diff%>%
  filter(herb_removal_other==1)%>%
  mutate(trt_type2="herb_rem_other")

CT_Sp_temp<-CT_diff%>%
  filter(temp==1)%>%
  mutate(trt_type2="temp")

CT_Sp_temp_other<-CT_diff%>%
  filter(temp_other==1)%>%
  mutate(trt_type2="temp_other")

CT_Sp_co2<-CT_diff%>%
  filter(CO2==1)%>%
  mutate(trt_type2="co2")


CT_Sp_co2_other<-CT_diff%>%
  filter(CO2_other==1)%>%
  mutate(trt_type2="co2_other")

CT_Sp_dist<-CT_diff%>%
  filter(dist==1)%>%
  mutate(trt_type2="disturbance")

CT_Sp_dist_other<-CT_diff%>%
  filter(dist_other==1)%>%
  mutate(trt_type2="dist_other")

CT_Sp_irg<-CT_diff%>%
  filter(irg==1)%>%
  mutate(trt_type2="irrigation")

CT_Sp_irg_other<-CT_diff%>%
  filter(irg_other==1)%>%
  mutate(trt_type2="irg_other")

CT_Sp_drt<-CT_diff%>%
  filter(drought==1)%>%
  mutate(trt_type2="drought")

CT_Sp_drt_other<-CT_diff%>%
  filter(drought_other==1)%>%
  mutate(trt_type2="drt_other")

CT_Sp_N<-CT_diff%>%
  filter(n==1)%>%
  mutate(trt_type2="n")

CT_Sp_P<-CT_diff%>%
  filter(p==1)%>%
  mutate(trt_type2="p")

CT_Sp_nuts_other<-CT_diff%>%
  filter(nuts_other==1)%>%
  mutate(trt_type2="nuts_other")

CT_Sp_allint<-CT_diff%>%
  filter(multtrts==1)%>%
  mutate(trt_type2="all mult")

Fulldataset_mixedmodels<-CT_Sp_allint%>%
  bind_rows(CT_Sp_co2, CT_Sp_dist, CT_Sp_drt, CT_Sp_herb, CT_Sp_irg, CT_Sp_N, CT_Sp_P, CT_Sp_temp)%>%
  select(site_code, project_name, community_type, species_matched, trt_type2, diff)

write.csv(Fulldataset_mixedmodels, paste(my.wd, "WinnersLosers paper/data/Species_DCiDiff_formixedmodels.csv", sep=""), row.names=F)
