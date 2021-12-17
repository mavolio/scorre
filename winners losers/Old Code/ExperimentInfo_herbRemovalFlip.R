setwd("C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\CleanedData\\Sites\\Species csv")

library(tidyverse)


# notes:
# nutrients, light, carbon, water, and other are binary variables for the entire experiment (whether one of these factors was manipulated)
# n-other_trt are variables describing the specific treatment amounts or categories
#### if all plots were burned (at any frequency (<20yrs) over the course of the experiment or at the site where the experiment takes place), fenced to remove herbivores (at any point in time where the experiment takes place), or mowed/clipped, then the variable gets a 0 for all plots, but management=1. SEE Experiment_List for Details. Otherwise they match treatments.
##herb_removal is fencing to remove herbivores, EXCEPT for MAERC fireplots, where it is grazers being added
#### successional and plant_mani are binary variables
####   if successional=1, then all plots were disturbed in some way prior to the experiment start
#   if plant_mani=1, then some or all species are planted into the plots at the start of the experiment; this can vary by treatment within a site (but most sites it is a 1 for all plots)
# plant_trt=1, the treatment manipulates the plants directly in some way; plant_trt is a binary variable for experiments that manipulate plants as a treatment (e.g., seed addition, removal of some species); that is, manipulation of the plants for some treatments but not others (not e.g., seeding all plots at the start of the experiment)
#### cessation is a binary variable for each treatment and year combination, and is 1 when the treatment has been stopped - DROPPING CESSATION
# plot_mani is the total number of factors manipulated compared to the control (i.e., if all plots are burned plot_mani does not increase, but if treatment plots are burned when controls are not, then plot_mani increases by 1)
# resource_mani is a binary variable to subset out the treatments that do not directly manipulate a resource (e.g., only increased temperature)
#   resource_mani=1 for all control plots and any treatment plot that directly manipulates a resource; resource_mani=0 for plots where no resource was directly manipulated (except for the controls)
# pulse is a binary variable for treatments that were a one time application that did not also get applied to the control plots (e.g., burning in CUL, which was a treatment just in the first year of the experiment, but was not applied to the controls)
#   a factor was not considered a pulse if all plots experienced it (e.g., solarization in ASGA clonal, wildfire in SEV WENNDex)
# public means the dataset is publically available and we don't need to ask permission to use it
# max_trt means that the treatments are the maximum magnitude of the treatment or controls (only applicable for experiments that have multiple levels of the same treatment, e.g., CDR e001 has many levels of N addition)
# factorial means that the treatments are factorially manipulated - can only be factorial with 2+ treatments, some experiments have some treatments that are factorially manipulated and others that are not, must have all levels (i.e., can't have 3 and 4 combos without 1 and 2 combos)
# trt_type is a categorical description of the treatments: mult_nutrient category is anything more than N*P (i.e., N*P is listed as a seperate category)

watering<-read.delim("ANG_watering.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='W', 20, ifelse(treatment=='S', 20, 0)), 
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=ifelse(treatment=='W', 'winter water addition', ifelse(treatment=='S', 'spring water addition', 0)), 
         trt_details=0,
         successional=0, 
         plant_mani=0,
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment %in% c('W','S'), 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment %in% c('W','S'), 'irr', 'control'))%>%
  unique()

fert1<-read.csv("ANR_Fert1.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n= ifelse(treatment %in% c('KNO3','NH4PO4', 'NH4NO3', 'full_nut'), 5, 0), 
         p= ifelse(treatment == 'full_nut', 1, 0), 
         k= ifelse(treatment == 'full_nut', 2.8,0), 
         CO2=0,
         precip=0, 
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=ifelse(treatment=='ACTIVATED CARBON', 'activated carbon addition', ifelse(treatment=='GLUCOS', 'glucose addition',
                   ifelse(treatment=='PROTEIN', 'BAS organic N addition', ifelse(treatment=='CACO3', 'lime addition',ifelse(treatment=='micronut', '0.25 g m2 of micronutrient solution added',0))))), 
         trt_details=0,
         successional=0, 
         plant_mani=ifelse(treatment == 'REDUCTION', 1, 0),
         plant_trt= ifelse(treatment == 'REDUCTION', 1, 0),
         pulse= ifelse(treatment == 'REDUCTION', 1, 0)) %>%
  mutate(plot_mani=ifelse(treatment == 'control', 0, 1))%>%
  mutate(resource_mani= ifelse(treatment == "REDUCTION", 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment == 'control', 'control', ifelse(treatment == 'REDUCTION', 'plant_mani', ifelse(treatment %in% c('micronut', 'full_nut'), 'mult_nutrient', ifelse(treatment == "ACTIVATED CARBON", 'C', ifelse(treatment == 'CACO3', 'lime', ifelse(treatment == 'PROTEIN', 'protein', 'N')))))))%>%
  unique()

fert2<-read.csv("ANR_Fert2.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n= ifelse(treatment %in% c("full_nut_Dflex_removal","full_nut_full_removal","full_nut_no_removal",
                                    "NH4NO3_Dflex_removal", "NH4NO3_full_removal", "NH4NO3_no_removal",
                                    "KNO3_Dflex_removal", "KNO3_full_removal", "KNO3_no_removal", 
                                    "NH4PO4_Dflex_removal", "NH4PO4_full_removal", "NH4PO4_no_removal",
                                    "full_nut_Eherm_removal", "NH4NO3_Eherm_removal", "KNO3_Eherm_removal",
                                    "NH4PO4_Eherm_removal"), 5, 0), 
         p= ifelse(treatment %in% c('full_nut_Dflex_removal', "full_nut_full_removal","full_nut_no_removal", "full_nut_Eherm_removal"), 1, 0), 
         k= ifelse(treatment %in% c('full_nut_Dflex_removal', "full_nut_full_removal","full_nut_no_removal", "full_nut_Eherm_removal"), 2.8, 0), 
         CO2=0,
         precip=0, 
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt= 0, 
         trt_details=0, 
         successional=0, 
         plant_mani=ifelse(treatment %in% c('full_nut_no_removal', 'NH4NO3_no_removal', 'KNO3_no_removal', 'NH4PO4_no_removal', 'control_no_removal'), 0 , 1),
         plant_trt= ifelse(treatment %in% c('full_nut_no_removal', 'NH4NO3_no_removal', 'KNO3_no_removal', 'NH4PO4_no_removal', 'control_no_removal'), 0 , 1),
         pulse= ifelse(treatment %in% c('full_nut_full_removal', 'NH4NO3_full_removal', 'KNO3_full_removal', 'NH4PO4_full_removal', 'control_full_removal'), 1, 0)) %>%
  mutate(plot_mani=ifelse(treatment == 'control_no_removal', 0,ifelse(treatment %in% c('control_Dflex_removal', 'control_Eherm_removal', 'control_full_removal'), 1, 2)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt= ifelse(treatment %in% c("full_nut_full_removal", "NH4NO3_full_removal", "KNO3_full_removal", "NH4PO4_full_removal", 'control_full_removal','control_no_removal'),1,0))%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment == 'control_no_removal', 'control', ifelse(treatment %in% c("full_nut_Eherm_removal", 'full_nut_Dflex_removal', 'full_nut_full_removal'), 'mult_nutrient*plant_mani',
                  ifelse(treatment %in% c("control_Eherm_removal", 'control_Dflex_removal', 'control_full_removal'), 'plant_mani', 'N*plant_mani'))))%>%
  unique()


mat2<-read.delim("ARC_mat2.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='NP', 10, 0), 
         p=ifelse(treatment=='NP', 5, 0), 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='NP', 2, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='NP', 'N*P', 'control'))%>%
  unique()

mnt<-read.delim("ARC_mnt.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='NP', 10, 0), 
         p=ifelse(treatment=='NP', 5, 0), 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0,
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='NP', 2, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='NP', 'N*P', 'control'))%>%
  unique()

clonal<-read.delim("ASGA_Clonal.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='mixed_CO', 0, ifelse(treatment=='non-clonal_CO', 0, 20.1)), 
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0,
         trt_details=ifelse(treatment=='non-clonal_CO', 'non-clonal species', ifelse(treatment=='mixed_LP', 'large nutrient patches', ifelse(treatment=='non-clonal_LP', 'non-clonal species, large nutrient patches', ifelse(treatment=='mixed_SP', 'small nutrient patches', ifelse(treatment=='non-clonal_SP', 'non-clonal species, small nutrient patches', ifelse(treatment=='non-clonal_UN', 'non-clonal species', 0)))))),
         successional=1, 
         plant_mani=1, 
         plant_trt=ifelse(treatment %in% c('non-clonal_CO','non-clonal_LP','non-clonal_SP','non-clonal_UN'), 1, 0),
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment %in% c('non-clonal_UN','non-clonal_LP','non-clonal_SP'), 2, ifelse(treatment %in% c('non-clonal_CO','mixed_LP','mixed_SP','mixed_UN'), 1, 0)))%>%
  mutate(resource_mani=ifelse(treatment=='non-clonal_CO', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='non-clonal_CO', 'plant_mani', ifelse(treatment %in% c('mixed_LP','mixed_SP','mixed_UN'), 'N', ifelse(treatment %in% c('non-clonal_LP','non-clonal_SP','non-clonal_UN'), 'N*plant_mani', 'control'))))%>%
  unique()

exp1<-read.delim("ASGA_Exp1.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='2_0_CO'|treatment=='1_0_CO'|treatment=='2_1_CO'|treatment=='1_1_CO', 0, 20.1),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0, 
         management=0,
         other_trt=0,
         trt_details=ifelse(treatment %in% c('2_0_PA','1_0_PA','2_1_PA','1_1_PA'), 'nutrient patches', 0),
         successional=ifelse(treatment %in% c('2_0_CO','2_0_PA','2_0_UN','2_1_CO','2_1_PA','2_1_UN'), 1, 0),
         plant_mani=ifelse(treatment %in% c('2_1_CO','2_1_PA','2_1_UN','1_1_CO','1_1_PA','1_1_UN'), 1, 0),
         plant_trt=ifelse(treatment %in% c('2_1_CO','2_1_PA','2_1_UN','1_1_CO','1_1_PA','1_1_UN'), 1, 0),
         pulse=ifelse(treatment %in% c('2_0_PA','2_0_UN','2_0_CO','2_1_PA','2_1_UN','2_1_CO'), 1, 0))%>%
  mutate(plot_mani=ifelse(treatment=='1_0_CO', 0, ifelse(treatment %in% c('1_0_PA','1_0_UN','1_1_CO','2_0_CO'), 1, ifelse(treatment %in% c('1_1_PA','1_1_UN','2_0_PA','2_0_UN','2_1_CO'), 2, 3))))%>%
  mutate(resource_mani=ifelse(treatment %in% c('2_0_CO','1_1_CO','2_1_CO'), 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment %in% c('1_0_PA','1_0_UN'), 'N', ifelse(treatment %in% c('2_0_PA','2_0_UN'), 'N*disturbance', ifelse(treatment %in% c('1_1_PA','1_1_UN'), 'N*plant_mani', ifelse(treatment=='2_1_CO', 'plant_mani*disturbance', ifelse(treatment %in% c('2_1_PA','2_1_UN'), 'N*plant_mani*disturbance', ifelse(treatment=='2_0_CO', 'disturbance', ifelse(treatment=='1_1_CO', 'plant_mani', 'control'))))))))%>%
  unique()

eelplot <- read.csv("AZI_EELplot.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1, 
         n=ifelse(treatment %in% c("N", "F+N", "F+N+W", "N+W"), 12, 0), 
         p= 0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp= ifelse(treatment %in% c("W", "F+N+W", "N+W", "F+W"), 1.47, 0), 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt= ifelse(treatment %in% c("F", "F+N", "F+W+N", "F+W"), 'fungicide', 0),
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment %in% c("N", "F", "W"), 1, ifelse(treatment == "F+N+W", 3, 2))))%>%
  mutate(resource_mani=ifelse(treatment == "F", 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='C', 'control', ifelse(treatment == "N", 'N', ifelse(treatment == "W", 'temp', ifelse(treatment == "F", 'fungicide',
                  ifelse(treatment == "F+N", "N*fungicide", ifelse(treatment == "N+W", "N*temp", ifelse(treatment == "F+W", "temp*fungicide", "N*temp*fungicide"))))))))%>%
  unique()


nitphos<-read.csv("AZI_NitPhos.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0, 
         n=ifelse(treatment=='N1P0', 5, ifelse(treatment=='N2P0', 10, ifelse(treatment=="N3P0", 15, ifelse(treatment=="N0P0", 0, 10)))), 
         p=ifelse(treatment=='N2P1', 2, ifelse(treatment=="N2P2", 4, ifelse(treatment=="N2P3", 8, 0))), 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0,
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N0P0', 0, ifelse(treatment=="N1P0", 1, ifelse(treatment=="N2P0", 1, ifelse(treatment=="N3P0", 1, 2)))))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=ifelse(treatment %in% c('N0P0','N2P0','N2P3','N3P0'), 1, 0))%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='N0P0', 'control', ifelse(treatment %in% c('N1P0','N2P0','N3P0'), 'N', 'N*P')))%>%
  unique()

#16 spp plots are controls
lind<-read.delim("BAY_LIND.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment %in% c('rain_rich1','rain_rich2','rain_rich4','rain_rich8','rain_rich16'), 8, 0),
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0, 
         management=1,
         pulse=0,
         other_trt=0,
         trt_details=ifelse(treatment %in% c('ref_rich1','rain_rich1'), '1 sp', ifelse(treatment %in% c('ref_rich2','rain_rich2'), '2 sp', ifelse(treatment %in% c('ref_rich4','rain_rich4'), '4 sp', ifelse(treatment %in% c('ref_rich8','rain_rich8'), '8 sp', '16 sp')))),
         successional=1, 
         plant_mani=1,  
         plant_trt=ifelse(treatment %in% c('ref_rich16','rain_rich16'), 0, 1),
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='ref_rich16', 0, ifelse(treatment %in% c('ref_rich1','ref_rich2','ref_rich4','ref_rich8','rain_rich16'), 1, 2)))%>%
  mutate(resource_mani=ifelse(treatment %in% c('ref_rich1','ref_rich2','ref_rich4','ref_rich8'), 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='ref_rich16', 'control', ifelse(treatment %in% c('ref_rich1','ref_rich2','ref_rich4','ref_rich8'), 'plant_mani', 'irr*plant_mani')))%>%
  unique()

events<-read.delim("Bt_EVENT2.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=ifelse(treatment=='CM-N1', 'reduced precip variability', ifelse(treatment=='D1-N1', 'early drought', ifelse(treatment=='D2-N1', 'late drought', 0))), 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='CA-N1', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='CA-N1', 'control', 'precip_vari'))%>%
  unique()

pq<-read.delim("BUX_PQ.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=0, 
         precip=ifelse(treatment %in% c('wet','warm wet'), 20, ifelse(treatment %in% c('dry','warm dry'), -20, 0)),
         temp=ifelse(treatment %in% c('warm','warm dry','warm wet'), 3, 0),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment %in% c('warm','dry','wet'), 1, ifelse(treatment=='control', 0, 2)))%>%
  mutate(resource_mani=ifelse(treatment=='warm', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='dry', 'drought', ifelse(treatment=='wet', 'irr', ifelse(treatment=='warm', 'temp', ifelse(treatment=='warm dry', 'drought*temp', ifelse(treatment=='warm wet', 'irr*temp', 'control'))))))%>%
  unique()

pennings<-read.delim("CAR_Pennings.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='NPK'&calendar_year>1999, 164, ifelse(treatment=='NPK'&calendar_year==1999, 84, 0)),
         p=ifelse(treatment=='NPK'&calendar_year>1999, 82, ifelse(treatment=='NPK'&calendar_year==1999, 42, 0)),
         k=ifelse(treatment=='NPK'&calendar_year>1999, 41, ifelse(treatment=='NPK'&calendar_year==1999, 21, 0)),
         CO2=0, 
         precip=0,
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='NPK', 3, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='NPK', 'mult_nutrient', 'control'))%>%
  unique()

rmapc<-read.delim("CAU_RMAPC.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='N', 9, ifelse(treatment=='NP', 9, 0)),
         p=ifelse(treatment=='P', 2.6, ifelse(treatment=='NP', 2.6, 0)),
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='H2O', 20, 0),
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=ifelse(treatment=='Ca', "lime added", 0), 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='Cont', 0, ifelse(treatment=='NP', 2, 1)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=ifelse(treatment=='H2O', 0, ifelse(treatment=='Ca', 0, 1)))%>%
  mutate(trt_type=ifelse(treatment=='H2O', 'irr', ifelse(treatment=='N', 'N', ifelse(treatment=='P', 'P', ifelse(treatment=='Ca', 'lime', ifelse(treatment=='NP', 'N*P', 'control'))))))%>%
  unique()

biocon<-read.csv("CDR_BioCON.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=1, water=0, other_manipulation=0,
         n=ifelse(treatment=='Camb_Nenrich', 4, ifelse(treatment=='Cenrich_Nenrich', 4, 0)),
         p=0, 
         k=0, 
         CO2=ifelse(treatment=='Cenrich_Namb', 160, ifelse(treatment=='Cenrich_Nenrich', 160, 0)),
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0,
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=1, 
         plant_mani=1,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='Camb_Namb', 0, ifelse(treatment=='Cenrich_Nenrich', 2, 1)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='Namb_Cenrich', 'CO2', ifelse(treatment=='Camb_Nenrich', 'N', ifelse(treatment=='Cenrich_Nenrich', 'N*CO2', 'control'))))%>%
  unique()

e001<-read.csv("CDR_e001.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(treatment=as.factor(treatment))%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='2', 1.02, ifelse(treatment=='3', 2.04, ifelse(treatment=='4', 3.40, ifelse(treatment=='5', 5.44, ifelse(treatment=='6', 9.52, ifelse(treatment=='7', 17, ifelse(treatment=='8', 27.2, 0))))))),
         p=ifelse(treatment=='9', 0, 4.6), 
         k=ifelse(treatment=='9', 0, 6.1),
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=ifelse(treatment=='9', 0, "micronutrients and lime added"), 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='1', 4, ifelse(treatment=='9', 0, 5)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=ifelse(treatment %in% c('8','1','9'), 1, 0))%>%
  mutate(public=1)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='9', 'control', 'mult_nutrient'))%>%
  unique()

e002<-read.delim("CDR_e002.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='2_f_u_n', 1.02, ifelse(treatment=='3_f_u_n', 2.04, ifelse(treatment=='4_f_u_n', 3.4, ifelse(treatment=='5_f_u_n', 5.44, ifelse(treatment=='6_f_u_n', 9.52, ifelse(treatment=='7_f_u_n', 17, ifelse(treatment=='8_f_u_n', 27.2, 0))))))),
         p=ifelse(treatment=='9_f_u_n', 0, 4.6),
         k=ifelse(treatment=='9_f_u_n', 0, 6.1),
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=ifelse(treatment=='9_f_u_n', 0, "micronutrients and lime added"),
         trt_details=0,
         successional=1, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='1_f_u_n', 4, ifelse(treatment=='9_f_u_n', 0, 5)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  mutate(max_trt=ifelse(treatment %in% c('1_f_u_n','8_f_u_n','9_f_u_n'), 1, 0))%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='9_f_u_n', 'control', 'mult_nutrient'))%>%
  filter(calendar_year<1992)%>%##drops everything once cessation starts
  unique()

megarich<-read.delim("CEH_Megarich.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=1, water=0, other_manipulation=1,
         n=10, 
         p=2, 
         k=2,
         CO2=ifelse(treatment=='EcAt', 280, ifelse(treatment=='EcEt', 280, 0)),
         precip=0,
         temp=ifelse(treatment=='AcEt', 2.9, ifelse(treatment=='EcEt', 2.9, 0)),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='AcAt', 0, ifelse(treatment=='EcEt', 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=='AcEt', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='AcEt', 'temp', ifelse(treatment=='EcAt', 'CO2', ifelse(treatment=='EcEt', 'CO2*temp', 'control'))))%>%
  unique()

imagine<-read.delim("CLE_imagine.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=1, water=1, other_manipulation=1,
         n=0,
         p=0, 
         k=0,
         CO2=ifelse(treatment=='TDCO2', 200, 0),
         precip=ifelse(treatment=='TD', -20, ifelse(treatment=='TDCO2', -20, 0)),
         temp=ifelse(treatment=='C', 0, 3.5),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment=='T', 1, ifelse(treatment=='TD', 2, 3))))%>%
  mutate(resource_mani=ifelse(treatment=='T', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='T', 'temp', ifelse(treatment=='TD', 'drought*temp', ifelse (treatment=='TDCO2', 'drought*CO2*temp', 'control'))))%>%
  unique()

culardoch<-read.delim("CUL_culardoch.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('N10','N10burn','N10clip','N10burnclip'), 1, ifelse(treatment %in% c('N20','N20burn','N20clip','N20burnclip'), 2, ifelse(treatment %in% c('N50','N50burn','N50clip','N50burnclip'), 5, 0))),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=ifelse(treatment %in% c('clip','burnclip','N10clip','N20clip','N50clip','N10burnclip','N20burnclip','N50burnclip'), 1, 0),
         burn=ifelse(treatment %in% c('N10burn','N20burn','N50burn','burn','burnclip','N10burnclip','N20burnclip','N50burnclip'), 1, 0),
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=ifelse(treatment %in% c('burn','burnclip','N10burn','N10burnclip','N20burn','N20burnclip','N50burn','N50burnclip'), 1, 0))%>%
  mutate(plot_mani=ifelse(treatment=='control', 0, ifelse(treatment %in% c('N10','N20','N50','burn','clip'), 1, ifelse(treatment %in% c('N10burnclip','N20burnclip','N50burnclip'), 3, 2))))%>%
  mutate(resource_mani=ifelse(treatment %in% c('burn','clip','burnclip'), 0, 1))%>%
  mutate(max_trt=ifelse(treatment %in% c('control','N50','burn','N50burn','clip','burnclip','N50clip','N50burnclip'), 1, 0))%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='control', 'control', ifelse(treatment=='burn', 'burn', ifelse(treatment=='clip', 'mow_clip', ifelse(treatment=='burnclip', 'burn*mow_clip', ifelse(treatment %in% c('N10','N20','N50'), 'N', ifelse(treatment %in% c('N10burn','N20burn','N50burn'), 'N*burn', ifelse(treatment %in% c('N10clip','N20clip','N50clip'), 'N*mow_clip', 'N*burn*mow_clip'))))))))%>%
  unique()

gap2<-read.csv("DCGS_gap.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=1, carbon=0, water=0, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0,
         burn=0, 
         herb_removal=0, 
         management=0,
         other_trt=ifelse(treatment=='_018', '18 ft opening', ifelse(treatment=='_033', '33 ft opening', ifelse(treatment=='_066', '66 ft opening', ifelse(treatment=='_100', '100 ft opening', ifelse(treatment=='_150', '150 ft opening', 0))))),
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=1)%>%
  mutate(plot_mani=ifelse(treatment=='_000', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=ifelse(treatment=='_000'|treatment=='_150', 1, 0))%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='_000', 'control', 'light'))%>%
  unique()

gcme<- read.csv("DCMIC_GCME.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1, 
         n=ifelse(treatment %in% c("N", "MNP", "NP", "MN"), 10, 0), 
         p=ifelse(treatment %in% c("MP", "MNP", "NP", "P"), 5, 0),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0, 
         mow_clip= ifelse(treatment %in% c("MP", "MNP", "M", "MN"), 1, 0), 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0,
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment %in% c("N", "M", "P"),1, ifelse(treatment=="MNP", 3, 2))))%>%
  mutate(resource_mani= ifelse(treatment == "M", 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment == "C", "control", ifelse(treatment == "N", "N", ifelse(treatment == "M", 'mow_clip', 
                  ifelse(treatment == "P", "P", ifelse(treatment == "MP", "P*mow_clip", ifelse(treatment == "MNP", "N*P*mow_clip",
                  ifelse(treatment == "NP", "N*P", "N*mow_clip"))))))))%>%
  unique()

gcme2<- read.csv("DCMIC_GCME2.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=1, 
         n=ifelse(treatment %in% c("N", "MNP", "NP", "MN"), 10, 0), 
         p=0,
         k=0, 
         CO2=0, 
         precip=ifelse(treatment %in% c("MP", "MNP", "NP", "P"), 120, 0), 
         temp=0, 
         mow_clip= ifelse(treatment %in% c("MP", "MNP", "M", "MN"), 1, 0), 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0,
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment %in% c("N", "M", "P"),1, ifelse(treatment=="MNP", 3, 2))))%>%
  mutate(resource_mani=ifelse(treatment == "M", 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment == "C", "control", ifelse(treatment == "N", "N", ifelse(treatment == "M", 'mow_clip', 
                  ifelse(treatment == "P", "irr", ifelse(treatment == "MP", "irr*mow_clip", ifelse(treatment == "MNP", "N*irr*mow_clip",
                  ifelse(treatment == "NP", "N*irr", "N*mow_clip"))))))))%>%
  unique()

precip<- read.csv("DCMIC_Precip.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=0, 
         n=0, 
         p=0,
         k=0, 
         CO2=0, 
         precip=ifelse(treatment == "P-6", -60, ifelse(treatment == "P-4", -40, ifelse(treatment =="P-2", -20,
                ifelse(treatment == "CK", 0, ifelse(treatment == "P+2", 20, ifelse(treatment == "P+4", 40, 60)))))), 
         temp=0, 
         mow_clip= 0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0,
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='CK', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt= ifelse(treatment %in% c("P-6", "P+6", "CK"), 1, 0))%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment == "CK", "control", ifelse(treatment %in% c("P-6", "P-4", "P-2"), "drought", "irr")))%>%
  unique()

nsfc<-read.delim("DL_NSFC.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='N', 10, ifelse(treatment=='WN', 10, 0)),
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='W', 49.8, ifelse(treatment=='WN', 49.8, 0)),
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment=='WN', 2, 1)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='N', 'N', ifelse(treatment=='W', 'irr', ifelse(treatment=='WN', 'N*irr', 'control'))))%>%
  unique()
nsfc2<-read.csv("DL_NSFC20132016.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='N10', 10, ifelse(treatment=='WN10', 10, 0)),
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='WCK', 49.8, ifelse(treatment=='WN10', 49.8, 0)),
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='CK', 0, ifelse(treatment=='WN10', 2, 1)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='N10', 'N', ifelse(treatment=='WCK', 'irr', ifelse(treatment=='WN10', 'N*irr', 'control'))))%>%
  unique()

warmnut<-read.delim("Finse_WarmNut.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('nutrient addition','warming + nutrient addition'), 10, 0),
         p=ifelse(treatment %in% c('nutrient addition','warming + nutrient addition'), 2, 0),
         k=ifelse(treatment %in% c('nutrient addition','warming + nutrient addition'), 8, 0),
         CO2=0, 
         precip=0,
         temp=ifelse(treatment %in% c('warming','warming + nutrient addition'), 1.5, 0),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='control', 0, ifelse(treatment=='warming', 1, ifelse(treatment=='nutrient addition', 3, 4))))%>%
  mutate(resource_mani=ifelse(treatment=='warming', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='nutrient addition', 'mult_nutrient', ifelse(treatment=='warming', 'temp', ifelse(treatment=='warming + nutrient addition', 'mult_nutrient*temp', 'control'))))%>%
  unique()

gfert<- read.csv("Glen_Fert.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0, 
         n=ifelse(treatment %in% c("N", "NP"), 7.5, 0), 
         p=ifelse(treatment %in% c("P", "NP"), 5, 0),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0, 
         mow_clip= 0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0,
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment == "NP", 2, 1)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment == "C", "control", ifelse(treatment == "N", "N", ifelse(treatment == "P", "P", "NP"))))%>%
  unique()

face<-read.delim("GVN_FACE.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=1, water=0, other_manipulation=0,
         n=0, 
         p=0, 
         k=0,
         CO2=ifelse(treatment=='A', 0, 160),
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='E', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='A', 'CO2', 'control'))%>%
  unique()

warmnit <- read.csv("Hayoka_WarmNit.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c("Nitrogen", "Winter warming + Nitrogen"), 0.5, 0), 
         p=0, 
         k=0,
         CO2=0,
         precip=0, 
         temp=ifelse(treatment %in% c("Winter warming", "Winter warming + Nitrogen"), 5, 0),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment == "Control", 0, ifelse(treatment == "Winter warming + Nitrogen", 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment == "Winter warming", 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=="Control", "control", ifelse(treatment == "Nitrogen", "N", ifelse(treatment == "Winter warming", "temp", "N*temp"))))%>%
  unique()

hprecip <- read.csv("HAYS_Precip.csv") %>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, 
         p=0, 
         k=0,
         CO2=0,
         precip= ifelse(treatment == "reduction", -50, ifelse(treatment == "add", 61, 0)), 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment == "Control", 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=="Control", "control", ifelse(treatment == "reduction", "drought", "irr"))) %>%
  unique()

phace <- read.csv("HPGRS_PHACE.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=1, water=0, other_manipulation=1,
         n=0, 
         p=0, 
         k=0,
         CO2= ifelse(treatment %in% c("Ct", "CT"), 600,0),
         precip=0, 
         temp=ifelse(treatment %in% c("cT", "CT"), 2, 0), # 1.5 in the day 3 at night, took average of these
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment == "ct", 0, ifelse(treatment == "CT", 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment == "cT", 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=="ct", "control", ifelse(treatment == "Ct", "CO2", ifelse(treatment == "cT", "temp", "CO2*temp")))) %>%
  unique()

nde<-read.csv("IMGERS_NDE.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c("N1M0","N1M1"), 1,ifelse(treatment %in% c("N2M0","N2M1"), 2, ifelse(treatment %in% c("N3M0","N3M1"), 3, ifelse(treatment %in% c("N4M0","N4M1"),5, ifelse(treatment %in% c("N5M0","N5M1"), 10, ifelse(treatment %in% c("N6M0","N6M1"),15, ifelse(treatment %in% c("N7M0","N7M1"), 20, ifelse(treatment %in% c("N8M0","N8M1"),50,0)))))))), 
         p=0, 
         k=0, 
         CO2=0,
         precip=0, 
         temp=0, 
         mow_clip=ifelse(treatment %in% c("N0M1","N1M1","N2M1","N3M1","N4M1","N5M1","N6M1","N7M1","N8M1"), 1,0), 
         burn=0, 
         herb_removal=0, 
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=="N0M0", 0, ifelse(treatment %in% c("N1M0",'N0M1',"N1M0","N2M0","N3M0","N4M0","N5M0","N6M0","N7M0","N8M0"), 1, 2)))%>%
  mutate(resource_mani=ifelse(treatment=="N0M1", 0, 1))%>%
  mutate(max_trt=ifelse(treatment %in% c('N0M0','N8M0','N0M1','N8M1'), 1, 0))%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='N0M0', 'control', ifelse(treatment=='N0M1', 'mow_clip', ifelse(treatment %in% c('N1M0','N2M0','N3M0','N4M0','N5M0','N6M0','N7M0','N8M0'), 'N', 'N*mow_clip'))))%>%
  unique()

yu<-read.delim("IMGERS_Yu.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N2', 5.6, ifelse(treatment=='N3', 11.2, ifelse(treatment=='N4', 22.4, ifelse(treatment=='N5', 39.2, ifelse(treatment=='N6', 56, 0))))),
         p=ifelse(treatment=='N0', 0, 1.55),
         k=ifelse(treatment=='N0', 0, 3.95),
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N0', 0, ifelse(treatment=='N1', 2, 3)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=ifelse(treatment=='N0'|treatment=='N1'|treatment=='N6', 1, 0))%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='N0', 'control', 'mult_nutrient'))%>%
  unique()

study119<-read.delim("JRN_Study119.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='T', 10, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='T', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=0)%>%
  filter(calendar_year<1986)%>%
  mutate(trt_type=ifelse(treatment=='T', 'N', 'control'))%>%
  unique()

study278<-read.delim("JRN_study278.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment %in% c('P1N1','P2N1','P3N1','P4N1','P5N1'), 10, 0),
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment %in% c('P1N0','P1N1'), -80, ifelse(treatment %in% c('P2N0','P2N1'), -50, ifelse(treatment %in% c('P4N0','P4N1'), 50, ifelse(treatment %in% c('P5N0','P5N1'), 80, 0)))),
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='P3N0', 0, ifelse(treatment %in% c('P1N0','P2N0','P3N1','P4N0','P5N0'), 1, 2)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=ifelse(treatment %in% c('P1N0','P1N1','P3N0','P3N1','P5N0','P5N1'), 1, 0))%>%
  mutate(public=1)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='P3N0', 'control', ifelse(treatment %in% c('P1N0','P2N0'), 'drought', ifelse(treatment %in% c('P4N0','P5N0'), 'irr', ifelse(treatment=='P3N1', 'N', ifelse(treatment %in% c('P1N1','P2N1'), 'N*drought', 'N*irr'))))))%>%
  unique()

gce<-read.delim("JSP_GCE2.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=1, water=1, other_manipulation=1,
         n=ifelse(treatment %in% c('N','RN','HN','HRN','CN','CRN','CHN','CHRN'), 7, 0),
         p=0, 
         k=0, 
         CO2=ifelse(treatment %in% c('C','CN','CR','CRN','CH','CHN','CHR','CHRN'), 300, 0),
         precip=ifelse(treatment %in% c('R','RN','HR','HRN','CR','CRN','CHR','CHRN'), 50, 0),
         temp=ifelse(treatment %in% c('H','HN','HR','HRN','CH','CHN','CHR','CHRN'), 1.5, 0),
         mow_clip=0,
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0,
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='amb', 0, ifelse(treatment %in% c('N','R','H','C'), 1, ifelse(treatment %in% c('HRN','CRN','CHN','CHR'), 3, ifelse(treatment=='CHRN', 4, 2)))))%>%
  mutate(resource_mani=ifelse(treatment=='H', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='amb', 'control', ifelse(treatment=='C', 'CO2', ifelse(treatment=='R', 'irr', ifelse(treatment=='N', 'N', ifelse(treatment=='H', 'temp', ifelse(treatment=='CR', 'irr*CO2', ifelse(treatment=='CHR', 'irr*CO2*temp', ifelse(treatment=='HR', 'irr*temp', ifelse(treatment=='CH', 'CO2*temp', ifelse(treatment=='CN', 'N*CO2', ifelse(treatment=='CHN', 'N*CO2*temp', ifelse(treatment=='RN', 'N*irr', ifelse(treatment=='CRN', 'N*irr*CO2', ifelse(treatment=='HN', 'N*temp', ifelse(treatment=='HRN', 'N*irr*temp', 'N*irr*CO2*temp'))))))))))))))))%>%
  unique()

wapaclip<-read.delim("KAEFS_WaPaClip.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment %in% c('U CH','U WH','C CH','C WH'), -50, ifelse(treatment %in% c('U CD','U WD','C CD','C WD'), 50, 0)),
         temp=ifelse(treatment %in% c('U WC','U WH','U WD','C WC','C WH','C WD'), 3, 0),
         mow_clip=ifelse(treatment %in% c('C CC','C CH','C CD','C WC','C WH','C WD'), 1, 0),
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='U CC', 0, ifelse(treatment %in% c('U CH','U CD','U WC','C CC'), 1, ifelse(treatment=='C WH', 3, ifelse(treatment=='C WD', 3, 2)))))%>%
  mutate(resource_mani=ifelse(treatment %in% c('U WC','C CC','C WC'), 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='C CC', 'mow_clip', ifelse(treatment=='C CD', 'irr*mow_clip', ifelse(treatment=='C CH', 'drought*mow_clip', ifelse(treatment=='C WC', 'temp*mow_clip', ifelse(treatment=='C WD', 'irr*temp*mow_clip', ifelse(treatment=='C WH', 'drought*temp*mow_clip', ifelse(treatment=='U CC', 'control', ifelse(treatment=='U CD', 'irr', ifelse(treatment=='U CH', 'drought', ifelse(treatment=='U WC', 'temp', ifelse(treatment=='U WD', 'irr*temp', 'drought*temp'))))))))))))%>%
  unique()

t7<-read.delim("KBS_T7.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('T0F1','T1F1'), 12.3, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0, 
         management=1,
         other_trt=ifelse(treatment %in% c('T1F0','T1F1'), 'tilled', 0),
         trt_details=0,
         successional=1, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='T0F0', 0, ifelse(treatment=='T1F1', 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=='T1F0', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='T0F0', 'control', ifelse(treatment=='T1F1', 'N*till', ifelse(treatment=='T0F1', 'N', 'till'))))%>%
  unique()

bffert<-read.delim("KLU_BFFert.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('N1F0','N1F1'), 17.5, 0),
         p=ifelse(treatment %in% c('N1F0','N1F1'), 5, 0),
         k=ifelse(treatment %in% c('N1F0','N1F1'), 1.5, 0),
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0,
         herb_removal=ifelse(treatment %in% c('N0F1','N1F1'), 1, 0),
         management=0,
         trt_details=0,
         other_trt=0, 
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(resource_mani=ifelse(treatment %in% c('N0F0','N0F1'), 0, 1))%>%
  mutate(plot_mani=ifelse(treatment=='N0F1', 0, ifelse(treatment=='N1F1',3, ifelse(treatment=='N1F0', 4, 1))))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='N0F1', 'control', ifelse(treatment=='N0F0', 'grazed', ifelse(treatment=='N1F1', 'mult_nutrient', 'mult_nutrient*grazed'))))%>%
  unique()

kgfert<-read.delim("KLU_KGFert.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('N1B0','N1B1'), 17.5, 0),
         p=ifelse(treatment %in% c('N1B0','N1B1'), 5.8, 0),
         k=ifelse(treatment %in% c('N1B0','N1B1'), 5.8, 0),
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0, 
         management=0,
         other_trt=ifelse(treatment %in% c('N0B1','N1B1'), "fungicide added", 0), 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N0B0', 0, ifelse(treatment=='N1B1', 4, ifelse(treatment=='N0B1', 1, 3))))%>%
  mutate(resource_mani=ifelse(treatment=='N1B0',1, ifelse(treatment=='N1B1', 1, 0)))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='N0B0', 'control', ifelse(treatment=='N0B1', 'fungicide', ifelse(treatment=='N1B0', 'mult_nutrient', 'mult_nutrient*fungicide'))))%>%
  unique()

bgp<-read.csv("KNZ_BGP.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('u_u_p','u_u_c','u_m_p','u_m_c','b_u_p','b_u_c','b_m_p','b_m_c'), 0, 10),
         p=ifelse(treatment %in% c('u_u_n','u_u_c','u_m_n','u_m_c','b_u_n','b_u_c','b_m_n','b_m_c'), 0, 1),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=ifelse(treatment %in% c('u_u_n','u_u_p','u_u_c','u_u_b','b_u_n','b_u_p','b_u_c','b_u_b'), 0, 1),
         burn=ifelse(treatment %in% c('u_u_n','u_u_p','u_u_c','u_u_b','u_m_n','u_m_p','u_m_c','u_m_b'), 0, 1),
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='u_u_c', 0, ifelse(treatment %in% c('u_u_n','u_u_p','u_m_c','b_u_c'), 1, ifelse(treatment %in% c('b_u_b','b_m_n','b_m_p','u_m_b'), 3, ifelse(treatment %in% c('b_u_n','b_u_p','u_u_b','u_m_n','u_m_p','b_m_c'), 2, 4)))))%>%
  mutate(resource_mani=ifelse(treatment %in% c('u_m_c','b_u_c','b_m_c'), 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='b_m_b', 'N*P*burn*mow_clip', ifelse(treatment=='b_m_c', 'burn*mow_clip', ifelse(treatment=='b_m_n', 'N*burn*mow_clip', ifelse(treatment=='b_m_p', 'P*burn*mow_clip', ifelse(treatment=='b_u_b', 'N*P*burn', ifelse(treatment=='b_u_c', 'burn', ifelse(treatment=='b_u_n', 'N*burn', ifelse(treatment=='b_u_p', 'P*burn', ifelse(treatment=='u_m_b', 'N*P*mow_clip', ifelse(treatment=='u_m_c', 'mow_clip', ifelse(treatment=='u_m_n', 'N*mow_clip', ifelse(treatment=='u_m_p', 'P*mow_clip', ifelse(treatment=='u_u_b', 'N*P', ifelse(treatment=='u_u_c', 'control', ifelse(treatment=='u_u_n', 'N', 'P'))))))))))))))))%>%
  unique()

irg<-read.delim("KNZ_IRG.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='i', 30, 0),
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='i', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='i', 'irr', 'control'))%>%
  unique()

gfp<-read.csv("KNZ_KNP_GFP.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment %in% c('Rainout_Grazed',"Rainout_Ungrazed"), -50, 0),
         temp=0, 
         mow_clip=ifelse(treatment %in% c("Rainout_Ungrazed","Open_Ungrazed"),0,1), 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='Open_Ungrazed', 0, ifelse(treatment=="Rainout_Grazed",2,1)))%>%
  mutate(resource_mani=ifelse(treatment=="Open_Grazed",0,1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='Open_Ungrazed', 'control', ifelse(treatment=='Open_Grazed', 'mow_clip', ifelse(treatment=='Rainout_Ungrazed', 'drought', 'drought*mow_clip'))))%>%
  unique()

pplots<-read.csv("KNZ_PPLOTS.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment %in% c('N1P0','N1P1','N1P2','N1P3'), 0, 10),
         p=ifelse(treatment %in% c('N1P1','N2P1'), 2.5, ifelse(treatment %in% c('N1P2','N2P2'), 5, ifelse(treatment %in% c('N1P3','N2P3'), 10, 0))),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N1P0', 0, ifelse(treatment %in% c('N1P1','N1P2','N1P3','N2P0'), 1, 2)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=ifelse(treatment %in% c('N1P0','N1P3','N2P0','N2P3'), 1, 0))%>%
  mutate(public=1)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment %in% c('N1P1','N1P2','N1P3'), 'P', ifelse(treatment=='N1P0', 'control', ifelse(treatment=='N2P0', 'N', 'N*P'))))%>%
  unique()

ramps<-read.csv("KNZ_Ramps.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=0,
         temp=ifelse(treatment %in% c('ambient_heated',"delayed_heated"), 1, 0),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=ifelse(treatment %in% c('delayed_control',"delayed_heated"),'increased precip vari', 'ambient'), 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='ambient_control', 0, ifelse(treatment=='delayed_heated', 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=='ambient_heated', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='ambient_control', 'control', ifelse(treatment=='ambient_heated', 'temp', ifelse(treatment=='delayed_control', 'precip_vari', 'precip_vari*temp'))))%>%
  unique()

rhps<-read.csv("KNZ_RHPs.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('N','stoneN'), 5, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         pulse=0,
         other_trt=ifelse(treatment %in% c('stone','stoneN', 'stoneC'), 'shallow soil', ifelse(treatment %in% c('C, stoneC'), 'carbon', 0)),
         trt_details=0,
         successional=1, 
         plant_mani=1,  
         plant_trt=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='control', 0, ifelse(treatment %in% c('stoneN', 'stoneC'), 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=='stone', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='control', 'control', 
                         ifelse(treatment=='N', 'N', 
                                ifelse(treatment=='stone', 'stone', 
                                       ifelse(treatment == 'stoneN', 'N*stone',
                                              ifelse(treatment == 'stoneC', 'C*stone', 'C'))))))%>%
  unique()

e2 <- read.csv("KUFS_E2.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type = 0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('N0S0H0', 'N0S0H1', 'N0S1H0', 'N0S1H1'), 0, 15),
         p=0,
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt= ifelse(treatment %in% c('N0S1H0', 'N0S1H1', 'N1S1H0', 'N1S1H1'), 'seed',
                                            ifelse(treatment %in% c('N0S0H1', 'N0S1H1', 'N1S0H1', 'N1S1H1'), 'mow_clip', 0)),
         trt_details=0,
         successional=1, 
         plant_mani=0,  
         plant_trt=0,
         pulse= ifelse(treatment %in% c('N0S1H0', 'N0S1H1', 'N1S1H0', 'N1S1H1'), 1, 0))%>%
  mutate(plot_mani= ifelse(treatment == 'N0S0H0', 0, 
                           ifelse(treatment == 'N1S1H1', 3,
                                  ifelse(treatment %in% c('N1S0H0', 'N0S1H0', 'N0S0H1'), 1, 2))))%>%
  mutate(resource_mani=ifelse(treatment %in% c("N0S1H0", "N0S1H1", "N0S1H1", "N0S0H1"), 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='N0S0H0', 'control',
                         ifelse(treatment =='N1S0H0', 'N', 
                                ifelse(treatment == 'N0S1H0', 'seed',
                                       ifelse(treatment == 'N0S0H1', 'mow_clip', 
                                              ifelse(treatment == 'N1S1H0', 'N*seed',
                                                     ifelse(treatment == 'N1S0H1', 'N*mow_clip',
                                                            ifelse(treatment == 'N0S1H1', 'seed*mow_clip', 'N*seed*mow_clip'))))))))%>%
  unique()

  
e6<-read.csv("KUFS_E6.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('N0P0S0','N0P8S0', "N0P8S1", "N0P0S1"), 0, 
                  ifelse(treatment %in% c('N4P0S0','N4P8S0', 'N4P0S1','N4P8S1'), 4, 
                         ifelse(treatment %in% c('N8P0S0','N8P8S0','N8P0S1','N8P8S1'), 8, 16))),
         p=ifelse(treatment %in% c('N0P0S0','N4P0S0','N8P0S0','N16P0S0',"N16P0S1","N8P0S1", "N4P0S1", "N0P0S1"), 0, 8),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt= ifelse(treatment %in% c("N0P0S1","N0P8S1","N4P0S1","N4P8S1","N8P0S1", "N8P8S1","N16P0S1", "N16P8S1"), "seed", 0), 
         trt_details=0,
         successional=1, 
         plant_mani=0,  
         plant_trt=0,
         pulse= ifelse(treatment %in% c("N0P0S1","N0P8S1","N4P0S1","N4P8S1","N8P0S1", "N8P8S1","N16P0S1", "N16P8S1"), 1, 0))%>%
  mutate(plot_mani=ifelse(treatment=='N0P0S0', 0, 
                          ifelse(treatment %in% c('N4P0S0','N8P0S0','N16P0S0','N0P8S0'), 1, 
                                 ifelse(treatment %in% c("N4P8S1", "N8P8S1", "N16P8S1"), 3, 2))))%>%
  mutate(resource_mani=ifelse(treatment == "N0P0S1", 0,1))%>%
  mutate(max_trt=ifelse(treatment %in% c('N0P0S0','N16P0S0','N0P8S0','N16P8S0'), 1, 0))%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='N0P0S0', 'control', ifelse(treatment=='N0P8S0', 'P', 
                                                                ifelse(treatment %in% c('N16P0S0','N8P0S0','N4P0S0'), 'N', 
                                                                       ifelse(treatment %in% c('N4P8S0', 'N8P8S0', 'N16P8S0'),'N*P',
                                                                              ifelse(treatment =='N0P0S1', 'seed',
                                                                                     ifelse(treatment %in% c('N4P0S1', 'N8P0S1', 'N16P0S1'), 'N*seed',
                                                                                            ifelse(treatment == 'N0P8S1', 'P*seed', "N*P*seed"))))))))%>%
  unique()

clip<-read.delim("LATNJA_CLIP.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('N','TN'), 5, 0),
         p=ifelse(treatment %in% c('N','TN'), 5, 0), 
         k=0, 
         CO2=0, 
         precip=0,
         temp=ifelse(treatment %in% c('T','TN'), 2, 0),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='CONTROL', 0, ifelse(treatment=='TN', 3, ifelse(treatment=='T',1,2))))%>%
  mutate(resource_mani=ifelse(treatment=='T', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='CONTROL', 'control', ifelse(treatment=='N', 'N*P', ifelse(treatment=='T', 'temp', 'N*P*temp'))))%>%
  unique()

pme<-read.csv("LEFT_PME.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0,
         p=0, 
         k=0, 
         CO2=0, 
         precip=ifelse(treatment=="winwet",50, ifelse(treatment=="winwet_sumwet", 100, 0)),
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=ifelse(treatment=='control', 'ambient precip', ifelse(treatment=='winwet', 'increase winter precip', ifelse(treatment=='winwet_sumdry','increase winter precip, decrease summer precip', ifelse(treatment=='winwet_sumwet', 'increase winter and summer precip', 'decrease winter precip increase summer precip')))), 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='control', 0, 1))%>% #we are considering this to be 1 manipulation, even when they manipulated precip in the winter and summer .
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='control', 'control', ifelse(treatment %in% c('windry_sumwet','winwet_sumdry'), 'precip_vari', 'irr')))%>%
  unique()

herbwood<-read.delim("LG_HerbWood.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='F'|treatment=='FW', 2.4, 0),
         p=ifelse(treatment=='F'|treatment=='FW', 0.66, 0),
         k=ifelse(treatment=='F'|treatment=='FW', 1.5, 0),
         CO2=0,
         precip=ifelse(treatment=='W'|treatment=='FW', 18, 0),
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment=='W', 1, ifelse(treatment=='F', 3, 4))))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='F', 'mult_nutrient', ifelse(treatment=='W', 'irr', ifelse(treatment=='FW', 'mult_nutrient*irr', 'control'))))%>%
  unique()

fireplots<-read.delim("MAERC_fireplots.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('snpu','snuu','unpu','unuu','wnpg','wnpu','wnug','wnuu'), 5, 0),
         p=ifelse(treatment %in% c('snpu','supu','unpu','uupu','wnpg','wnpu','wupg','wupu'), 2, 0),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0,
         burn=ifelse(treatment %in% c('uuuu','uupu','unpu','unuu'), 0, 1),
         herb_removal=ifelse(treatment %in% c('wnpg','wnug','wupg','wuug'), 1, 0), 
         management=0,
         other_trt=0,
         trt_details=ifelse(treatment %in% c('snpu','snuu','supu','suuu'), 'summer burn', ifelse(treatment %in% c('wnpg','wnpu','wnug','wnuu','wupg','wupu','wuug','wuuu'), 'winter burn', 0)),
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='wnpg', 4, ifelse(treatment %in% c('unpu','snpu','wnpu','wupg','wnug'), 3, ifelse(treatment %in% c('uupu','unuu','suuu','wuuu'), 1, ifelse(treatment=='uuuu', 0, 2)))))%>%
  mutate(resource_mani=ifelse(treatment %in% c('uuuu','wuuu','suuu','wuug'), 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='uuuu', 'control', ifelse(treatment=='wuug', 'burn*graze', ifelse(treatment=='unuu', 'N', ifelse(treatment=='wnug', 'N*burn*graze', ifelse(treatment=='unpu', 'N*P', ifelse(treatment=='wnpg', 'N*P*burn*graze', ifelse(treatment=='uupu', 'P', ifelse(treatment=='wupg', 'P*burn*graze', ifelse(treatment %in% c('suuu','wuuu'), 'burn', ifelse(treatment %in% c('snuu','wnuu'), 'N*burn', ifelse(treatment %in% c('supu','wupu'), 'P*burn', 'N*P*burn'))))))))))))%>%
  unique()

mwatfer<-read.csv("MNR_watfer.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment %in% c('F','FW'), 10, 0),
         p=ifelse(treatment %in% c('F','FW'), 10, 0),
         k=ifelse(treatment %in% c('F','FW'), 10, 0),
         CO2=0,
         precip=ifelse(treatment %in% c('W','FW'), 18, 0),
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment=='W', 1, ifelse(treatment=='F', 3, 4))))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='C', 'control', ifelse(treatment=='F', 'mult_nutrient', ifelse(treatment=='W', 'irr', 'mult_nutrient*irr'))))%>%
  unique()

wet<-read.delim("NANT_wet.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment %in% c('1N0P','1N1P'), 67.2, 0),
         p=ifelse(treatment %in% c('0N0P','1N0P'), 0, 33.6),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='0N0P', 0, ifelse(treatment=='1N1P', 2, 1)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='0N0P', 'control', ifelse(treatment=='0N1P', 'P', ifelse(treatment=='1N0P', 'N', 'N*P'))))%>%
  unique()

gb<-read.delim("NGBER_gb.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  filter(treatment!='AMBIENT')%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=ifelse(treatment=='AMBIENT', 'ambient rainfall', ifelse(treatment=='CURRENT', 'current pattern', ifelse(treatment=='SPRING', 'spring addition', 'winter addition'))), 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='CURRENT', 0, 1))%>%
  mutate(resource_mani=ifelse(treatment=='CURRENT', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment %in% c('CURRENT'), 'control', 'precip_vari'))%>%
  unique()

herbdiv<-read.csv("NIN_herbdiv.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('1NF','2NF','3NF','4NF','5NF'), 0, 12),
         p=ifelse(treatment %in% c('1NF','2NF','3NF','4NF','5NF'), 0, 3.3),
         k=ifelse(treatment %in% c('1NF','2NF','3NF','4NF','5NF'), 0, 8),
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=ifelse(treatment %in% c('1NF','1F'), 0, 1),
         management=1,
         other_trt=0,
         trt_details=ifelse(treatment %in% c('2NF','2F'), 'aboveground exclosure', ifelse(treatment %in% c('3NF','3F'), 'insecticide', ifelse(treatment %in% c('4NF','4F'), 'aboveground exclosure/insecticide', ifelse(treatment %in% c('5NF','5F'), 'above/below exclosure/insecticide', 0)))), 
         successional=0, 
         plant_mani=1,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='1NF', 0, ifelse(treatment %in% c('2NF','3NF'), 1, ifelse(treatment=='4NF', 2, ifelse(treatment %in% c('2F','3F'), 4, ifelse(treatment=='4F', 5, ifelse(treatment %in% c('1F','5NF'), 3, 6)))))))%>%
  mutate(resource_mani=ifelse(treatment %in% c('2NF','3NF','4NF','5NF'), 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='5NF', 'control', ifelse(treatment=='1F', 'mult_nutrient', ifelse(treatment  %in% c('2F','3F','4F','5F'), 'mult_nutrient*herb_removal', 'herb_removal'))))%>%
  unique()

ccd<-read.delim("NTG_CCD.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment %in% c('CH-','CL-','CN-','WH-','WL-','WN-'), -60, ifelse(treatment %in% c('CH+','CL+','CN+','WH+','WL+','WN+'), 60, 0)),
         temp=ifelse(site_code=='Alberta'&treatment %in% c('WH-','WHA','WH+','WL-','WLA','WL+','WN-','WNA','WN+'), 2.9,ifelse(site_code=='Manitoba'&treatment %in% c('WH-','WHA','WH+','WL-','WLA','WL+','WN-','WNA','WN+'), 1.4,ifelse(site_code=='Saskatchewan'&treatment %in% c('WH-','WHA','WH+','WL-','WLA','WL+','WN-','WNA','WN+'), 1.3,0))),
         mow_clip=ifelse(treatment %in% c('CN-','CNA','CN+','WN-','WN+','WNA'), 0, 1),
         burn=0, 
         herb_removal=0, 
         management=0,
         pulse=0,
         other_trt=0,
         trt_details=ifelse(treatment %in% c('CH-','CHA','WH-','CH+','WHA','WH+'), 'high intensity defoliation', ifelse(treatment %in% c('CL-','CLA','CL+','WL-','WLA','WL+'), 'low intensity defoliation', 0)),
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment %in% c('CHA','CLA','CN-','CN+','WNA'), 1,ifelse(treatment=='CNA', 0,ifelse(treatment %in% c('CH-','CH+','CL-','CL+','WHA','WLA','WNA','WN-','WN+'), 2, 3))))%>%
  mutate(resource_mani=ifelse(treatment %in% c('CHA','CLA','WHA','WLA','WNA'), 0, 1))%>%
  mutate(max_trt=ifelse(treatment %in% c('CH-','CHA','CH+','CN-','CNA','CN+','WH-','WHA','WH+','WN-','WNA','WN+'), 1, 0))%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='CNA', 'control', ifelse(treatment=='CN-', 'drought', ifelse(treatment %in% c('CH-','CL-'), 'drought*mow_clip', ifelse(treatment=='WN-', 'drought*temp', ifelse(treatment %in% c('WH-','WL-'),'drought*temp*mow_clip', ifelse(treatment=='CN+', 'irr', ifelse(treatment %in% c('CH+','CL+'), 'irr*mow_clip', ifelse(treatment=='WN+', 'irr*temp', ifelse(treatment %in% c('WH+','WL+'), 'irr*temp*mow_clip', ifelse(treatment %in% c('CHA','CLA'), 'mow_clip', ifelse(treatment=='WNA', 'temp', 'temp*mow_clip'))))))))))))%>%
  unique()

nutnet <- read.csv("NutNet.csv") %>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c("NPK+Fence", "N", "NK", "NP", "NPK"), 10, 0),
         p=ifelse(treatment %in% c("NPK+Fence", "PK", "P", "NP", "NPK"), 10, 0), 
         k=ifelse(treatment %in% c("NPK+Fence", "PK", "NK", "K", "NPK"), 10, 0), 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=ifelse(treatment %in% c("NPK+Fence", "Fence"), 1, 0),
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=="Control", 0, 
                          ifelse(treatment == "NPK+Fence", 4, 
                                 ifelse(treatment %in% c("N", "P", "K", "Fence"), 1, 
                                        ifelse(treatment == "NPK", 3, 2)))))%>%
  mutate(resource_mani=ifelse(treatment == "Fence", 0,1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment== "Control", "control", 
                         ifelse(treatment == "N", "N", 
                                ifelse(treatment == "P", "P", 
                                       ifelse(treatment == "K", "K", 
                                              ifelse(treatment == "NP", "N*P", 
                                                     ifelse(treatment == "Fence", "herb_removal", 
                                                            ifelse(treatment == "NPK+Fence", "mult_nutrient*herb_removal", "mult_nutrient"))))))))%>%
  unique()

nfert<-read.delim("NWT_246NFert.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='x', 0, ifelse(treatment=='low', 2, ifelse(treatment=='med', 4, 6))),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='x', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=ifelse(treatment %in% c('x','high'), 1, 0))%>%
  mutate(public=1)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='x', 'control', 'N'))%>%
  unique()

bowman<-read.delim("NWT_bowman.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment, community_type)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N'&calendar_year<=1991, 25, ifelse(treatment=='NP'&calendar_year<=1991, 25, ifelse(treatment %in% c('Control','P'), 0, 10))),
         p=ifelse(treatment=='P'&calendar_year<=1991, 25, ifelse(treatment=='NP'&calendar_year<=1991, 25, ifelse(treatment %in% c('Control','N'), 0, 10))),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0,
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='Control', 0, ifelse(treatment=='NP', 2, 1)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='Control', 'control', ifelse(treatment=='N', 'N', ifelse(treatment=='P', 'P', 'N*P'))))%>%
  unique()

snow<-read.csv("NWT_snow.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=1,
         n=ifelse(treatment=='XNX'&calendar_year<2011, 28, ifelse(treatment=='XNW'&calendar_year<2011, 28, ifelse(treatment=='PNX'&calendar_year<2011, 28, ifelse(treatment=='PNW'&calendar_year<2011, 28, ifelse(treatment=='XNX'&calendar_year>2010, 10, ifelse(treatment=='XNW'&calendar_year>=2011, 10, ifelse(treatment=='PNX'&calendar_year>=2011, 10, ifelse(treatment=='PNW'&calendar_year>=2011, 10, 0)))))))),
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment %in% c('XXX','XXW','XNX','XNW'), 0, 116),
         temp=ifelse(treatment %in% c('XXW','XNW','PXW','PNW'), 1, 0),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=1,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='XXX', 0, ifelse(treatment %in% c('XXW','XNX','PXX'), 1, ifelse(treatment %in% c('XNW','PXW','PNX'), 2, 3))))%>%
  mutate(resource_mani=ifelse(treatment=='XXW', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='XXX', 'control', ifelse(treatment=='XXW', 'temp', ifelse(treatment=='XNX', 'N', ifelse(treatment=='XNW', 'N*temp', ifelse(treatment=='PXX', 'irr', ifelse(treatment=='PXW', 'irr*temp', ifelse(treatment=='PNX', 'N*irr', 'N*irr*temp'))))))))%>%
  unique()

oface<-read.delim("ORNL_FACE.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=1, water=0, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=ifelse(treatment=='elevated', 170, 0),
         precip=0, 
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0, 
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='elevated', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='ambient', 'control', 'CO2'))%>%
  unique()

tide<-read.delim("PIE_Tide.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='Enriched', 37.5, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='Enriched', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='Reference', 'control', 'N'))%>%
  unique()

nut <- read.csv("Rengen_Nut.csv") %>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment %in% c("Control", "Ca"), 0, 10),
         p=ifelse(treatment %in% c("CaNP", "CaNP-KCl", "CaNP-K2SO4"), 3.5, 0), 
         k=ifelse(treatment %in% c("CaNP-KCl", "CaNP-K2SO4"), 13.3, 0), 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=ifelse(treatment == "Control", 0, ifelse(treatment == "Ca", "Ca addition", "Ca + Mg addition")), 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment== "Control", 0, ifelse(treatment == "Ca", 1, ifelse(treatment == "CaN", 3, ifelse(treatment == "CaNP", 4, 5)))))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=ifelse(treatment == "Ca", 0, 1))%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment== "Control", "control", ifelse(treatment == "Ca","Ca addition", "mult_nutrient")))%>%
  unique()

interaction<-read.delim("RIO_interaction.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment %in% c('N1W1','N1W2','N1W0'), 5, 0),
         p=0, 
         k=0,
         CO2=0,
         precip=ifelse(treatment %in% c('N0W0','N1W0','control'), 0, 27),
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0,
         trt_details=ifelse(treatment %in% c('N0W1','N1W1'), 'small precip pulse', ifelse(treatment %in% c('N0W2','N1W2'), 'large precip pulse', 0)),
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='control', 0, ifelse(treatment %in% c('N1W0','N0W1','N0W2'), 1, 2)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='control', 'control', ifelse(treatment=='N1W0', 'N', ifelse(treatment %in% c('N0W1','N0W2'), 'irr', 'N*irr'))))%>%
  unique()

lucero<-read.csv("SCL_Lucero.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N1', 20, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N1', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='N1', 'control', 'N'))%>%
  unique()

ter<-read.csv("SCL_TER.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('OF','CF'), 20, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=ifelse(treatment %in% c('CO','CF'), 1, 0), 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='OO', 0, ifelse(treatment=="CF",2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=="CO", 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='OO', 'control', ifelse(treatment=='OF', 'N', ifelse(treatment=='CO', 'mow_clip', 'N*mow_clip'))))%>%
  unique()

cxn<-read.csv("SERC_CXN.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=1, water=0, other_manipulation=0,
         n=ifelse(treatment %in% c('t2','t4'), 25, 0),
         p=0, 
         k=0, 
         CO2=ifelse(treatment %in% c('t3','t4'), 340,0),
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='t1', 0, ifelse(treatment=='t4',2,1)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='t1', 'control', ifelse(treatment=='t2', 'N', ifelse(treatment=='t3', 'CO2', 'N*CO2'))))%>%
  unique()

tmece<-read.csv("SERC_TMECE.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment, community_type)%>%
  mutate(nutrients=0, light=0, carbon=1, water=0, other_manipulation=0,
         n=0,
         p=0, 
         k=0, 
         CO2=ifelse(treatment=='E', 340,0), 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='A', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='E', 'CO2', 'control'))%>%
  unique()

snfert<-read.csv("SEV_NFert20.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='F', 10, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='F', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='F', 'N', 'control'))%>%
  unique()

wenndex<-read.csv("SEV_WENNDEx20.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=1, other_manipulation=1,
         n=ifelse(treatment %in% c('C','P','T','TP'), 0, 2),
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment %in% c('C','N','T','TN'), 0, 50),
         temp=ifelse(treatment %in% c('C','N','P','PN'), 0, 1),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment=='N'|treatment=='P'|treatment=='T', 1, ifelse(treatment=='TPN', 3, 2))))%>%
  mutate(resource_mani=ifelse(treatment=='T', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='C', 'control', ifelse(treatment=='N', 'N', ifelse(treatment=='P', 'irr', ifelse(treatment=='PN', 'N*irr', ifelse(treatment=='T', 'temp', ifelse(treatment=='TN', 'N*temp', ifelse(treatment=='TP', 'irr*temp', 'N*irr*temp'))))))))%>%
  unique()

grazeprecip<-read.csv("SFREC_GrazePrecip.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment, community_type)%>%
  mutate(nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0,
         p=0, 
         k=0, 
         CO2=0, 
         precip=ifelse(treatment=="C", 0, ifelse(treatment=='W',80,-80)), 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='C', 'control', ifelse(treatment=='W', 'irr', 'drought')))%>%
  unique()

sprecip <- read.csv("SGS_Precip.csv") %>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0,
         p=0, 
         k=0, 
         CO2=0, 
         precip=ifelse(treatment == "reduction", -60, ifelse(treatment == "add", 50,0)), 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=="control", 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=1)%>%
  mutate(public=1)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=="control", "control", ifelse(treatment == "add", "irr", "drought")))%>%
  unique()

nash <- read.csv("Sil_NASH.csv")%>% 
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(grepl("all.nutr", treatment, fixed = T),10,
                  ifelse(grepl("plus.N", treatment, fixed = T),10,
                         ifelse(grepl("min.pk", treatment, fixed = T),10,
                                ifelse(grepl("min.p", treatment, fixed = T),10,
                                       ifelse(grepl("min.k", treatment, fixed = T),10,
                                              ifelse(grepl("min.mg", treatment, fixed = T),10, 0)))))),
         p=ifelse(grepl("all.nutr", treatment, fixed = T),3.5,
                  ifelse(grepl("plus.p", treatment, fixed = T),3.5,
                         ifelse(grepl("plus.pk", treatment, fixed = T),3.5,
                                ifelse(grepl("min.n", treatment, fixed = T),3.5,
                                       ifelse(grepl("min.k", treatment, fixed = T),3.5,
                                              ifelse(grepl("min.mg", treatment, fixed = T),3.5,0)))))), 
         k=ifelse(grepl("all.nutr", treatment, fixed = T),22.5,
                  ifelse(grepl("plus.k", treatment, fixed = T),22.5,
                         ifelse(grepl("plus.pk", treatment, fixed = T),22.5,
                                ifelse(grepl("min.n", treatment, fixed = T),22.5,
                                       ifelse(grepl("min.p", treatment, fixed = T),22.5,
                                              ifelse(grepl("min.mg", treatment, fixed = T),22.5,0)))))), 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=ifelse(grepl("spray", treatment, fixed = T),1,
                             ifelse(grepl("pellets", treatment, fixed = T),1,
                                    ifelse(grepl("fenced", treatment, fixed = T),1,0))),
         management=0,
         other_trt=ifelse(grepl("limed.control", treatment, fixed = T),"limed",
                          ifelse(grepl("unlimed.min.herb", treatment, fixed = T),"herb-specific herbicide",
                                 ifelse(grepl("unlimed.min.grass", treatment, fixed = T), "grass-specific herbicide",
                                        ifelse(grepl("limed.min.herb", treatment, fixed = T),"lime + herb-specific herbicide",
                                               ifelse(grepl("limed.min.grass", treatment, fixed = T),"lime + grass-specific herbicide",0))))), 
         trt_details=0,
         successional=0, 
         plant_mani=0,  
         plant_trt=0,
         pulse= ifelse(grepl("min.herb", treatment, fixed = T),1,ifelse(grepl("min.grass", treatment, fixed = T), 1, 0)))%>%
  mutate(plot_mani=ifelse(treatment== "insects.molluscs.rabbits.unlimed.control.no.nutr", 0,
                          ifelse(treatment %in% c("insects.molluscs.rabbits.unlimed.control.plus.n", "insects.molluscs.rabbits.unlimed.control.plus.p",
                                                  "insects.molluscs.rabbits.unlimed.control.plus.k", "insects.molluscs.rabbits.unlimed.control.plus.mg",
                                                  "insects.molluscs.rabbits.unlimed.min.herb.no.nutr", "insects.molluscs.rabbits.unlimed.min.grass.no.nutr",
                                                  "insects.molluscs.rabbits.limed.control.no.nutr","insects.molluscs.fence.unlimed.control.no.nutr",
                                                  "insects.pellets.rabbits.unlimed.control.no.nutr","spray.molluscs.rabbits.unlimed.control.no.nutr"), 1,
                                 ifelse(treatment %in% c("spray.pellets.rabbits.unlimed.control.no.nutr","spray.molluscs.fenced.unlimed.control.no.nutr",
                                                         "spray.molluscs.rabbits.limed.control.no.nutr","spray.molluscs.rabbits.unlimed.min.herb.no.nutr",
                                                         "spray.molluscs.rabbits.unlimed.min.grass.no.nutr","spray.molluscs.rabbits.unlimed.control.plus.n",
                                                         "spray.molluscs.rabbits.unlimed.control.plus.p","spray.molluscs.rabbits.unlimed.control.plus.k",
                                                         "spray.molluscs.rabbits.unlimed.control.plus.mg","insects.pellets.fenced.unlimed.control.no.nutr",
                                                         "insects.pellets.rabbits.limed.control.no.nutr","insects.pellets.rabbits.unlimed.min.herb.no.nutr",
                                                         "insects.pellets.rabbits.unlimed.min.grass.no.nutr", "insects.pellets.rabbits.unlimed.control.plus.n",
                                                         "insects.pellets.rabbits.unlimed.control.plus.p", "insects.pellets.rabbits.unlimed.control.plus.k",
                                                         "insects.pellets.rabbits.unlimed.control.plus.mg","insects.molluscs.fenced.limed.control.no.nutr",
                                                         "insects.molluscs.fenced.unlimed.min.herb.no.nutr","insects.molluscs.fenced.unlimed.min.grass.no.nutr",
                                                         "insects.molluscs.fenced.unlimed.control.plus.n","insects.molluscs.fenced.unlimed.control.plus.p",
                                                         "insects.molluscs.fenced.unlimed.control.plus.k","insects.molluscs.fenced.unlimed.control.plus.mg",
                                                         "insects.molluscs.rabbits.limed.min.herb.no.nutr", "insects.molluscs.rabbits.limed.min.grass.no.nutr",
                                                         "insects.molluscs.rabbits.limed.control.plus.n","insects.molluscs.rabbits.limed.control.plus.p",
                                                         "insects.molluscs.rabbits.limed.control.plus.k","insects.molluscs.rabbits.limed.control.plus.mg",
                                                         "insects.molluscs.rabbits.unlimed.min.herb.plus.n","insects.molluscs.rabbits.unlimed.min.herb.plus.p",
                                                         "insects.molluscs.rabbits.unlimed.min.herb.plus.k","insects.molluscs.rabbits.unlimed.min.herb.plus.mg",
                                                         "insects.molluscs.rabbits.unlimed.min.grass.plus.n","insects.molluscs.rabbits.unlimed.min.grass.plus.p",
                                                         "insects.molluscs.rabbits.unlimed.min.grass.plus.k","insects.molluscs.rabbits.unlimed.min.grass.plus.mg",
                                                         "insects.molluscs.rabbits.unlimed.control.plus.pk","insects.molluscs.rabbits.unlimed.control.min.pk"),2,
                                        ifelse(treatment %in% c("spray.pellets.fenced.unlimed.control.no.nutr","spray.pellets.rabbits.limed.control.no.nutr",
                                                                "spray.pellets.rabbits.unlimed.min.grass.no.nutr","spray.pellets.rabbits.unlimed.min.herb.no.nutr",
                                                                "spray.pellets.rabbits.unlimed.control.plus.n","spray.pellets.rabbits.unlimed.control.plus.p",
                                                                "spray.pellets.rabbits.unlimed.control.plus.k","spray.pellets.rabbits.unlimed.control.plus.mg",
                                                                "spray.molluscs.fenced.limed.control.no.nutr","spray.molluscs.fenced.unlimed.min.grass.no.nutr",
                                                                "spray.molluscs.fenced.unlimed.min.herb.no.nutr","spray.molluscs.fenced.unlimed.control.plus.n",
                                                                "spray.molluscs.fenced.unlimed.control.plus.p","spray.molluscs.fenced.unlimed.control.plus.k",
                                                                "spray.molluscs.fenced.unlimed.control.plus.mg","spray.molluscs.rabbits.limed.min.grass.no.nutr",
                                                                "spray.molluscs.rabbits.limed.min.herb.no.nutr","spray.molluscs.rabbits.limed.control.plus.n",
                                                                "spray.molluscs.rabbits.limed.control.plus.p","spray.molluscs.rabbits.limed.control.plus.k",
                                                                "spray.molluscs.rabbits.limed.control.plus.mg","spray.molluscs.rabbits.unlimed.min.grass.plus.n",
                                                                "spray.molluscs.rabbits.unlimed.min.grass.plus.p","spray.molluscs.rabbits.unlimed.min.grass.plus.k",
                                                                "spray.molluscs.rabbits.unlimed.min.grass.plus.mg","spray.molluscs.rabbits.unlimed.min.herb.plus.n",
                                                                "spray.molluscs.rabbits.unlimed.min.herb.plus.p","spray.molluscs.rabbits.unlimed.min.herb.plus.k",
                                                                "spray.molluscs.rabbits.unlimed.min.herb.plus.mg","spray.molluscs.rabbits.unlimed.control.plus.pk",
                                                                "spray.molluscs.rabbits.unlimed.control.min.pk","insects.pellets.fenced.limed.control.no.nutr",
                                                                "insects.pellets.fenced.unlimed.min.grass.no.nutr", "insects.pellets.fenced.unlimed.min.herd.no.nutr",
                                                                "insects.pellets.fenced.unlimed.control.plus.n", "insects.pellets.fenced.unlimed.control.plus.p",
                                                                "insects.pellets.fenced.unlimed.control.plus.k","insects.pellets.fenced.unlimed.control.plus.mg",
                                                                "insects.pellets.rabbits.limed.min.grass.no.nutr", "insects.pellets.rabbits.limed.min.herb.no.nutr",
                                                                "insects.pellets.rabbits.limed.control.plus.n","insects.pellets.rabbits.limed.control.plus.p",
                                                                "insects.pellets.rabbits.limed.control.plus.k","insects.pellets.rabbits.limed.control.plus.mg",
                                                                "insects.pellets.rabbits.unlimed.min.grass.plus.n","insects.pellets.rabbits.unlimed.min.grass.plus.p",
                                                                "insects.pellets.rabbits.unlimed.min.grass.plus.k","insects.pellets.rabbits.unlimed.min.grass.plus.mg",
                                                                "insects.pellets.rabbits.unlimed.min.herb.plus.n","insects.pellets.rabbits.unlimed.min.herb.plus.p",
                                                                "insects.pellets.rabbits.unlimed.min.herb.plus.k","insects.pellets.rabbits.unlimed.min.herb.plus.mg",
                                                                "insects.pellets.rabbits.unlimed.control.plus.pk","insects.pellets.rabbits.unlimed.control.min.pk",
                                                                "insects.molluscs.fenced.limed.min.grass.no.nutr", "insects.molluscs.fenced.limed.min.herb.no.nutr",
                                                                "insects.molluscs.fenced.limed.control.plus.n","insects.molluscs.fenced.limed.control.plus.p",
                                                                "insects.molluscs.fenced.limed.control.plus.k","insects.molluscs.fenced.limed.control.plus.mg",
                                                                "insects.molluscs.fenced.unlimed.min.grass.plus.n","insects.molluscs.fenced.unlimed.min.grass.plus.p",
                                                                "insects.molluscs.fenced.unlimed.min.grass.plus.k","insects.molluscs.fenced.unlimed.min.grass.plus.mg",
                                                                "insects.molluscs.fenced.unlimed.min.herb.plus.n","insects.molluscs.fenced.unlimed.min.herb.plus.p",
                                                                "insects.molluscs.fenced.unlimed.min.herb.plus.k","insects.molluscs.fenced.unlimed.min.herb.plus.mg",
                                                                "insects.molluscs.fenced.unlimed.control.plus.pk","insects.molluscs.fenced.unlimed.control.min.pk",
                                                                "insects.molluscs.rabbits.limed.min.grass.plus.n","insects.molluscs.rabbits.limed.min.grass.plus.p",
                                                                "insects.molluscs.rabbits.limed.min.grass.plus.k","insects.molluscs.rabbits.limed.min.grass.plus.mg",
                                                                "insects.molluscs.rabbits.limed.min.herb.plus.n","insects.molluscs.rabbits.limed.min.herb.plus.p",
                                                                "insects.molluscs.rabbits.limed.min.herb.plus.k","insects.molluscs.rabbits.limed.min.herb.plus.mg",
                                                                "insects.molluscs.rabbits.limed.control.plus.pk","insects.molluscs.rabbits.unlimed.control.min.n",
                                                                "insects.molluscs.rabbits.unlimed.control.min.p","insects.molluscs.rabbits.unlimed.control.min.k",
                                                                "insects.molluscs.rabbits.unlimed.control.min.mg"),3, 
                                               ifelse(treatment %in% c("spray.pellets.fenced.limed.min.grass.no.nutr","spray.pellets.fenced.limed.min.herb.no.nutr",
                                                                       "spray.pellets.fenced.limed.control.plus.n","spray.pellets.fenced.limed.control.plus.p",
                                                                       "spray.pellets.fenced.limed.control.plus.k", "spray.pellets.fenced.limed.control.plus.mg",
                                                                       "spray.pellets.fenced.unlimed.min.grass.plus.n","spray.pellets.fenced.unlimed.min.grass.plus.p",
                                                                       "spray.pellets.fenced.unlimed.min.grass.plus.k","spray.pellets.fenced.unlimed.min.grass.plus.mg",
                                                                       "spray.pellets.fenced.unlimed.min.herb.plus.n","spray.pellets.fenced.unlimed.min.herb.plus.p",
                                                                       "spray.pellets.fenced.unlimed.min.herb.plus.k", "spray.pellets.fenced.unlimed.min.herb.plus.mg",
                                                                       "spray.pellets.rabbits.limed.min.grass.plus.n", "spray.pellets.rabbits.limed.min.grass.plus.p",
                                                                       "spray.pellets.rabbits.limed.min.grass.plus.k","spray.pellets.rabbits.limed.min.grass.plus.mg",
                                                                       "spray.pellets.rabbits.limed.min.herb.plus.n","spray.pellets.rabbits.limed.min.herb.plus.p",
                                                                       "spray.pellets.rabbits.limed.min.herb.plus.k","spray.pellets.rabbits.limed.min.herb.plus.mg",
                                                                       "insects.pellets.fenced.limed.min.grass.plus.n","insects.pellets.fenced.limed.min.grass.plus.p",
                                                                       "insects.pellets.fenced.limed.min.grass.plus.k","insects.pellets.fenced.limed.min.grass.plus.mg",
                                                                       "insects.pellets.fenced.limed.min.herb.plus.n", "insects.pellets.fenced.limed.min.herb.plus.p",
                                                                       "insects.pellets.fenced.limed.min.herb.plus.k", "insects.pellets.fenced.limed.min.herb.plus.mg",
                                                                       "spray.pellets.fenced.unlimed.control.plus.pk", "spray.pellets.fenced.unlimed.control.min.pk",
                                                                       "spray.molluscs.fenced.limed.control.plus.pk","spray.molluscs.fenced.limed.control.min.pk",
                                                                       "spray.molluscs.rabbits.limed.min.grass.plus.pk","spray.molluscs.rabbits.limed.min.grass.min.pk",
                                                                       "spray.molluscs.rabbits.limed.min.herb.plus.pk","spray.molluscs.rabbits.limed.min.herb.min.pk",
                                                                       "insects.pellets.fenced.limed.control.plus.pk","insects.pellets.fenced.limed.control.min.pk",
                                                                       "insects.pellets.fenced.unlimed.min.grass.plus.pk","insects.pellets.fenced.unlimed.min.grass.min.pk",
                                                                       "insects.pellets.fenced.unlimed.min.herb.plus.pk","insects.pellets.fenced.unlimed.min.herb.min.pk",
                                                                       "insects.pellets.rabbits.limed.min.grass.plus.pk","insects.pellets.rabbits.limed.min.grass.min.pk",
                                                                       "insects.pellets.rabbits.limed.min.herb.plus.pk","insects.pellets.rabbits.limed.min.herb.min.pk",
                                                                       "insects.molluscs.fenced.limed.min.grass.plus.pk","insects.molluscs.fenced.limed.min.grass.min.pk",
                                                                       "insects.molluscs.fenced.limed.min.herb.plus.pk","insects.molluscs.fenced.limed.min.herb.min.pk",
                                                                       "spray.pellets.rabbits.unlimed.control.min.n","spray.pellets.rabbits.unlimed.control.min.p",
                                                                       "spray.pellets.rabbits.unlimed.control.min.k", "spray.pellets.rabbits.unlimed.control.min.mg",
                                                                       "spray.molluscs.fenced.unlimed.control.min.n","spray.molluscs.fenced.unlimed.control.min.p",
                                                                       "spray.molluscs.fenced.unlimed.control.min.k","spray.molluscs.fenced.unlimed.control.min.mg",
                                                                       "spray.molluscs.rabbits.limed.control.min.n","spray.molluscs.rabbits.limed.control.min.p",
                                                                       "spray.molluscs.rabbits.limed.control.min.k","spray.molluscs.rabbits.limed.control.min.mg",
                                                                       "spray.molluscs.rabbits.unlimed.min.grass.min.n","spray.molluscs.rabbits.unlimed.min.grass.min.p",
                                                                       "spray.molluscs.rabbits.unlimed.min.grass.min.k","spray.molluscs.rabbits.unlimed.min.grass.min.mg",
                                                                       "spray.molluscs.rabbits.unlimed.min.herb.min.n","spray.molluscs.rabbits.unlimed.min.herb.min.p",
                                                                       "spray.molluscs.rabbits.unlimed.min.herb.min.k","spray.molluscs.rabbits.unlimed.min.herb.min.mg",
                                                                       "insects.pellets.fenced.unlimed.control.min.n","insects.pellets.rabbits.limed.control.min.n",
                                                                       "insects.pellets.rabbits.unlimed.min.grass.min.n","insects.pellets.rabbits.unlimed.min.herb.min.n",
                                                                       "insects.molluscs.fenced.limed.control.min.n","insects.molluscs.fenced.unlimed.min.grass.min.n",
                                                                       "insects.molluscs.fences.unlimed.min.herb.min.n","insects.molluscs.rabbits.limed.min.grass.min.n", 
                                                                       "insects.molluscs.rabbits.limed.min.herb.min.n","insects.pellets.fenced.unlimed.control.min.p",
                                                                       "insects.pellets.rabbits.limed.control.min.p","insects.pellets.rabbits.unlimed.min.grass.min.p",
                                                                       "insects.pellets.rabbits.unlimed.min.herb.min.p","insects.molluscs.fenced.limed.control.min.p",
                                                                       "insects.molluscs.fenced.unlimed.min.grass.min.p","insects.molluscs.fences.unlimed.min.herb.min.p",
                                                                       "insects.molluscs.rabbits.limed.min.grass.min.p", "insects.molluscs.rabbits.limed.min.herb.min.p",
                                                                       "insects.pellets.fenced.unlimed.control.min.k","insects.pellets.rabbits.limed.control.min.k",
                                                                       "insects.pellets.rabbits.unlimed.min.grass.min.k","insects.pellets.rabbits.unlimed.min.herb.min.k",
                                                                       "insects.molluscs.fenced.limed.control.min.k","insects.molluscs.fenced.unlimed.min.grass.min.k",
                                                                       "insects.molluscs.fences.unlimed.min.herb.min.k","insects.molluscs.rabbits.limed.min.grass.min.k", 
                                                                       "insects.molluscs.rabbits.limed.min.herb.min.k","insects.pellets.fenced.unlimed.control.min.mg",
                                                                       "insects.pellets.rabbits.limed.control.min.mg","insects.pellets.rabbits.unlimed.min.grass.min.mg",
                                                                       "insects.pellets.rabbits.unlimed.min.herb.min.mg","insects.molluscs.fenced.limed.control.min.mg",
                                                                       "insects.molluscs.fenced.unlimed.min.grass.min.mg","insects.molluscs.fences.unlimed.min.herb.min.mg",
                                                                       "insects.molluscs.rabbits.limed.min.grass.min.mg", "insects.molluscs.rabbits.limed.min.herb.min.mg", 
                                                                       "spray.molluscs.rabbits.unlimed.control.all.nutr", "insects.pellets.rabbits.unlimed.control.all.nutr",
                                                                       "insects.molluscs.fenced.unlimed.control.all.nutr", "insects.molluscs.rabbits.limed.control.all.nutr",
                                                                       "insects.molluscs.rabbits.unlimed.min.grass.all.nutr", "insects.molluscs.rabbits.unlimed.min.grass.all.nutr"),5,
                                                      ifelse(treatment %in% c("spray.pellets.fenced.limed.min.grass.plus.n","spray.pellets.fenced.limed.min.grass.plus.p",
                                                                              "spray.pellets.fenced.limed.min.grass.plus.k","spray.pellets.fenced.limed.min.grass.plus.mg",
                                                                              "spray.pellets.fenced.limed.min.herb.plus.n","spray.pellets.fenced.limed.min.herb.plus.p",
                                                                              "spray.pellets.fenced.limed.min.herb.plus.k","spray.pellets.fenced.limed.min.herb.plus.mg",
                                                                              "spray.pellets.fenced.limed.control.plus.pk",  "spray.pellets.fenced.unlimed.min.grass.plus.pk", 
                                                                              "spray.pellets.fenced.unlimed.min.herb.plus.pk","spray.pellets.rabbits.limed.min.grass.plus.pk",
                                                                              "spray.pellets.rabbits.limed.min.herb.plus.pk","spray.molluscs.fenced.limed.min.grass.plus.pk",
                                                                              "spray.molluses.fenced.limed.min.herb.plus.pk","insects.pellets.fenced.limed.min.grass.plus.pk", 
                                                                              "insects.pellets.fenced.limed.min.herb.plus.pk","spray.pellets.fenced.limed.control.min.pk", 
                                                                              "spray.pellets.fenced.unlimed.min.grass.min.pk", "spray.pellets.fenced.unlimed.min.herb.min.pk",
                                                                              "spray.pellets.rabbits.limed.min.grass.min.pk","spray.pellets.rabbits.limed.min.herb.min.pk",
                                                                              "spray.molluscs.fenced.limed.min.grass.min.pk","spray.molluses.fenced.limed.min.herb.min.pk", 
                                                                              "insects.pellets.fenced.limed.min.grass.min.pk", "insects.pellets.fenced.limed.min.herb.min.pk",
                                                                              "spray.pellets.fenced.unlimed.control.min.n","spray.pellets.rabbits.limed.control.min.n", 
                                                                              "spray.pellets.rabbits.unlimed.min.grass.min.n", "spray.pellets.rabbits.unlimed.min.herb.min.n",
                                                                              "spray.molluscs.fenced.limed.control.min.n","spray.molluscs.fenced.unlimed.min.grass.min.n", 
                                                                              "spray.molluscs.fenced.unlimed.min.herb.min.n","spray.molluscs.fenced.limed.control.min.n",
                                                                              "spray.molluscs.fenced.unlimed.min.grass.min.n","spray.molluscs.fenced.unlimed.min.herb.min.n",
                                                                              "spray.molluscs.rabbits.limed.min.grass.min.n","spray.molluscs.rabbits.limed.min.herb.min.n",
                                                                              "insects.pellets.fenced.limed.control.min.n","insects.pellets.fenced.unlimed.min.grass.min.n",
                                                                              "insects.pellets.fenced.unlimed.min.herb.min.n","insects.pellets.rabbits.limed.min.grass.min.n",
                                                                              "insects.pellets.rabbits.limed.min.herb.min.n","insects.molluscs.fenced.limed.min.grass.min.n",
                                                                              "insects.molluscs.fenced.limed.min.herb.min.n", "spray.pellets.fenced.unlimed.control.min.p",
                                                                              "spray.pellets.rabbits.limed.control.min.p", "spray.pellets.rabbits.unlimed.min.grass.min.p",
                                                                              "spray.pellets.rabbits.unlimed.min.herb.min.p","spray.molluscs.fenced.limed.control.min.p",
                                                                              "spray.molluscs.fenced.unlimed.min.grass.min.p", "spray.molluscs.fenced.unlimed.min.herb.min.p",
                                                                              "spray.molluscs.fenced.limed.control.min.p","spray.molluscs.fenced.unlimed.min.grass.min.p",
                                                                              "spray.molluscs.fenced.unlimed.min.herb.min.p","spray.molluscs.rabbits.limed.min.grass.min.p",
                                                                              "spray.molluscs.rabbits.limed.min.herb.min.p","insects.pellets.fenced.limed.control.min.p",
                                                                              "insects.pellets.fenced.unlimed.min.grass.min.p","insects.pellets.fenced.unlimed.min.herb.min.p",
                                                                              "insects.pellets.rabbits.limed.min.grass.min.p", "insects.pellets.rabbits.limed.min.herb.min.p",
                                                                              "insects.molluscs.fenced.limed.min.grass.min.p","insects.molluscs.fenced.limed.min.herb.min.p",
                                                                              "spray.pellets.fenced.unlimed.control.min.k","spray.pellets.rabbits.limed.control.min.k", 
                                                                              "spray.pellets.rabbits.unlimed.min.grass.min.k","spray.pellets.rabbits.unlimed.min.herb.min.k",
                                                                              "spray.molluscs.fenced.limed.control.min.k","spray.molluscs.fenced.unlimed.min.grass.min.k", 
                                                                              "spray.molluscs.fenced.unlimed.min.herb.min.k","spray.molluscs.fenced.limed.control.min.k",
                                                                              "spray.molluscs.fenced.unlimed.min.grass.min.k","spray.molluscs.fenced.unlimed.min.herb.min.k",
                                                                              "spray.molluscs.rabbits.limed.min.grass.min.k","spray.molluscs.rabbits.limed.min.herb.min.k",
                                                                              "insects.pellets.fenced.limed.control.min.k","insects.pellets.fenced.unlimed.min.grass.min.k",
                                                                              "insects.pellets.fenced.unlimed.min.herb.min.k","insects.pellets.rabbits.limed.min.grass.min.k",
                                                                              "insects.pellets.rabbits.limed.min.herb.min.k","insects.molluscs.fenced.limed.min.grass.min.k",
                                                                              "insects.molluscs.fenced.limed.min.herb.min.k","spray.pellets.fenced.unlimed.control.min.mg",
                                                                              "spray.pellets.rabbits.limed.control.min.mg","spray.pellets.rabbits.unlimed.min.grass.min.mg",
                                                                              "spray.pellets.rabbits.unlimed.min.herb.min.mg","spray.molluscs.fenced.limed.control.min.mg",
                                                                              "spray.molluscs.fenced.unlimed.min.grass.min.mg", "spray.molluscs.fenced.unlimed.min.herb.min.mg",
                                                                              "spray.molluscs.fenced.limed.control.min.mg","spray.molluscs.fenced.unlimed.min.grass.min.mg",
                                                                              "spray.molluscs.fenced.unlimed.min.herb.min.mg","spray.molluscs.rabbits.limed.min.grass.min.mg",
                                                                              "spray.molluscs.rabbits.limed.min.herb.min.mg","insects.pellets.fenced.limed.control.min.mg",
                                                                              "insects.pellets.fenced.unlimed.min.grass.min.mg","insects.pellets.fenced.unlimed.min.herb.min.mg",
                                                                              "insects.pellets.rabbits.limed.min.grass.min.mg","insects.pellets.rabbits.limed.min.herb.min.mg",
                                                                              "insects.molluscs.fenced.limed.min.grass.min.mg", "insects.molluscs.fenced.limed.min.herb.min.mg",
                                                                              "spray.pellets.rabbits.unlimed.control.all.nutr", "spray.molluscs.fenced.unlimed.control.all.nutr", 
                                                                              "spray.molluscs.rabbits.limed.control.all.nutr", "spray.molluscs.rabbits.unlimed.min.grass.all.nutr",
                                                                              "spray.molluscs.rabbits.unlimed.min.her.all.nutr", "insects.pellets.fenced.unlimed.control.all.nutr",
                                                                              "insects.pellets.rabbits.limed.control.all.nutr", "insects.pellets.rabbits.unlimed.min.grass.all.nutr", 
                                                                              "insects.pellets.rabbits.unlimed.min.herb.all.nutr", "insects.molluscs.fenced.limed.control.all.nutr", 
                                                                              "insects.molluscs.fenced.unlimed.min.grass.all.nutr", "insects.molluscs.fenced.unlimed.min.herb.all.nutr",
                                                                              "insects.molluscs.rabbits.limed.min.grass.all.nutr", "insects.molluscs.rabbits.limed.min.herb.all.nutr"),6,
                                                             ifelse(treatment %in% c("spray.pellets.fenced.limed.min.herb.plus.pk","spray.pellets.fenced.limed.min.herb.min.pk",
                                                                                     "spray.pellets.fenced.limed.min.grass.plus.pk","spray.pellets.fenced.limed.min.grass.min.pk",
                                                                                     "spray.pellets.fenced.limed.control.min.n","spray.pellets.fenced.unlimed.min.grass.min.n",
                                                                                     "spray.pellets.rabbits.limed.min.grass.min.n","spray.molluscs.fenced.limed.min.grass.min.n",
                                                                                     "insects.pellets.fenced.limed.min.grass.min.n","spray.pellets.fenced.unlimed.min.herb.min.n",
                                                                                     "spray.pellets.rabbits.limed.min.herb.min.n","spray.molluscs.fenced.limed.min.herb.min.n",
                                                                                     "insects.pellets.fenced.limed.min.herb.min.n","spray.pellets.fenced.limed.control.min.p",
                                                                                     "spray.pellets.fenced.unlimed.min.grass.min.p","spray.pellets.rabbits.limed.min.grass.min.p",
                                                                                     "spray.molluscs.fenced.limed.min.grass.min.p","insects.pellets.fenced.limed.min.grass.min.p",
                                                                                     "spray.pellets.fenced.unlimed.min.herb.min.p","spray.pellets.rabbits.limed.min.herb.min.p",
                                                                                     "spray.molluscs.fenced.limed.min.herb.min.p","insects.pellets.fenced.limed.min.herb.min.p",
                                                                                     "spray.pellets.fenced.limed.control.min.k","spray.pellets.fenced.unlimed.min.grass.min.k",
                                                                                     "spray.pellets.rabbits.limed.min.grass.min.k","spray.molluscs.fenced.limed.min.grass.min.k",
                                                                                     "insects.pellets.fenced.limed.min.grass.min.k","spray.pellets.fenced.unlimed.min.herb.min.k",
                                                                                     "spray.pellets.rabbits.limed.min.herb.min.k","spray.molluscs.fenced.limed.min.herb.min.k",
                                                                                     "insects.pellets.fenced.limed.min.herb.min.k","spray.pellets.fenced.limed.control.min.mg",
                                                                                     "spray.pellets.fenced.unlimed.min.grass.min.mg","spray.pellets.rabbits.limed.min.grass.min.mg",
                                                                                     "spray.molluscs.fenced.limed.min.grass.min.mg","insects.pellets.fenced.limed.min.grass.min.mg",
                                                                                     "spray.pellets.fenced.unlimed.min.herb.min.mg","spray.pellets.rabbits.limed.min.herb.min.mg",
                                                                                     "spray.molluscs.fenced.limed.min.herb.min.mg", "insects.pellets.fenced.limed.min.herb.min.mg", 
                                                                                     "spray.pellets.fenced.unlimed.control.all.nutr", "spray.pellets.rabbits.limed.control.all.nutr", 
                                                                                     "spray.pellets.rabbits.unlimed.min.grass.all.nutr", "spray.pellets.rabbits.unlimed.min.herb.all.nutr",
                                                                                     "spray.molluscs.fenced.limed.control.all.nutr", "spray.molluscs.fenced.unlimed.min.grass.all.nutr",
                                                                                     "spray.molluscs.fenced.unlimed.min.herb.all.nutr", "spray.molluscs.rabbits.limed.min.grass.all.nutr",
                                                                                     "spray.molluscs.rabbits.limed.min.herb.all.nutr", "insects.pellets.fenced.limed.control.all.nutr",
                                                                                     "insects.pellets.fenced.unlimed.min.grass.all.nutr", "insects.pellets.fenced.unlimed.min.herb.all.nutr",
                                                                                     "insects.molluscs.fenced.limed.min.grass.all.nutr", "insects.molluscs.fenced.limed.min.herb.all.nutr"),7,
                                                                    ifelse(treatment %in% c("spray.pellets.fenced.limed.min.herb.min.n","spray.pellets.fenced.limed.min.grass.min.n",
                                                                                            "spray.pellets.fenced.limed.min.herb.min.p","spray.pellets.fenced.limed.min.grass.min.p",
                                                                                            "spray.pellets.fenced.limed.min.herb.min.k","spray.pellets.fenced.limed.min.grass.min.k",
                                                                                            "spray.pellets.fenced.limed.min.herb.min.mg","spray.pellets.fenced.limed.min.grass.min.mg",
                                                                                            "spray.pellets.fenced.limed.control.all.nutr", "spray.pellets.fenced.unlimed.min.grass.all.nutr", 
                                                                                            "spray.pellets.fenced.unlimed.min.herb.all.nutr", "spray.pellets.rabbits.limed.min.grass.all.nutr",
                                                                                            "spray.pellets.rabbits.limed.min.herb.all.nutr","spary.molluscs.fenced.limed.min.grass.all.nutr",
                                                                                            "spray.molluscs.fenced.limed.min.herb.all.nutr", "insects.pellets.fenced.limed.min.grass.all.nutr",
                                                                                            "insects.pellets.fenced.limed.min.herb.all.nutr"),8,
                                                                           ifelse(treatment %in% c("spray.pellets.fenced.limed.min.grass.all.nutr", "spray.pellets.fenced.limed.min.herb.all.nutr"), 9, 4 ))))))))))%>%
  mutate(resource_mani=ifelse(grepl("no.nutr", treatment, fixed = T),0,1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment== "insects.molluscs.rabbits.unlimed.control.no.nutr", 'control', 
                         ifelse(treatment %in% c("spray.molluscs.rabbits.unlimed.control.no.nutr","spray.pellets.rabbits.unlimed.control.no.nutr","spray.molluscs.fenced.unlimed.control.no.nutr","spray.pellets.fenced.unlimed.control.no.nutr", "insects.pellets.rabbits.unlimed.control.no.nutr","insects.pellets.fenced.unlimed.control.no.nutr","insects.molluscs.fenced.unlimed.control.no.nutr"),"herb_removal",
                                ifelse(treatment == "insects.molluscs.rabbits.limed.control.no.nutr", "lime", 
                                       ifelse(treatment %in% c("insects.molluscs.rabbits.unlimed.min.herb.no.nutr", "insects.molluscs.rabbits.unlimed.min.grass.no.nutr"), "herbicide",
                                              ifelse(treatment  == "insects.molluscs.rabbits.unlimed.control.plus.n", "N", 
                                                     ifelse(treatment == "insects.molluscs.rabbits.unlimed.control.plus.p", "P", 
                                                            ifelse(treatment == "insects.molluscs.rabbits.unlimed.control.plus.k", "K", 
                                                                   ifelse(treatment == "insects.molluscs.rabbits.unlimed.control.plus.mg", "Mg", 
                                                                          ifelse(treatment %in% c("insects.molluscs.rabbits.unlimed.control.min.n", "insects.molluscs.rabbits.unlimed.control.min.p", "insects.molluscs.rabbits.unlimed.control.min.k", "insects.molluscs.rabbits.unlimed.control.min.mg", "insects.molluscs.rabbits.unlimed.control.all.nutr"), "mult_nutrient", 
                                                                                 ifelse(treatment %in% c("spray.pellets.fenced.limed.control.no.nutr", "spray.pellets.rabbits.limed.control.no.nutr", "spray.molluscs.rabbits.limed.control.no.nutr", "spray.molluscs.fenced.limed.no.nutr", 
                                                                                                         "insects.pellets.fenced.limed.control.no.nutr", "insects.pellets.rabbits.limed.control.no.nutr", "insects.molluscs.fenced.limed.control.no.nutr"),"herb_removal*lime",
                                                                                        ifelse(treatment %in% c("spray.pellets.fenced.unlimed.min.herb.no.nutr", "spray.pellets.rabbits.unlimed.min.herb.no.nutr", "spray.molluscs.rabbits.unlimed.min.herb.no.nutr", "spray.molluscs.unlimed.min.herb.no.nutr", 
                                                                                                                "insects.pellets.fenced.unlimed.min.herb.no.nutr", "insects.pellets.rabbits.unlimed.min.herb.no.nutr", "insects.molluscs.fenced.unlimed.min.herb.no.nutr",
                                                                                                                "spray.pellets.fenced.unlimed.min.grass.no.nutr", "spray.pellets.rabbits.unlimed.min.grass.no.nutr", "spray.molluscs.rabbits.unlimed.min.grass.no.nutr", "spray.molluscs.unlimed.min.grass.no.nutr", 
                                                                                                                "insects.pellets.fenced.unlimed.min.grass.no.nutr", "insects.pellets.rabbits.unlimed.min.grass.no.nutr", "insects.molluscs.fenced.unlimed.min.grass.no.nutr"), "herb_removal*herbicide",
                                                                                               ifelse(treatment %in% c("spray.pellets.fenced.unlimed.control.plus.n", "spray.pellets.rabbits.unlimed.control.plus.n", "spray.molluscs.rabbits.unlimed.control.plus.n", "spray.molluscs.unlimed.control.plus.n", 
                                                                                                      "insects.pellets.fenced.unlimed.control.plus.n", "insects.pellets.rabbits.unlimed.control.plus.n", "insects.molluscs.fenced.unlimed.control.plus.n"),"N*herb_removal",
                                                                                               ifelse(treatment %in% c("spray.pellets.fenced.unlimed.control.plus.p", "spray.pellets.rabbits.unlimed.control.plus.p", "spray.molluscs.rabbits.unlimed.control.plus.p", "spray.molluscs.unlimed.control.plus.p", 
                                                                                                                       "insects.pellets.fenced.unlimed.control.plus.p", "insects.pellets.rabbits.unlimed.control.plus.p", "insects.molluscs.fenced.unlimed.control.plus.p"), "P*herb_removal",
                                                                                                      ifelse(treatment %in% c("spray.pellets.fenced.unlimed.control.plus.k", "spray.pellets.rabbits.unlimed.control.plus.k", "spray.molluscs.rabbits.unlimed.control.plus.k", "spray.molluscs.unlimed.control.plus.k", 
                                                                                                                              "insects.pellets.fenced.unlimed.control.plus.k", "insects.pellets.rabbits.unlimed.control.plus.k", "insects.molluscs.fenced.unlimed.control.plus.k"), "K*herb_removal",
                                                                                                             ifelse(treatment %in% c("spray.pellets.fenced.unlimed.control.plus.mg", "spray.pellets.rabbits.unlimed.control.plus.mg", "spray.molluscs.rabbits.unlimed.control.plus.mg", "spray.molluscs.unlimed.control.plus.mg", 
                                                                                                      "insects.pellets.fenced.unlimed.control.plus.mg", "insects.pellets.rabbits.unlimed.control.plus.mg", "insects.molluscs.fenced.unlimed.control.plus.mg"), "Mg*herb_removal",
                                                                                                      ifelse(treatment %in% c("spray.pellets.fenced.unlimed.control.all.nutr", "spray.pellets.rabbits.unlimed.control.plus.pk", "spray.molluscs.rabbits.unlimed.control.min.pk", "spray.molluscs.unlimed.control.min.n", 
                                                                                                      "insects.pellets.fenced.unlimed.control.min.p", "insects.pellets.rabbits.unlimed.control.min.k", "insects.molluscs.fenced.unlimed.control.min.mg"), "mult_nutrient*herb_removal",
                                                                                                      ifelse(treatment %in% c("insects.molluscs.rabbits.limed.min.grass.no.nutr", "insects.molluscs.rabbits.limed.min.herb.no.nutr"), "lime*herbicide",
                                                                                                      ifelse(treatment %in% c("insects.molluscs.rabbits.limed.control.all.nutr", "insects.molluscs.rabbits.limed.control.plus.pk", "insects.molluscs.rabbits.limed.control.min.pk", "insects.molluscs.rabbits.limed.control.min.n", "insects.molluscs.rabbits.limed.control.min.p", "insects.molluscs.rabbits.limed.control.min.k", "insects.molluscs.rabbits.limed.control.min.mg"), "mult_nutrient*lime",
                                                                                                             ifelse(treatment == "insects.molluscs.rabbits.limed.control.plus.n", "N*lime",
                                                                                                                    ifelse(treatment == "insects.molluscs.rabbits.limed.control.plus.p", "P*lime",
                                                                                                                           ifelse(treatment == "insects.molluscs.rabbits.limed.control.plus.K", "K*lime",
                                                                                                                                  ifelse(treatment == "insects.molluscs.rabbits.limed.control.plus.mg", "Mg*lime",
                                                                                                                                         ifelse(treatment %in% c("insects.molluscs.rabbits.unlimed.min.grass.all.nutr", "insects.molluscs.rabbits.unlimed.min.grass.plus.pk", "insects.molluscs.rabbits.unlimed.min.grass.min.pk", "insects.molluscs.rabbits.unlimed.min.grass.min.n", "insects.molluscs.rabbits.unlimed.min.grass.min.p", "insects.molluscs.rabbits.unlimed.min.grass.min.k", "insects.molluscs.rabbits.unlimed.min.grass.min.mg",
                                                                                                                                                                 "insects.molluscs.rabbits.unlimed.min.herb.all.nutr", "insects.molluscs.rabbits.unlimed.min.herb.plus.pk", "insects.molluscs.rabbits.unlimed.min.herb.min.pk", "insects.molluscs.rabbits.unlimed.min.herb.min.n", "insects.molluscs.rabbits.unlimed.min.herb.min.p", "insects.molluscs.rabbits.unlimed.min.herb.min.k", "insects.molluscs.rabbits.unlimed.min.herb.min.mg"),"mult_nutrient*herbicide",
                                                                                                                                                ifelse(treatment %in% c("insects.molluscs.rabbits.unlimed.min.grass.plus.n", "insects.molluscs.rabbits.unlimed.min.herb.plus.n"), "N*herbicide",
                                                                                                                                                       ifelse(treatment %in% c("insects.molluscs.rabbits.unlimed.min.grass.plus.p", "insects.molluscs.rabbits.unlimed.min.herb.plus.p"), "P*herbicide",
                                                                                                                                                              ifelse(treatment %in% c("insects.molluscs.rabbits.unlimed.min.grass.plus.k", "insects.molluscs.rabbits.unlimed.min.herb.plus.k"), "K*herbicide",
                                                                                                                                                                     ifelse(treatment %in% c("insects.molluscs.rabbits.unlimed.min.grass.plus.mg", "insects.molluscs.rabbits.unlimed.min.herb.plus.mg"), "Mg*herbicide",
                                                                                                                                                                            ifelse(treatment %in% c("spray.pellets.fenced.limed.min.grass.no.nutr","spray.pellets.rabbits.limed.min.grass.no.nutr",
                                                                                                                                                                                                    "spray.molluscs.fenced.limed.min.grass.no.nutr", "spray.molluscs.fenced.limed.min.grass.no.nutr",
                                                                                                                                                                                                    "insects.pellets.fenced.limed.min.grass.no.nutr","insects.pellets.rabbits.limed.min.grass.no.nutr",
                                                                                                                                                                                                    "insects.molluscs.rabbits.limed.min.grass.no.nutr",
                                                                                                                                                                                                    "spray.pellets.fenced.limed.min.herb.no.nutr","spray.pellets.rabbits.limed.min.herb.no.nutr",
                                                                                                                                                                                                    "spray.molluscs.fenced.limed.min.herb.no.nutr", "spray.molluscs.fenced.limed.min.herb.no.nutr",
                                                                                                                                                                                                    "insects.pellets.fenced.limed.min.herb.no.nutr","insects.pellets.rabbits.limed.min.herb.no.nutr",
                                                                                                                                                                                                    "insects.molluscs.rabbits.limed.min.herb.no.nutr"), "herb_removal*lime*herbicide",
                                                                                                                                                                                   ifelse(treatment %in% c("spray.pellets.fenced.limed.control.plus.n","spray.pellets.rabbits.limed.control.plus.n",
                                                                                                                                                                                          "spray.molluscs.fenced.limed.control.plus.n", "spray.molluscs.fenced.limed.control.plus.n",
                                                                                                                                                                                          "insects.pellets.fenced.limed.control.plus.n","insects.pellets.rabbits.limed.control.plus.n",
                                                                                                                                                                                          "insects.molluscs.rabbits.limed.control.plus.n"), "N*herb_removal*lime",
                                                                                                                                                                                   ifelse(treatment %in% c("spray.pellets.fenced.limed.control.plus.p","spray.pellets.rabbits.limed.control.plus.p",
                                                                                                                                                                                                           "spray.molluscs.fenced.limed.control.plus.p", "spray.molluscs.fenced.limed.control.plus.p",
                                                                                                                                                                                                           "insects.pellets.fenced.limed.control.plus.p","insects.pellets.rabbits.limed.control.plus.p",
                                                                                                                                                                                                           "insects.molluscs.rabbits.limed.control.plus.p"), "P*herb_removal*lime",
                                                                                                                                                                                          ifelse(treatment %in% c("spray.pellets.fenced.limed.control.plus.k","spray.pellets.rabbits.limed.control.plus.k",
                                                                                                                                                                                                                  "spray.molluscs.fenced.limed.control.plus.k", "spray.molluscs.fenced.limed.control.plus.k",
                                                                                                                                                                                                                  "insects.pellets.fenced.limed.control.plus.k","insects.pellets.rabbits.limed.control.plus.k",
                                                                                                                                                                                                                  "insects.molluscs.rabbits.limed.control.plus.k"), "K*herb_removal*lime",
                                                                                                                                                                                                 ifelse(treatment %in% c("spray.pellets.fenced.limed.control.plus.mg","spray.pellets.rabbits.limed.control.plus.mg",
                                                                                                                                                                                                                         "spray.molluscs.fenced.limed.control.plus.mg", "spray.molluscs.fenced.limed.control.plus.mg",
                                                                                                                                                                                                                         "insects.pellets.fenced.limed.control.plus.mg","insects.pellets.rabbits.limed.control.plus.mg",
                                                                                                                                                                                                                         "insects.molluscs.rabbits.limed.control.plus.mg"), "Mg*herb_removal*lime",
                                                                                                                                                                                                        ifelse(treatment %in% c("spray.pellets.fenced.limed.control.min.n","spray.pellets.rabbits.limed.control.min.n",
                                                                                                                                                                                                                                "spray.molluscs.fenced.limed.control.min.n", "spray.molluscs.fenced.limed.control.min.n",
                                                                                                                                                                                                                                "insects.pellets.fenced.limed.control.min.n","insects.pellets.rabbits.limed.control.min.n",
                                                                                                                                                                                                                                "insects.molluscs.rabbits.limed.control.min.n",
                                                                                                                                                                                                                                "spray.pellets.fenced.limed.control.plus.pk","spray.pellets.rabbits.limed.control.plus.pk",
                                                                                                                                                                                                                                "spray.molluscs.fenced.limed.control.plus.pk", "spray.molluscs.fenced.limed.control.plus.pk",
                                                                                                                                                                                                                                "insects.pellets.fenced.limed.control.plus.pk","insects.pellets.rabbits.limed.control.plus.pk",
                                                                                                                                                                                                                                "insects.molluscs.rabbits.limed.control.plus.pk",
                                                                                                                                                                                                                                "spray.pellets.fenced.limed.control.min.pk","spray.pellets.rabbits.limed.control.min.pk",
                                                                                                                                                                                                                                "spray.molluscs.fenced.limed.control.min.pk", "spray.molluscs.fenced.limed.control.min.pk",
                                                                                                                                                                                                                                "insects.pellets.fenced.limed.control.min.pk","insects.pellets.rabbits.limed.control.min.pk",
                                                                                                                                                                                                                                "insects.molluscs.rabbits.limed.control.min.pk",
                                                                                                                                                                                                                                "spray.pellets.fenced.limed.control.min.p","spray.pellets.rabbits.limed.control.min.p",
                                                                                                                                                                                                                                "spray.molluscs.fenced.limed.control.min.p", "spray.molluscs.fenced.limed.control.min.p",
                                                                                                                                                                                                                                "insects.pellets.fenced.limed.control.min.p","insects.pellets.rabbits.limed.control.min.p",
                                                                                                                                                                                                                                "insects.molluscs.rabbits.limed.control.min.p",
                                                                                                                                                                                                                                "spray.pellets.fenced.limed.control.min.k","spray.pellets.rabbits.limed.control.min.k",
                                                                                                                                                                                                                                "spray.molluscs.fenced.limed.control.min.k", "spray.molluscs.fenced.limed.control.min.k",
                                                                                                                                                                                                                                "insects.pellets.fenced.limed.control.min.k","insects.pellets.rabbits.limed.control.min.k",
                                                                                                                                                                                                                                "insects.molluscs.rabbits.limed.control.min.k",
                                                                                                                                                                                                                                "spray.pellets.fenced.limed.control.min.mg","spray.pellets.rabbits.limed.control.min.mg",
                                                                                                                                                                                                                                "spray.molluscs.fenced.limed.control.min.mg", "spray.molluscs.fenced.limed.control.min.mg",
                                                                                                                                                                                                                                "insects.pellets.fenced.limed.control.min.mg","insects.pellets.rabbits.limed.control.min.mg",
                                                                                                                                                                                                                                "insects.molluscs.rabbits.limed.control.min.mg",
                                                                                                                                                                                                                                "spray.pellets.fenced.limed.control.all.nutr","spray.pellets.rabbits.limed.control.all.nutr",
                                                                                                                                                                                                                                "spray.molluscs.fenced.limed.control.all.nutr", "spray.molluscs.fenced.limed.control.all.nutr",
                                                                                                                                                                                                                                "insects.pellets.fenced.limed.control.all.nutr","insects.pellets.rabbits.limed.control.all.nutr",
                                                                                                                                                                                                                                "insects.molluscs.rabbits.limed.control.all.nutr"), "mult_nutrient*herb_removal*lime",
                                                                                                                                                                                                               ifelse(treatment %in% c("spray.pellets.fenced.unlimed.min.grass.plus.n",
                                                                                                                                                                                                                      "spray.pellets.rabbits.unlimed.min.grass.plus.n",
                                                                                                                                                                                                                      "spray.molluscs.fenced.unlimed.min.grass.plus.n",
                                                                                                                                                                                                                      "spray.molluscs.fenced.unlimed.min.grass.plus.n",
                                                                                                                                                                                                                      "insects.pellets.fenced.unlimed.min.grass.plus.n",
                                                                                                                                                                                                                      "insects.pellets.rabbits.unlimed.min.grass.plus.n",
                                                                                                                                                                                                                      "insects.molluscs.rabbits.unlimed.min.grass.plus.n",
                                                                                                                                                                                                                      "spray.pellets.fenced.unlimed.min.herb.plus.n",
                                                                                                                                                                                                                      "spray.pellets.rabbits.unlimed.min.herb.plus.n",
                                                                                                                                                                                                                      "spray.molluscs.fenced.unlimed.min.herb.plus.n",
                                                                                                                                                                                                                      "spray.molluscs.fenced.unlimed.min.herb.plus.n",
                                                                                                                                                                                                                      "insects.pellets.fenced.unlimed.min.herb.plus.n",
                                                                                                                                                                                                                      "insects.pellets.rabbits.unlimed.min.herb.plus.n",
                                                                                                                                                                                                                      "insects.molluscs.rabbits.unlimed.min.herb.plus.n"), "N*herb_removal*herbicide",
                                                                                                                                                                                                               ifelse(treatment %in% c("spray.pellets.fenced.unlimed.min.grass.plus.p",
                                                                                                                                                                                                                                       "spray.pellets.rabbits.unlimed.min.grass.plus.p",
                                                                                                                                                                                                                                       "spray.molluscs.fenced.unlimed.min.grass.plus.p",
                                                                                                                                                                                                                                       "spray.molluscs.fenced.unlimed.min.grass.plus.p",
                                                                                                                                                                                                                                       "insects.pellets.fenced.unlimed.min.grass.plus.p",
                                                                                                                                                                                                                                       "insects.pellets.rabbits.unlimed.min.grass.plus.p",
                                                                                                                                                                                                                                       "insects.molluscs.rabbits.unlimed.min.grass.plus.p",
                                                                                                                                                                                                                                       "spray.pellets.fenced.unlimed.min.herb.plus.p",
                                                                                                                                                                                                                                       "spray.pellets.rabbits.unlimed.min.herb.plus.p",
                                                                                                                                                                                                                                       "spray.molluscs.fenced.unlimed.min.herb.plus.p",
                                                                                                                                                                                                                                       "spray.molluscs.fenced.unlimed.min.herb.plus.p",
                                                                                                                                                                                                                                       "insects.pellets.fenced.unlimed.min.herb.plus.p",
                                                                                                                                                                                                                                       "insects.pellets.rabbits.unlimed.min.herb.plus.p",
                                                                                                                                                                                                                                       "insects.molluscs.rabbits.unlimed.min.herb.plus.p"), "P*herb_removal*herbicide",
                                                                                                                                                                                                                      ifelse(treatment %in% c("spray.pellets.fenced.unlimed.min.grass.plus.k",
                                                                                                                                                                                                                                              "spray.pellets.rabbits.unlimed.min.grass.plus.k",
                                                                                                                                                                                                                                              "spray.molluscs.fenced.unlimed.min.grass.plus.k",
                                                                                                                                                                                                                                              "spray.molluscs.fenced.unlimed.min.grass.plus.k",
                                                                                                                                                                                                                                              "insects.pellets.fenced.unlimed.min.grass.plus.k",
                                                                                                                                                                                                                                              "insects.pellets.rabbits.unlimed.min.grass.plus.k",
                                                                                                                                                                                                                                              "insects.molluscs.rabbits.unlimed.min.grass.plus.k",
                                                                                                                                                                                                                                              "spray.pellets.fenced.unlimed.min.herb.plus.k",
                                                                                                                                                                                                                                              "spray.pellets.rabbits.unlimed.min.herb.plus.k",
                                                                                                                                                                                                                                              "spray.molluscs.fenced.unlimed.min.herb.plus.k",
                                                                                                                                                                                                                                              "spray.molluscs.fenced.unlimed.min.herb.plus.k",
                                                                                                                                                                                                                                              "insects.pellets.fenced.unlimed.min.herb.plus.k",
                                                                                                                                                                                                                                              "insects.pellets.rabbits.unlimed.min.herb.plus.k",
                                                                                                                                                                                                                                              "insects.molluscs.rabbits.unlimed.min.herb.plus.k"), "K*herb_removal*herbicide",
                                                                                                                                                                                                                             ifelse(treatment %in% c("spray.pellets.fenced.unlimed.min.grass.plus.mg",
                                                                                                                                                                                                                                                     "spray.pellets.rabbits.unlimed.min.grass.plus.mg",
                                                                                                                                                                                                                                                     "spray.molluscs.fenced.unlimed.min.grass.plus.mg",
                                                                                                                                                                                                                                                     "spray.molluscs.fenced.unlimed.min.grass.plus.mg",
                                                                                                                                                                                                                                                     "insects.pellets.fenced.unlimed.min.grass.plus.mg",
                                                                                                                                                                                                                                                     "insects.pellets.rabbits.unlimed.min.grass.plus.mg",
                                                                                                                                                                                                                                                     "insects.molluscs.rabbits.unlimed.min.grass.plus.mg",
                                                                                                                                                                                                                                                     "spray.pellets.fenced.unlimed.min.herb.plus.mg",
                                                                                                                                                                                                                                                     "spray.pellets.rabbits.unlimed.min.herb.plus.mg",
                                                                                                                                                                                                                                                     "spray.molluscs.fenced.unlimed.min.herb.plus.mg",
                                                                                                                                                                                                                                                     "spray.molluscs.fenced.unlimed.min.herb.plus.mg",
                                                                                                                                                                                                                                                     "insects.pellets.fenced.unlimed.min.herb.plus.mg",
                                                                                                                                                                                                                                                     "insects.pellets.rabbits.unlimed.min.herb.plus.mg",
                                                                                                                                                                                                                                                     "insects.molluscs.rabbits.unlimed.min.herb.plus.mg"), "Mg*herb_removal*herbicide",
                                                                                                                                                                                                                                    ifelse(treatment %in% c("spray.pellets.fenced.unlimed.min.grass.plus.pk",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.grass.plus.pk",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.plus.pk",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.plus.pk",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.grass.plus.pk",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.grass.plus.pk",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.grass.plus.pk",
                                                                                                                                                                                                                                                            "spray.pellets.fenced.unlimed.min.herb.plus.pk",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.herb.plus.pk",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.plus.pk",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.plus.pk",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.herb.plus.pk",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.herb.plus.pk",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.herb.plus.pk",
                                                                                                                                                                                                                                                            "spray.pellets.fenced.unlimed.min.grass.min.n",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.grass.min.n",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.min.n",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.min.n",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.grass.min.n",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.grass.min.n",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.grass.min.n",
                                                                                                                                                                                                                                                            "spray.pellets.fenced.unlimed.min.herb.min.n",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.herb.min.n",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.min.n",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.min.n",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.herb.min.n",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.herb.min.n",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.herb.min.n",
                                                                                                                                                                                                                                                            "spray.pellets.fenced.unlimed.min.grass.min.p",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.grass.min.p",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.min.p",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.min.p",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.grass.min.p",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.grass.min.p",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.grass.min.p",
                                                                                                                                                                                                                                                            "spray.pellets.fenced.unlimed.min.herb.min.p",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.herb.min.p",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.min.p",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.min.p",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.herb.min.p",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.herb.min.p",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.herb.min.p",
                                                                                                                                                                                                                                                            "spray.pellets.fenced.unlimed.min.grass.min.k",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.grass.min.k",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.min.k",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.min.k",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.grass.min.k",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.grass.min.k",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.grass.min.k",
                                                                                                                                                                                                                                                            "spray.pellets.fenced.unlimed.min.herb.min.k",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.herb.min.k",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.min.k",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.min.k",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.herb.min.k",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.herb.min.k",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.herb.min.k",
                                                                                                                                                                                                                                                            "spray.pellets.fenced.unlimed.min.grass.min.mg",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.grass.min.mg",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.min.mg",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.min.mg",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.grass.min.mg",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.grass.min.mg",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.grass.min.mg",
                                                                                                                                                                                                                                                            "spray.pellets.fenced.unlimed.min.herb.min.mg",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.herb.min.mg",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.min.mg",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.min.mg",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.herb.min.mg",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.herb.min.mg",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.herb.min.mg",
                                                                                                                                                                                                                                                            "spray.pellets.fenced.unlimed.min.grass.min.pk",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.grass.min.pk",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.min.pk",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.min.pk",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.grass.min.pk",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.grass.min.pk",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.grass.min.pk",
                                                                                                                                                                                                                                                            "spray.pellets.fenced.unlimed.min.herb.min.pk",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.herb.min.pk",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.min.pk",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.min.pk",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.herb.min.pk",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.herb.min.pk",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.herb.min.pk",
                                                                                                                                                                                                                                                            "spray.pellets.fenced.unlimed.min.grass.all.nutr",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.grass.all.nutr",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.all.nutr",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.grass.all.nutr",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.grass.all.nutr",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.grass.all.nutr",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.grass.all.nutr",
                                                                                                                                                                                                                                                            "spray.pellets.fenced.unlimed.min.herb.all.nutr",
                                                                                                                                                                                                                                                            "spray.pellets.rabbits.unlimed.min.herb.all.nutr",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.all.nutr",
                                                                                                                                                                                                                                                            "spray.molluscs.fenced.unlimed.min.herb.all.nutr",
                                                                                                                                                                                                                                                            "insects.pellets.fenced.unlimed.min.herb.all.nutr",
                                                                                                                                                                                                                                                            "insects.pellets.rabbits.unlimed.min.herb.all.nutr",
                                                                                                                                                                                                                                                            "insects.molluscs.rabbits.unlimed.min.herb.all.nutr"), "mult_nutrient*herb_removal*herbicide",
                                                                                                                                                                                                                                           ifelse(treatment %in% c("insects.molluscs.rabbits.limed.min.grass.plus.n", "insects.molluscs.rabbits.limed.min.herb.plus.n"), "N*lime*herbicide",
                                                                                                                                                                                                                                           ifelse(treatment %in% c("insects.molluscs.rabbits.limed.min.grass.plus.p", "insects.molluscs.rabbits.limed.min.herb.plus.p"), "P*lime*herbicide",
                                                                                                                                                                                                                                                  ifelse(treatment %in% c("insects.molluscs.rabbits.limed.min.grass.plus.k", "insects.molluscs.rabbits.limed.min.herb.plus.k"), "K*lime*herbicide",
                                                                                                                                                                                                                                                         ifelse(treatment %in% c("insects.molluscs.rabbits.limed.min.grass.plus.mg", "insects.molluscs.rabbits.limed.min.herb.plus.mg"), "Mg*lime*herbicide",
                                                                                                                                                                                                                                                                ifelse(treatment %in% c("insects.molluscs.rabbits.limed.min.grass.all.nutr", "insects.molluscs.rabbits.limed.min.herb.all.nutr",
                                                                                                                                                                                                                                                                                        "insects.molluscs.rabbits.limed.min.grass.plus.pk", "insects.molluscs.rabbits.limed.min.herb.plus.pk",
                                                                                                                                                                                                                                                                                        "insects.molluscs.rabbits.limed.min.grass.min.pk", "insects.molluscs.rabbits.limed.min.herb.min.pk",
                                                                                                                                                                                                                                                                                        "insects.molluscs.rabbits.limed.min.grass.min.n", "insects.molluscs.rabbits.limed.min.herb.min.n",
                                                                                                                                                                                                                                                                                        "insects.molluscs.rabbits.limed.min.grass.min.p", "insects.molluscs.rabbits.limed.min.herb.min.p",
                                                                                                                                                                                                                                                                                        "insects.molluscs.rabbits.limed.min.grass.min.k", "insects.molluscs.rabbits.limed.min.herb.min.k",
                                                                                                                                                                                                                                                                                        "insects.molluscs.rabbits.limed.min.grass.min.mg", "insects.molluscs.rabbits.limed.min.herb.min.mg"), "mult_nutrient*lime*herbicide",
                                                                                                                                                                                                                                                                       ifelse(treatment %in% c("spray.pellets.fenced.limed.min.grass.plus.n",
                                                                                                                                                                                                                                                                                               "spray.pellets.rabbits.limed.min.grass.plus.n",
                                                                                                                                                                                                                                                                                               "spray.molluscs.fenced.limed.min.grass.plus.n",
                                                                                                                                                                                                                                                                                               "spray.molluscs.fenced.limed.min.grass.plus.n",
                                                                                                                                                                                                                                                                                               "insects.pellets.fenced.limed.min.grass.plus.n",
                                                                                                                                                                                                                                                                                               "insects.pellets.rabbits.limed.min.grass.plus.n",
                                                                                                                                                                                                                                                                                              "insects.molluscs.rabbits.limed.min.grass.plus.n"), "N*herb_removal*lime*herbicide",
                                                                                                                                                                                                                                                                              ifelse(treatment %in% c("spray.pellets.fenced.limed.min.grass.plus.p",
                                                                                                                                                                                                                                                                                                      "spray.pellets.rabbits.limed.min.grass.plus.p",
                                                                                                                                                                                                                                                                                                      "spray.molluscs.fenced.limed.min.grass.plus.p",
                                                                                                                                                                                                                                                                                                      "spray.molluscs.fenced.limed.min.grass.plus.p",
                                                                                                                                                                                                                                                                                                      "insects.pellets.fenced.limed.min.grass.plus.p",
                                                                                                                                                                                                                                                                                                      "insects.pellets.rabbits.limed.min.grass.plus.p",
                                                                                                                                                                                                                                                                                                      "insects.molluscs.rabbits.limed.min.grass.plus.p"), "P*herb_removal*lime*herbicide",
                                                                                                                                                                                                                                                                                     ifelse(treatment %in% c("spray.pellets.fenced.limed.min.grass.plus.k",
                                                                                                                                                                                                                                                                                                             "spray.pellets.rabbits.limed.min.grass.plus.k",
                                                                                                                                                                                                                                                                                                             "spray.molluscs.fenced.limed.min.grass.plus.k",
                                                                                                                                                                                                                                                                                                             "spray.molluscs.fenced.limed.min.grass.plus.k",
                                                                                                                                                                                                                                                                                                             "insects.pellets.fenced.limed.min.grass.plus.k",
                                                                                                                                                                                                                                                                                                             "insects.pellets.rabbits.limed.min.grass.plus.k",
                                                                                                                                                                                                                                                                                                             "insects.molluscs.rabbits.limed.min.grass.plus.k"), "K*herb_removal*lime*herbicide",
                                                                                                                                                                                                                                                                                            ifelse(treatment %in% c("spray.pellets.fenced.limed.min.grass.plus.mg",
                                                                                                                                                                                                                                                                                                                    "spray.pellets.rabbits.limed.min.grass.plus.mg",
                                                                                                                                                                                                                                                                                                                    "spray.molluscs.fenced.limed.min.grass.plus.mg",
                                                                                                                                                                                                                                                                                                                    "spray.molluscs.fenced.limed.min.grass.plus.mg",
                                                                                                                                                                                                                                                                                                                    "insects.pellets.fenced.limed.min.grass.plus.mg",
                                                                                                                                                                                                                                                                                                                    "insects.pellets.rabbits.limed.min.grass.plus.mg",
                                                                                                                                                                                                                                                                                                                    "insects.molluscs.rabbits.limed.min.grass.plus.mg"), "Mg*herb_removal*lime*herbicide", "mult_nutrient*herb_removal_lime*herbicide"))))))))))))))))))))))))))))))))))))))))))))))))%>%
  unique()

ton <- read.csv("SIU_TON.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c("AC", "AS", "AB"), 9.8, 
                  ifelse(treatment %in% c("1C", "1S", "1B") & treatment_year %in% c(1,6,11,17,22),2.8,0)),
         p=ifelse(treatment %in% c("AC", "AS", "AB"), 3.97, 
                  ifelse(treatment %in% c("1C", "1S", "1B") & treatment_year %in% c(1,6,11,17,22),2.8,0)), 
         k=ifelse(treatment %in% c("AC", "AS", "AB"), 8.74, 
                  ifelse(treatment %in% c("1C", "1S", "1B") & treatment_year %in% c(1,6,11,17,22),2.8,0)), 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=ifelse(treatment %in% c("CC", "1C", "AC"), 0, 1), 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=ifelse(treatment == "1C", "fertilized every 5 years",
                            ifelse(treatment %in% c("CS","AS"), "spring mowing", 
                                   ifelse(treatment == "1S", "fertilized every 5 years + spring mowing", 
                                          ifelse(treatment %in% c("CB", "AB"), "spring and fall mowing", 
                                                 ifelse(treatment == "1B", "fertilized every 5 years + spring and fall mowing", 0))))),
         successional=0,
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment %in% c('1C', 'AC', 'CS', 'CB'), 1, ifelse(treatment == "CC", 0, 2)))%>%
  mutate(resource_mani=ifelse(treatment %in% c("CS", "CB"), 0, 1))%>%
  mutate(max_trt= ifelse(treatment %in% c('CC', 'AC', 'AS', 'CB', '1B', 'AB'), 1, 0))%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment== "CC", "control", 
                         ifelse(treatment %in% c("CS", "CB"), "mow_clip",
                                ifelse(treatment %in% c("1C, AC"), "mult_nutrient", "mult_nutrient*mow_clip"))))%>%
  unique()

uk<-read.delim("SKY_UK.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment %in% c('C','H'), 0, 30),
         temp=ifelse(treatment %in% c('C','P'), 0, 3),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=1, 
         plant_mani=1,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment=='HP', 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=='H', 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='C', 'control', ifelse(treatment=='H', 'temp', ifelse(treatment=='P', 'irr', 'irr*temp'))))%>%
  unique()

clima <- read.csv("SORBAS_CLIMARID.csv") %>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type = 0,
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0,
         p=0, 
         k=0, 
         CO2=0, 
         precip=ifelse(treatment %in% c("RR", "W+RR"), -30, 0), 
         temp=ifelse(treatment %in% c("W", "W+RR"), 3, 0),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=1, 
         plant_mani=0,  
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=="C", 0, ifelse(treatment=='W+RR',2,1)))%>%
  mutate(resource_mani= ifelse(treatment == "W", 0,1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=="C", "control", ifelse(treatment == "RR", "drought", ifelse(treatment == "W", "temp", "drought*temp"))))%>%
  unique()
  
  
nitrogen<-read.csv("SR_Nitrogen.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment, community_type)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('1_NITROGEN','0_NITROGEN'), 4, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=1, 
         plant_mani=ifelse(treatment %in% c('1_CONTROL','1_NITROGEN'),1,0),  
         plant_trt=ifelse(treatment %in% c('1_CONTROL','1_NITROGEN'),1,0),
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='0_CONTROL', 0, ifelse(treatment=='1_NITROGEN',2,1)))%>%
  mutate(resource_mani=ifelse(treatment=='1_CONTROL',0,1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='0_CONTROL', 'control', ifelse(treatment=='0_NITROGEN', 'N', ifelse(treatment=='1_CONTROL', 'plant_mani', 'N*plant_mani'))))%>%
  unique()

water<-read.csv("SR_Water.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment, community_type)%>%
  mutate(nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0,
         p=0,
         k=0, 
         CO2=0, 
         precip=ifelse(treatment %in% c('0_WATER_0','0_WATER_1','1_WATER_0','1_WATER_1'), 34.1,0), 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=ifelse(treatment %in% c('0_CONTROL_0','1_CONTROL_0',"0_WATER_0",'1_WATER_0'), 1,0),
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=1, 
         plant_mani=ifelse(treatment %in% c('1_CONTROL_0','1_CONTROL_1',"1_WATER_0",'1_WATER_1'),1,0), 
         plant_trt=ifelse(treatment %in% c('1_CONTROL_0','1_CONTROL_1',"1_WATER_0",'1_WATER_1'),1,0),
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='0_CONTROL_0', 0, ifelse(treatment=='1_WATER_1',3,ifelse(treatment %in% c('1_CONTROL_0','0_CONTROL_1','0_WATER_0'),1,2))))%>%
  mutate(resource_mani=ifelse(treatment %in% c('1_CONTROL_1','0_CONTROL_0',"1_CONTROL_0"),0,1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='0_CONTROL_0', 'control', ifelse(treatment=='0_CONTROL_1', 'grazed', ifelse(treatment=='0_WATER_0', 'irr', ifelse(treatment=='1_CONTROL_1', 'plant_mani*grazed', ifelse(treatment=='1_CONTROL_0', 'plant_mani', ifelse(treatment=='1_WATER_1', 'irr*plant_mani*grazed', ifelse(treatment=='1_WATER_0', 'irr*plant_mani', 'irr*grazed'))))))))%>%
  unique()

gane<-read.delim("SVA_GANE.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment %in% c('C','P'), 0, ifelse(treatment %in% c('LN','LNP'), 0.5, 5)),
         p=ifelse(treatment %in% c('P','LNP','HNP'), 1, 0),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment %in% c('LN','HN','P'), 1, 2)))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=ifelse(treatment %in% c('C','HN','P','HNP'), 1, 0))%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='C', 'control', ifelse(treatment=='P', 'P', ifelse(treatment %in% c('HN','LN'), 'N', 'N*P'))))%>%
  unique()

tface<-read.csv("TAS_FACE.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0,
         nutrients=0, light=0, carbon=1, water=0, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=ifelse(treatment %in% c('UnwarmedFACE',"WarmedFACE"), 170,0), 
         precip=0, 
         temp=ifelse(treatment %in% c('WarmedControl',"WarmedFACE"),2,0),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='UnwarmedControl', 0, ifelse(treatment=='WarmedFACE',2,1)))%>%
  mutate(resource_mani=ifelse(treatment=='WarmedControl',0,1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='UnwarmedControl', 'control', ifelse(treatment=='UnwarmedFACE', 'CO2', ifelse(treatment=='WarmedControl', 'temp', 'CO2*temp'))))%>%
  unique()

lovegrass<-read.csv("TRA_Lovegrass.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0,
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment %in% c('gcc','ghc','gsc',"ncc",'nhc','nsc'), 0, 0.432), 
         p=ifelse(treatment %in% c('gcc','ghc','gsc',"ncc",'nhc','nsc'), 0, 0.022), 
         k=ifelse(treatment %in% c('gcc','ghc','gsc',"ncc",'nhc','nsc'), 0, 0.082), 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=ifelse(treatment %in% c('gsc','gsn','nsc',"nsn"),1, 0), 
         burn=0, 
         herb_removal=ifelse(treatment %in% c('ncc','ncn','nhc',"nhn","nsc",'nsn'),1,0),
         management=0,
         other_trt=0,
         trt_details=0, 
         successional=0, 
         plant_mani=ifelse(treatment %in% c('ghc','ghn','nhc','nhn'),1,0), 
         plant_trt=ifelse(treatment %in% c('ghc','ghn','nhc','nhn'),1,0),
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='ncc',0,ifelse(treatment %in% c('nhc','nsc','gcc'),1, ifelse(treatment %in% c('ghc','gsc'),2,ifelse(treatment=="ncn",3,ifelse(treatment %in% c('ghn','gsn'),5,4))))))%>%
  mutate(resource_mani=ifelse(treatment %in% c('ghc','gsc',"ncc",'nhc','nsc'), 0, 1))%>%
  mutate(max_trt=1)%>%
  mutate(public=0)%>%
  mutate(factorial=1)%>%
  mutate(trt_type=ifelse(treatment=='ncc', 'control', ifelse(treatment=='ncn', 'mult_nutrient', ifelse(treatment=='nhc', 'plant_mani', ifelse(treatment=='nhn', 'mult_nutrient*plant_mani', ifelse(treatment=='nsc', 'mow_clip', ifelse(treatment=='nsn', 'mult_nutrient*mow_clip', ifelse(treatment=='gcc', 'grazed', ifelse(treatment=='gcn', 'mult_nutrient*grazed', ifelse(treatment=='ghc', 'plant_mani*grazed', ifelse(treatment=='ghn', 'mult_nutrient*plant_mani*grazed', ifelse(treatment=='gsc', 'grazed*mow_clip', 'mult_nutrient*grazed*mow_clip'))))))))))))%>%
  unique()

edge<-read.csv('USA_EDGE.csv')%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment, community_type)%>%
  mutate(nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=0, 
         precip=ifelse(treatment=="chr", -66, ifelse(treatment=="int", -100, 0)), 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=ifelse(treatment=="del", "monsoon rain applied later", 0), 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='con', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=0)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='con', 'control', ifelse(treatment=='del', 'precip_vari', 'drought')))%>%
  unique()
edge <- edge[edge$site_code != "SEV",]

sedge <- read.csv('SEV_EDGE20.csv')%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment, community_type)%>%
  mutate(nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=0, 
         precip=ifelse(treatment=="E", -66, 0), 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=ifelse(treatment=="D", "monsoon rain applied later", 0), 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=0)%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='C', 'control', ifelse(treatment=='D', 'precip_vari', 'drought')))%>%
  unique()

nitadd<-read.csv("YMN_NitAdd.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0,
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N0',0, ifelse(treatment=='N5', 2.3, ifelse(treatment=='N10',4.7, ifelse(treatment=='N20',9.3, ifelse(treatment=='N40',18.7,37.3))))), 
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         management=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         plant_trt=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N0', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(max_trt=ifelse(treatment=='N0'|treatment=='N80', 1, 0))%>%
  mutate(public=0)%>%
  mutate(factorial=0)%>%
  mutate(trt_type=ifelse(treatment=='N0', 'control', 'N'))%>%
  unique()

###merge all datasets
combine<-rbind(bffert, bgp, biocon, bowman, ccd, clip, clima, clonal, culardoch, cxn, e001, e002, 
               e2, e6, edge, eelplot, events, exp1, face, fert1, fert2, fireplots, gane, gap2, gb, 
               gce, gcme, gcme2, gfp, grazeprecip, herbdiv, herbwood, hprecip, imagine, interaction, 
               irg, kgfert, lind, lovegrass, lucero, mat2, megarich, mnt, mwatfer, nash, nde, nfert, 
               nitadd, nitphos,  nitrogen,nsfc, nsfc2, nut, nutnet, oface, pennings, phace, pme, precip, 
               pplots, pq, ramps, rhps, rmapc, sedge, snfert, snow, sprecip, study119, study278, t7, 
               ter, tface,tide,tmece,ton, uk, wapaclip, warmnit, warmnut, water, watering, wenndex, wet, yu)

write.csv(combine, "~/Dropbox/CoRRE_database/Data/CompiledData/ExperimentInfo.csv")


temp_df <- unique(combine[,c(1,2,6,33)])
trt_sum <- as.data.frame(table(temp_df$trt_type))
trt_sum <- trt_sum[order(trt_sum$Freq, decreasing = TRUE),]
names(trt_sum) <- c("treatment_type", "Number_of_experiments")
write.csv(trt_sum, "~/Dropbox/CoRRE_database/Data/CompiledData/treatment_type_summary.csv")

# Creating a summary of treatment types by location - all treatment types
df <- unique(combine[,c(1,2,29,33)])
temp <- df %>% count(trt_type, site_code)
temp1<- temp[,c(1,2)] %>% count(trt_type)
names(temp1) <- c("treatment_type","Number_of_locations")
test <- merge(trt_sum, temp1)
write.csv(test, "~/Dropbox/CoRRE_database/Data/CompiledData/treatment_type_summary_location.csv")

# Larger grouping of treatment types
# CO2, N, P, Drought, Irr, Var
df1 <- df %>% mutate(trt_type = ifelse(trt_type %in% c("N*P"), "mult_nutrient",
                                       ifelse(trt_type %in% c("mow_clip", "herb_removal", "plant_mani", "lime", "burn", "seed","fungicide", "disturbance", "herbicide", "stone","till"), "other_non_resource", 
                                              ifelse(trt_type=="N*irr*CO2", "3 resources", 
                                                     ifelse(trt_type == "N*P", "mult_nutrient", 
                                                            ifelse(trt_type %in% c("irr*herb_removal", "irr*mow_clip", "irr*plant_mani", "irr*plant_mani*herb_removal"), "irr*non_resource", 
                                                                   ifelse(trt_type %in% c("mult_nutrient*herb_removal", "mult_nutrient*mow_clip", "mult_nutrient*plant_mani", "N*P*burn", "N*P*mow_clip", "N*P*seed", "mult_nutient*fungicide", "mult_nutrient*herb_removal_lime*herbicide", "mult_nutrient*herb_removal*lime", "mult_nutrient*herb_removal*lime", "mult_nutrient*herb_removal*mow_clip", "mult_nutrient*herbicide", "mult_nutrient*lime", "mult_nutrient*lime*herbicide", "mult_nutrient*plant_mani*herb_removal", "N*P*burn*graze", "N*P*burn*mow_clip", "mult_nutrient*herb_removal*herbicide", "mult_nutrient*fungicide"), "mult_nutrient*non_resource",
                                                                          ifelse(trt_type %in% c("N*irr*temp", "N*P*temp", "drought*CO2*temp", "irr*CO2*temp", "mult_nutrient*temp", "N*CO2*temp", "N*irr*CO2*temp"), "mult_resources*temp", 
                                                                                 ifelse(trt_type %in% c("N*till", "N*temp*fungicide", "N*stone", "N*seed*mow_clip", "N*plant_mani*disturbance", "N*lime*herbicide", "N*lime", "H*herbicide", "N*herb_removal*lime*herbicide", "N*herb_removal*lime", "N*herb_removal)herbicide", "N*herb_removal", "N*fungicide", "N*disturbance", "N*burn*graze", "N*burn*mow_clip", "N*seed", "N*burn", "N*mow_clip", "N*plant_mani", "N*herbicide", "N*herb_removal*herbicide"), "N*other_non_resource(s)", 
                                                                                        ifelse(trt_type %in% c("P*burn", "P*seed", "P*burn*graze", "P*burn*mow_clip", "P*herb_removal", "P*herb_removal*herbicide", "P*herb_removal*lime", "P*herb_removal*lime*herbicide", "P*herbicide", "P*lime", "P*lime*herbicide", "P*mow_clip"), "P*other_non_resource(s)",
                                                                                               ifelse(trt_type %in% c("plant_mani*herb_removal", "burn*mow_clip", "burn*graze", "herb_removal*herbicide","herb_removal*lime", "herb_removal*lime*herbicide", "herb_removal*mow_clip", "lime*herbicide", "seed*mow_clip", "plant_mani*disturbance"), "mult_non_resource", 
                                                                                                      ifelse(trt_type %in% c("K", "C", "Ca addition", "Mg", "protein"), "other_nutrient", 
                                                                                                             ifelse(trt_type %in% c("temp*mow_clip", "temp*fungicide"), "temp*non_resource", 
                                                                                                                    ifelse(trt_type %in% c("C*stone", "*herb_removal", "K*herb_removal", "K*herb_removal*herbicide", "K*herb_removal*lime", "K*herb_removal*lime*herbicide", "K*herbicide", "K*lime*herbicide","Mg*herb_removal", "Mg*herb_removal*herbicide", "Mg*herb_removal*lime", "Mg*herb_removal*lime*herbicide", "Mg*herbicide", "Mg*lime","Mg*lime*herbicide"), "other_nutrient*non_resource(s)", 
                                                                                                                           trt_type))))))))))))))
df2 <- temp_df %>% mutate(trt_type = ifelse(trt_type %in% c("N*P"), "mult_nutrient",
                                       ifelse(trt_type %in% c("mow_clip", "herb_removal", "plant_mani", "lime", "burn", "seed","fungicide", "disturbance", "herbicide", "stone","till"), "other_non_resource", 
                                              ifelse(trt_type=="N*irr*CO2", "3 resources", 
                                                     ifelse(trt_type == "N*P", "mult_nutrient", 
                                                            ifelse(trt_type %in% c("irr*herb_removal", "irr*mow_clip", "irr*plant_mani", "irr*plant_mani*herb_removal"), "irr*non_resource", 
                                                                   ifelse(trt_type %in% c("mult_nutrient*herb_removal", "mult_nutrient*mow_clip", "mult_nutrient*plant_mani", "N*P*burn", "N*P*mow_clip", "N*P*seed", "mult_nutient*fungicide", "mult_nutrient*herb_removal_lime*herbicide", "mult_nutrient*herb_removal*lime", "mult_nutrient*herb_removal*lime", "mult_nutrient*herb_removal*mow_clip", "mult_nutrient*herbicide", "mult_nutrient*lime", "mult_nutrient*lime*herbicide", "mult_nutrient*plant_mani*herb_removal", "N*P*burn*graze", "N*P*burn*mow_clip", "mult_nutrient*herb_removal*herbicide", "mult_nutrient*fungicide"), "mult_nutrient*non_resource",
                                                                          ifelse(trt_type %in% c("N*irr*temp", "N*P*temp", "drought*CO2*temp", "irr*CO2*temp", "mult_nutrient*temp", "N*CO2*temp", "N*irr*CO2*temp"), "mult_resources*temp", 
                                                                                 ifelse(trt_type %in% c("N*till", "N*temp*fungicide", "N*stone", "N*seed*mow_clip", "N*plant_mani*disturbance", "N*lime*herbicide", "N*lime", "H*herbicide", "N*herb_removal*lime*herbicide", "N*herb_removal*lime", "N*herb_removal)herbicide", "N*herb_removal", "N*fungicide", "N*disturbance", "N*burn*graze", "N*burn*mow_clip", "N*seed", "N*burn", "N*mow_clip", "N*plant_mani", "N*herbicide", "N*herb_removal*herbicide"), "N*other_non_resource(s)", 
                                                                                        ifelse(trt_type %in% c("P*burn", "P*seed", "P*burn*graze", "P*burn*mow_clip", "P*herb_removal", "P*herb_removal*herbicide", "P*herb_removal*lime", "P*herb_removal*lime*herbicide", "P*herbicide", "P*lime", "P*lime*herbicide", "P*mow_clip"), "P*other_non_resource(s)",
                                                                                               ifelse(trt_type %in% c("plant_mani*herb_removal", "burn*mow_clip", "burn*graze", "herb_removal*herbicide","herb_removal*lime", "herb_removal*lime*herbicide", "herb_removal*mow_clip", "lime*herbicide", "seed*mow_clip", "plant_mani*disturbance"), "mult_non_resource", 
                                                                                                      ifelse(trt_type %in% c("K", "C", "Ca addition", "Mg", "protein"), "other_nutrient", 
                                                                                                             ifelse(trt_type %in% c("temp*mow_clip", "temp*fungicide"), "temp*non_resource", 
                                                                                                                    ifelse(trt_type %in% c("C*stone", "*herb_removal", "K*herb_removal", "K*herb_removal*herbicide", "K*herb_removal*lime", "K*herb_removal*lime*herbicide", "K*herbicide", "K*lime*herbicide","Mg*herb_removal", "Mg*herb_removal*herbicide", "Mg*herb_removal*lime", "Mg*herb_removal*lime*herbicide", "Mg*herbicide", "Mg*lime","Mg*lime*herbicide"), "other_nutrient*non_resource(s)", 
                                                                                                                           trt_type))))))))))))))



temp <- df1 %>% count(trt_type, site_code)
temp1<- temp[,c(1,2)] %>% count(trt_type)
names(temp1)<- c("treatment_type", "Number_of_locations")
trt_sum2 <- as.data.frame(table(df2$trt_type))
names(trt_sum2) <- c("treatment_type", "Number_of_experiments")
trt_sum3 <- merge(trt_sum2, temp1)
#trt_sum3 <- trt_sum3[trt_sum3$treatment_type != "control",]
write.csv(trt_sum3, "~/Dropbox/CoRRE_database/Data/CompiledData/treatment_type_summary_broad_groups.csv")



