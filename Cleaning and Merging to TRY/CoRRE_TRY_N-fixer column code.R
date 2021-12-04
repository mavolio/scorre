### IDENTIFYING N FIXERS IN CORRE-TRY DATABASE ####
### Ben Taylor 12/10/2019 ###

################################################################################
### Read in data sheets and get them into shape ################################
################################################################################

setwd("C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\CoRRE_database\\Data\\TRYCoRREMerge\\CORRE-TRY N Fixers")

library(tidyverse)

### CORRE-TRY DATABASE ###
ct_orig<-read.csv("CoRRE_new spp.csv")%>%
  rbind(read.csv("CoRRE_new spp_2021add.csv"))

ct<-ct_orig
ct.nms<-as.data.frame(unique(ct[,1])) #isolates just the scientific names for matching
names(ct.nms)[1]<-"scientific_name" #renames the species column to match other dataframes

################################################################################
### GRIN/SPRENT DATABASE ###
grin<-read.csv("GRIN_nodulation_reports.csv") #Database of species 
grin[,1]<-as.character(grin[,1])
grin[,2]<-as.character(grin[,2])
grin$scientific_name<-paste(grin$Genus, grin$Species)

#making the Conspecs list for the GRIN database
gn<-unique(grin$Genus)
conspp<-NULL
for(i in 1:length(gn)){
  print(i)
  pass<-grin[grin$Genus==gn[i],]
  pass$Conspecs<-(sum(pass$fixer)/nrow(pass))
  conspp<-rbind(conspp,pass)
}
grin<-conspp
#Making the fixer/non-fixer column for GRIN
fixcut <- 0.6
grin$conspp.fixer<-ifelse(grin$Conspecs>=fixcut, 1,0)
grin$fix.lib<-ifelse(grin$fixer+grin$conspp.fixer<1, 0, 1)
grin.lib<-grin[,c("scientific_name", "fix.lib")]

################################################################################
### WERNER ET AL. 2014 DATABASE ####
wern<-read.csv("Werner_Nfix.csv")
wern<-wern[,c(1,3,6,12)]
colnames(wern)<-c("scientific_name","Legume","fixer","genus")
wern<-wern[!duplicated(wern$scientific_name),]
wern$fixer<-ifelse(wern$fixer=="Yes", 1,0)
wern$Legume<-ifelse(wern$Legume=="Yes", 1,0)

################################################################################
### Isolating Lists of N-fixer Species names from each database ###
grn.fix.nms<-as.factor(as.character(unique(grin[grin$fixer==1,]$scientific_name)))
grn.fix.nms.lib<-as.factor(as.character(unique(grin[grin$fixer==1,]$scientific_name)))
wern.fix.nms<-as.factor(as.character(unique(wern[wern$fixer==1,]$scientific_name)))
wern.act.nms<-as.factor(as.character(unique(wern[wern$fixer==1&wern$Legume==0,]$scientific_name)))

################################################################################
### CONCATENATED LIST OF ALL N-FIXING SPECIES (MOST LIBERAL DEFINITION) ###
fix.nms<-sort(unique(c(as.character(grn.fix.nms),as.character(grn.fix.nms.lib), 
           as.character(wern.fix.nms))))

################################################################################
### Adding N-fixer columns to CoRRE-TRY Database ###
ct_orig$fixer<-ifelse(ct_orig$species_matched%in%fix.nms, yes=1,no=0) #Includes a 1 for all N-fixers from both sources                          
ct_orig$act<-ifelse(ct_orig$species_matched%in%wern.act.nms, yes=1,no=0) #Includes a 1 for actinorhizal N-fixers
ct_orig$grin<-ifelse(ct_orig$species_matched%in%grn.fix.nms, yes=1,no=0) #Includes a 1 for N-fixers from GRIN
ct_orig$wern<-ifelse(ct_orig$species_matched%in%wern.fix.nms, yes=1,no=0) #Includes a 1 for N-fixers from WERNER

write.csv(ct_orig, file="CoRRE_TRY_species_list_N-fixers.csv", row.names=F) #write file to .csv

#find just new species from 2021 update to CoRRE 2.0

new<-read.csv("CoRRE_new spp_2021add.csv")%>%
  select(species_matched)%>%
  mutate(version=2)%>%
  full_join(ct_orig)%>%
  filter(version==2)
write.csv(new, file="CoRRE_TRY_species_list_N-fixers_2021add.csv", row.names=F)


################################################################################
### Simple Summary Stats ###
sum(ct_orig$fixer) #446 total N-fixers in the CoRRE-TRY database
sum(ct_orig$act) #22 actinorhizal N-fixers in the CoRRE-TRY database
sum(ct_orig$grin) #422 of the N-fixers were reported in the GRIN database
sum(ct_orig$wern) #435 of the N-fixers were reported in the WERNER 2014 database

#Species named as fixers in GRIN but not in WERNER 2014
grin.not.wern<-ct_orig[ct_orig$grin==1&ct_orig$wern==0,]

#Species named as fixers in WERNER 2014 but not in GRIN
wern.not.grin<-ct_orig[ct_orig$grin==0&ct_orig$wern==1,]
