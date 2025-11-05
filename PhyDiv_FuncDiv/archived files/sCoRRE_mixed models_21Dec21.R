### Running mixed models comparing FDis and PD across treatments and sites
###
### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Created: Dec 16, 2021; last updated: Dec 16, 2021

### R version: 4.1.1

#library(FD)
library(tidyverse)
library(nlme)
#library(beepr)

# rm(list=ls())


#setwd("~/Dropbox/sDiv_sCoRRE_shared/")
#setwd("/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/") #Padu's wd
setwd( "C:\\Users\\wilco\\Dropbox\\shared working groups\\sDiv_sCoRRE_shared\\paper 2_PD and FD responses\\data\\") # Kevin's laptop wd

# Standard Error Function:
# se <- function(x, na.rm=na.rm){
#   SE=sd(x,na.rm=TRUE)/sqrt(length(x))
#   return(SE)
# }

###
### Read in and prep all data
###

### Read in and prep data
div_metrics_raw <- read.csv("CoRRE_allDiversityMetrics_phyFunAnalysis.csv") %>%
  filter(treatment_year==experiment_length) %>%
  mutate(site_proj_comm = paste(site_code, project_name, community_type, sep="::")) %>%
  dplyr::select(site_code:plot_id, site_proj_comm, treatment, trt_type, trt_type2, plot_mani, experiment_length:anpp, MAP, MAT,
                FDis, mntd.ses) %>%
  mutate(trt_binary = ifelse(plot_mani==0,"Control","Treatment"))

fdis_model_full <- lme(FDis ~ trt_binary*trt_type2
    , data=div_metrics_raw
    , random = ~1 |site_proj_comm
    , na.action = na.omit)

### Error in MEEM(object, conLin, control$niterEM) : 
### Singularity in backsolve at level 0, block 1


