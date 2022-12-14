################################################################################
##  7_sCoRRE_phylogeneticDiversity_CWMtraits.R: Examining community weighted mean trait responses to treatments in the CoRRE database.
##
##  Author: Kimberly Komatsu
##  Date created: December 13, 2022
################################################################################

library(data.table)
library(codyn)
library(fixest)
library(lme4)
library(tidyverse)


#functional diversity data
fDiv <- read.csv('CoRRE_functionalDiversity_2022-12-13.csv')