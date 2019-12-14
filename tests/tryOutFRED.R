library(tidyverse)
library(data.table)

my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/"
sp <- read.csv("FRED (fine root database)/FRED-species.csv", header=F, sep=";")
sp$sp <- paste(sp[,1],sp[,2])
we <- read.csv("CoRRE data/CoRRE_TRY_species_list.csv", header=T, sep=",")
length(intersect(we$species_matched, sp$sp))
