library(picante)
library(tidyverse)
library(gridExtra)
library(data.table)

my.wd <- "C:/Users/mavolio2/Dropbox/sDiv_sCoRRE_shared/"
my.wd <- "~/Dropbox/sDiv_sCoRRE_shared/"

#read in data and phylogeny
my.dat.first <- read.csv(paste(my.wd, "CoRRE data/CoRRE_raw_abundance_Nov2019.csv",sep=""))
my.tree <- read.tree(paste(my.wd, "CoRRE data/Phylogenies/scorre.tree.S3.tre",sep=""))
my.species.list <- read.csv(paste(my.wd, "CoRRE data/CoRRE_TRY_species_list.csv",sep=""))


#--------------------------------
# preparation
#--------------------------------

#select a sub-dataset
my.data <- my.dat.first %>%
  filter(site_code=="KNZ" & project_name=="pplots") %>%
  select(-X) %>%
  filter(treatment=="N1P0"|treatment=="N2P3")

#add a column with names that match phylogeny names
my.data$species_match <- unlist(lapply(my.data$genus_species, function(x) 
  unique(my.species.list[which(x == my.species.list$genus_species), "species_matched"])
  ))
head(my.data)

#check if it worked
str(my.tree)
my.tree$tip.label <- unlist(lapply(my.tree$tip.label, function(x){
  y <- strsplit(x, split="_")
  paste(y[[1]][1], y[[1]][2])
}))

my.species <- unique(my.data$species_match)
matches <- intersect(my.species, my.tree$tip.label)

# all good if TRUE
length(matches) == length(my.species)

#--------------------------------
# analyses
#--------------------------------
head(my.data)
unique(my.data$community_type)
my.data$community <- paste(my.data$treatment_year, my.data$treatment, my.data$plot_id, sep="_")


#shift from long to wide data.table format
wide.my.data <- dcast(my.data, community ~ species_match, value.var = "abundance")
rownames(wide.my.data) <- wide.my.data$community
wide.my.data <- as.matrix(wide.my.data[,-1])
wide.my.data <- apply(wide.my.data, 1:2, function(x) ifelse(is.na(x), 0, x))

my.tree.matrix <- cophenetic(my.tree)
colnames <- rownames <- my.tree$tip.label


#calculate mpd with null model per community
my.mpd <- ses.mpd(wide.my.data, my.tree.matrix)

new_columns <- t(sapply(rownames(my.mpd), function(x){
  z <- strsplit(x, split="_")
  c(treatment_year=z[[1]][1], treatment=z[[1]][2], plot_id=z[[1]][3])
}))
my.results <- cbind(new_columns, my.mpd)
my.results$treatment_year <- factor(my.results$treatment_year, levels = as.character(0:12))


#plot results (all communities are phylogenetically clustered)
ggplot(data=my.results, aes(x=treatment_year , y=mpd.obs.z)) +
  geom_boxplot(aes(color = as.factor(treatment))) 

ggplot(data=my.results, aes(x=treatment_year , y=mpd.obs.z)) +
  geom_violin(aes(fill = as.factor(treatment))) 

