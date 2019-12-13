library(picante)

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
  unique(my.species.list[which(x == my.species.list$genus_species), 2])
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















