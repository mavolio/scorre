
my.wd <- "~/Dropbox/sCoRRE/sDiv_sCoRRE_shared"


#  data 
ab_data <- read.csv(paste(my.wd, "/CoRRE data/CoRRE data/community composition/CoRRE_RawAbundance_Dec2021.csv",sep=""))
trait_data_cont <- read.csv(paste(my.wd, "/Trait Data/TRY Data/TRY Continuous data/TRY_trait_data_continuous_Nov2021.csv",sep=""))
trait_data_cont$species_matched <- tolower(trait_data_cont$species_matched)
trait_data_cont[1:10,1:10]


# only experimennts that are at least 5 years

##############################################
# identify winners and losers per experiment #
##############################################
# DCi though time for control and treatment and then slopes in Adam's framework + extinct and immigration species
# important!remove too rare species; maybe identify via CSi values, do histograms to identify threshlds but it could be around 0.1


# Emily developed code to generate DCi through time: “emily ambient v2.R”







##############################################
# calculate mdns, dnns, wmdns
##############################################

















#