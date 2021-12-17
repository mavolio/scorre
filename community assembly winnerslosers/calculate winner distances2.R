
#load packages
library(gawdis)
library(ape)
library(dplyr) 

#set directory
my.wd <- "~/Dropbox/sCoRRE/sDiv_sCoRRE_shared"
my.wd <- "/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/" #Padu's wd
# my.wd <- "~/Dropbox/sCoRRE/sDiv_sCoRRE_shared"

##########################################
# Get functional tree (functional dissimilarities)
##########################################

#load trait data:
trait<-read.table(paste(my.wd, "/CoRRE data/CoRRE data/trait data/Final Cleaned Traits/Continuous_Traits/Backtrans_GapFilled_sCorre.csv", sep=""), header=T, sep=",", fill = TRUE)

#get mean trait values:
trait <- trait %>%  group_by(species_matched) %>% summarise_at(vars("seed_dry_mass","stem_spec_density", "leaf_N",                 
                                                                  "leaf_P", "LDMC", "leaf_C","leaf_dry_mass",          
                                                                   "plant_height_vegetative", "leaf_C.N", "SLA","water_content",          
                                                                   "rooting_depth", "SRL"), median, na.rm=T) %>% as.data.frame() 
rownames(trait)<-trait$species_matched
trait$species_matched<-NULL

#calculate functional dissimilarities between species using GAWDIS:
#I haven't log-transformed traits. Pending!!
#trait.dis<-gawdis(trait, w.type="optimized") #This took forever and I had to stop it. Can remove "w.type="optiimized" but a warning pops up)
#trait.dis<-gawdis(trait) #with "analytical" method
trait.dis<-cluster::daisy(trait, metric="gower") #or with gower distance

#convert distance object into phylotree using the Ward's algorithm (https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust):
ftree<-as.phylo(hclust(trait.dis, method="average"))

#save output:
write.tree(ftree, paste(my.wd, "ftree.scorre.gowdis.log.upgma.tre", sep=""))

#get dissimilarity matrix between species:
trait.dis.matrix<-as.matrix(trait.dis)
diag(trait.dis.matrix) <- NA

##############################################
# identify winners and losers per experiment #
##############################################
# DCi though time for control and treatment and then slopes in Adam's framework + extinct and immigration species
# important!remove too rare species; maybe identify via CSi values, do histograms to identify threshlds but it could be around 0.1


# Emily developed code to generate DCi through time: “emily ambient v2.R”






##############################################
# calculate mdns, dnns, wmdns
##############################################
# only experimennts that are at least 5 years

# first step in one experiment
#ftree<-read.tree(paste(my.wd, "/CoRRE data/CoRRE data/trait data/ftree.scorre.gowdis.log.upgma.tre", sep=""))
dis<-trait.dis.matrix
#as.matrix(cophenetic.phylo(ftree))


corre<-read.table(paste(my.wd, "/CoRRE data/CoRRE data/community composition/CoRRE_RelativeCover_Dec2021.csv", sep=""), header=T, sep=",", fill = TRUE)
#take example experiment suggested by Adam
CDR_tmp<-subset(corre, site_code=="CDR" & project_name=="e001" & community_type%in%c("C") & treatment%in%c("1", "8"))

#load species to use:
spp<-read.table(paste(my.wd, "/CoRRE data/CoRRE data/trait data/FullList_Nov2021.csv", sep=""), header=T, sep=",", fill = TRUE)[,c(2,5)]

CDR<-merge(CDR_tmp, spp, by="genus_species", all.x=T)

dim(CDR)


# here we go with a visual selection with help from Adam for a focal experiment
win_los<-c("Elymus repens", "Schizachyrium scoparium")
# for each
CDR<-unique(CDR$plot_id)
treat<-unique(CDR[,c(7,9)])$treatment

species<- unique(CDR$species_matched)

# calculate for all species distance to the rest of the community 
out2<-NULL
for(i in 1:length(plots)){
  print(plots[i])
  corre.filt<-subset(CDR, plot_id==plots[i])
  years<-unique(corre.filt$calendar_year)
  out<-NULL
  for(j in 1:length(years)){
    print(years[j])
    corre.filtered<-subset(corre.filt, calendar_year==years[j])
    species.names<-unique(corre.filtered$species_matched)
    #add winner/loser in case it is not yet in the community
    species.names<-unique(c(species.names, win_los))
    
    if (length(species.names)>1) {
      trait.dis.sub<-dis[(rownames(dis) %in% species.names), (colnames(dis) %in% species.names)]
      
      a<-as.data.frame(rowMeans(trait.dis.sub, na.rm=T))
      a$species_matched<-rownames(a)
      names(a)[1]<-"mpd"
      a$year<-years[j]
      
      out<-rbind(out, a)
    }
  }
  out$plot_id<-plots[i]
  out$treatment<-treat[i]
  out2<-rbind(out2, out)
}
rownames(out2)<-NULL





write.table(out2, paste(my.wd, "mpd_win_los.csv", sep=""))
















#