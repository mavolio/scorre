
#load libraries:
library(ggtree)
library(ggplot2)
library(stringr)
library(plyr)
library(ape)
library(scico)
library(phytools)

#set directory.
my.wd<-"/Users/padulles/Documents/PD_MasarykU/sCoRRE/sCoRre/"

###
# Load data
###

#load species data (filtered after removing mosses and species missing from the phylogeny):
species.data<-read.table(paste(my.wd,"Species_DCiDiff_newtrts_filtered.csv",sep=""), header=T, sep=" ")

#load phylogenetic tree:
tree<-read.tree(paste(my.wd, "scorre.tree.win.los.tre", sep=""))


####
# Create empty data frame to store results:
###

result<-as.data.frame(tree$tip.label) #data frame with all species in our complete tree
names(result)[1]<-paste("species") #change name of first column


####
# Create loop to get significant DCi for all species in each treatment (can take a while to run)
####

trt<-levels(species.data$trt_type2) #get list of treatments
k<-3
for(k in 1:length(trt)){
  print(k/length(trt)) #keep track of the process
  
  # Filter by treatment
  dat<-subset(species.data, trt_type2==trt[k])[,c(1,4)] #select treatment
  dat<-aggregate(dat[, 2], list(dat$species_matched), mean, na.rm=T) #get mean DCi value per species
  dat<-dat[dat$Group.1 %in% tree$tip.label, ] #make sure all species in the data are on the tree
  rownames(dat)<-dat$Group.1 #set species names as rownames
  dat$Group.1<-NULL #and delete column with species names
  
  # Prune tree:
  tree2<-keep.tip(tree, rownames(dat))
  
  # Run function to calculate if each node has significantly higher or lower median DCi than
  # expected if phylogenetic relationships were at random (laod function in the "DCi_nodes_scorre.R" script)
  res<-node.mean(tree2, dat, 999)
 
  # Create table of significant nodes with direction of the effect
  significant<-res #create a copy of the main result
  significant$P_value[significant$P_value>0.05 | significant$SD_Exp<0.001]<-NA #replace non-significant with NA (and remove cases with very low SD; basal nodes)
  significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
  significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
  significant<-significant[complete.cases(significant), ] #remove non-significant nodes
  
  # Get list of species and assign if DCi was significantly greater or lower than random
  spp<-as.data.frame(tree2$tip.label) #create table
  names(spp)[1]<-paste("species") #change name of first column
  spp[,"DCi"] <- NA #create empty column
  
  # Go through each significant node, extract species and assign values in the "spp" table:
  for(i in 1:nrow(significant)){
    taxa<-tips(tree2, significant$Node[i]) #get list of taxa for a specific node
    #replace values in "spp" based on value
    for(j in 1:length(taxa)){
      spp$DCi[spp$species==taxa[j]] <- significant$P_value[i]
    }
  }
  spp$DCi[is.na(spp$DCi)] <- "random" #assign random to the other species
  names(spp)[2]<-trt[k]
  
  # Merge with result
  result<-merge(result, spp, by="species", all.x=T)
}
write.table(result, "species_sig_dci2.csv")



###
# Re-load result
##

result<-read.table("species_sig_dci2.csv") #re-load result
rownames(result)<-result$species #set species as rownames
result$species<-NULL #and delete column
colnames(result)<-c("All treatments", "+ CO2", "Disturbance", "Drought", "Herb removal", "Irrigation",
                    "N addition", "P addition", "+ Temperature")

###
# Assign families to nodes
###

#load family data and clean-up:
fam<-read.table(paste(my.wd, "species_families_2021.csv",sep=""), header=T, sep=",", fill = TRUE) #load data
fam$species_matched<-gsub(" ", "_", fam$species_matched) #adapt species nomenclature
fam$family[fam$family=="Compositae"]<-"Asteraceae" #replace family name
fam$family[fam$family=="Leguminosae"]<-"Fabaceae" #replace family name
fam<-fam[which(fam$species_matched %in% rownames(result)),] #subset only species included in our treatment

#get table with ranked families based on their number of species:
famf<-as.data.frame(table(fam$family)) #create dataframe
famf<-famf[order(-famf$Freq),] #order it
row.names(famf) <- NULL #get rid of rownames

#create vector with unique names of species:
list.r <- unique(fam$family)

#get nodes for families:
list.nod<-NULL
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$species_matched) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- findMRCA(tree, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

#clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)




###
# Plot rectangular phylogenetic tree
###

# Get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=50)) #select top families based on their number of taxa

# Plot tree 
p <- ggtree(tree, layout="rectangular", size=0.5, branch.letngth="none")+ # build circular tree
  geom_cladelabel(node=subset(famf, Var1==toplot[1])$num, label=toplot[1],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[2])$num, label=toplot[2],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[3])$num, label=toplot[3],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[4])$num, label=toplot[4],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[5])$num, label=toplot[5],  fontsize=4, barsize = 0.5, angle=360) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[6])$num, label=toplot[6],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[7])$num, label=toplot[7],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[8])$num, label=toplot[8],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[9])$num, label=toplot[9],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[10])$num, label=toplot[10],  fontsize=4, barsize = 0.5, angle=360) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[11])$num, label=toplot[11],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[12])$num, label=toplot[12],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[13])$num, label=toplot[13],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[14])$num, label=toplot[14],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[15])$num, label=toplot[15],  fontsize=4, barsize = 0.5, angle=360) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[16])$num, label=toplot[16],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[17])$num, label=toplot[17],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[18])$num, label=toplot[18],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[19])$num, label=toplot[19],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[20])$num, label=toplot[20],  fontsize=4, barsize = 0.5, angle=360) +

  geom_cladelabel(node=subset(famf, Var1==toplot[21])$num, label=toplot[21],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[22])$num, label=toplot[22],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[23])$num, label=toplot[23],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[24])$num, label=toplot[24],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[25])$num, label=toplot[25],  fontsize=4, barsize = 0.5, angle=360) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[26])$num, label=toplot[26],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[27])$num, label=toplot[27],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[28])$num, label=toplot[28],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[29])$num, label=toplot[29],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[30])$num, label=toplot[30],  fontsize=4, barsize = 0.5, angle=360) +

  geom_cladelabel(node=subset(famf, Var1==toplot[31])$num, label=toplot[31],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[32])$num, label=toplot[32],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[33])$num, label=toplot[33],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[34])$num, label=toplot[34],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[35])$num, label=toplot[35],  fontsize=4, barsize = 0.5, angle=360) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[36])$num, label=toplot[36],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[37])$num, label=toplot[37],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[38])$num, label=toplot[38],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[39])$num, label=toplot[39],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[40])$num, label=toplot[40],  fontsize=4, barsize = 0.5, angle=360) +

  #geom_cladelabel(node=subset(famf, Var1==toplot[41])$num, label=toplot[41],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[42])$num, label=toplot[42],  fontsize=4, barsize = 0.5, angle=360) +
  #geom_cladelabel(node=subset(famf, Var1==toplot[43])$num, label=toplot[43],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[44])$num, label=toplot[44],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[45])$num, label=toplot[45],  fontsize=4, barsize = 0.5, angle=360) +
  
  #geom_cladelabel(node=subset(famf, Var1==toplot[46])$num, label=toplot[46],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[47])$num, label=toplot[47],  fontsize=4, barsize = 0.5, angle=360) +
  #geom_cladelabel(node=subset(famf, Var1==toplot[48])$num, label=toplot[48],  fontsize=4, barsize = 0.5, angle=360) +
  #geom_cladelabel(node=subset(famf, Var1==toplot[49])$num, label=toplot[49],  fontsize=4, barsize = 0.5, angle=360) +
  geom_cladelabel(node=subset(famf, Var1==toplot[50])$num, label=toplot[50],  fontsize=4, barsize = 0.5, angle=360) +
  
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        #legend.position="bottom"
        legend.key.size = unit(1, "cm"),
        )

# Add heatmap for treatment:
p <- gheatmap(p, result, offset=86, width=0.7, 
              colnames = T, 
              colnames_position = "top",
              colnames_offset_y = 10,
              color = NULL,
              font.size=5,
              hjust = 0,
              #trace = "none",
              colnames_angle=45) +
  scale_fill_manual(labels=c("Lower DCi", "Higher DCi", "Random DCi", "Missing"),
                    values=c("red", "blue", "grey"), na.value="white",
                    na.translate=FALSE)+
  theme(legend.position="bottom") +
  ylim(NA, 2150)+
  xlim(NA, 800)


#save output:
png("phylo_ring_rect.png",
    res=300,height=20,width=12,units="in"); 
p
dev.off()


#clean-up:
#rm(list = ls())
