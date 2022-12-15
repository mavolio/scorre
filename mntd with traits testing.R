#spcomp_wide <- read.csv("kruger_comp_ex.csv")
spcomp_wide<-read.csv("kruger_comp_ex.csv")
comp_matrix <- as.matrix(spcomp_wide[,c("ARICON", "BOTRAD","DIGERI","ERASUP","PANCOL","PANMAX","SCHPAP","SETINC","THETRI","UROMOS")])

#traits_sp_means <- read.csv("kruger_traits_ex.csv")
traits_sp_means<-read.csv("kruger_traits_ex.csv")
rownames(traits_sp_means)<-colnames(comp_matrix)
trait_dist <- dist(traits_sp_means, "euclidean", diag=T, upper=T)

picante::mntd(comp_matrix, as.matrix(trait_dist))


