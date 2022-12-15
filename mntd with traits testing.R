spcomp_wide <- read.csv("kruger_comp_ex.csv")
traits_sp_means <- read.csv("kruger_traits_ex.csv")

trait_dist <- dist(traits_sp_means, "euclidean", diag=T, upper=T)
comp_matrix <- as.matrix(spcomp_wide[,c("ARICON", "BOTRAD","DIGERI","ERASUP","PANCOL","PANMAX","SCHPAP","SETINC","THETRI","UROMOS")])

mntd(samp=comp_matrix, dis=trait_dist)

