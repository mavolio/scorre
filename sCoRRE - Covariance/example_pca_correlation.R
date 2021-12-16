## get data and run PCA
require(vegan)
# standardise to zero mean and unit sd before analysis
D1 = apply(iris[1:75,-5],2,function(x) (x-mean(x))/sd(x))
D2 = apply(iris[76:150,-5],2,function(x) (x-mean(x))/sd(x))

iris.pca1 <- rda(D1)
iris.pca2 <- rda(D2)

## pca1
iris.pca1$CA$v # the PCA vectors - these are the formula for relationships among traits

# overall explanatory power of each axis
iris.pca1$CA$eig/sum(iris.pca1$CA$eig) # explanatory power of axes based on eigenvalues

get_pca_exp = function(pca, data, dimensions, observations = 1:nrow(data)) {
  # r2-like index for PCA
  cor(c(pca$CA$u[observations,dimensions,drop=FALSE] %*% t(pca$CA$v[,dimensions,drop=FALSE])), c(data[observations,]))^2
}

# run across all data in D1
get_pca_exp(iris.pca1, D1, 1) # dimension 1
get_pca_exp(iris.pca1, D1, 2) # dimension 2
get_pca_exp(iris.pca1, D1, 3) # dimension 3
get_pca_exp(iris.pca1, D1, 4) # dimension 4
# matches explanatory power from eigenvalues

# run on subsets of data (just first 10 observations)
get_pca_exp(iris.pca1, D1, 1, observations = 1:10) # dimension 1
get_pca_exp(iris.pca1, D1, 2, observations = 1:10) # dimension 2
get_pca_exp(iris.pca1, D1, 3, observations = 1:10) # dimension 3
get_pca_exp(iris.pca1, D1, 4, observations = 1:10) # dimension 4


## calculate correspondence between pca1 and pca2
require(MatrixCorrelation)
PSI(iris.pca1$CA$v, iris.pca2$CA$v, center = TRUE) # very high correlation
