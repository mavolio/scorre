
## funtion to compute mean trait values for each clade in the phylogeny and
## compare it to a random expectation in which observed values are shuffled 
## across the phylogeny.
## 
## It returns 6 columns:
## Node: the number of the node in the phylogeny.
## SR: the number of taxa in each node.
## Obs: observed mean trait value for that clade
## Mean_Exp: mean expected mean trait value after randomizations
## SD_Exp: standard deviation of mean trait value after randomizations
## P_value: associated P-value
##
## Input data requires:
## tree: an object of clas "phylo".
## samp: a dataframe object with 1 single column with the trait. Rownames are species names.
##       The nomenclature between "tree" and "samp" must match.
## N: number of randomizations.

library(geiger)

node.mean <- function(tree, samp, N) {
  n.internal.nodes <- tree$Nnode # number of nodes in tree
  n.tips <- length(tree$tip.label) # number of tips in tree
  
  ## prepare dataframe to store results:
  out<-as.data.frame(matrix(ncol = N+1, nrow = n.internal.nodes)) #create object to store results
  rownames(out)<-seq((n.tips+1), (n.tips+n.internal.nodes)) #rownames are the names (numbers) of internal nodes

  ## calculating mean observed value for each node:
  sr.node<-list() #create object to store the number of species  found in each node.
  for(i in (n.tips+1):(n.tips+n.internal.nodes)){
    node <- i
    spp.node <- tips(tree, node) # get taxa in node
    out[i-n.tips,1] <- mean(samp[rownames(samp) %in% spp.node, ], na.rm = T) # subset taxa and get the mean
    sr.node[[i]]<-length(samp[rownames(samp) %in% spp.node, ])
  }
  
  ## calculating mean expected values for each node:
  for(j in 1:N){
    tree2<-tree #create copy of the tree
    set.seed(123+j) #for reproducibility
    tree2$tip.label<-sample(tree2$tip.label) #shuffle tips in the phylogeny
    
    for(i in (n.tips+1):(n.tips+n.internal.nodes)){
      node <- i
      spp.node <- tips(tree2, node) # get taxa in node
      out[i-n.tips, j+1] <- mean(samp[rownames(samp) %in% spp.node, ], na.rm = T) # subset taxa and get the mean
    }
  }
  
  ## prepare data frame to store results:
  out2<-as.data.frame(matrix(ncol = 6, nrow = n.internal.nodes)) #create object to store results
  colnames(out2)<-c("Node", "SR", "Obs", "Mean_Exp", "SD_Exp", "P_value") #set names of columns
  out2$Node<-seq((n.tips+1), (n.tips+n.internal.nodes)) #rownames are the names (numbers) of internal nodes
  out2$SR<-unlist(sr.node)

  ## calculate P-value:  
  for(j in 1:nrow(out)){
    mean.h <- mean(as.numeric(out[j, 2:(N+1)]))
    sd.h <- sd(as.numeric(out[j, 2:(N+1)]))
    expect.s <- as.numeric(out[j, 2:(N+1)]) - mean.h
    obs.s <- out[j,1] - mean.h
    
    #store values:
    out2$Obs[j]<-out[j,1]
    out2$Mean_Exp[j]<-mean.h
    out2$SD_Exp[j]<-sd.h
    out2$P_value[j] <- mean(abs(expect.s)>abs(obs.s))
  }
  return(out2)
}
