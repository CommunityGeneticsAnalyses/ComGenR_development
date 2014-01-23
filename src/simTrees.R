simTrees <- function(tree.gpm='tree genotype-phenotype map',VeT=2){
  if (length(tree.gpm)==1){tree.gpm <- gpmTrees()}
  T <- nrow(tree.gpm) #number of trees
  trees <- tree.gpm[,2] + runif(T,0,1) * VeT - VeT/2 #scores are in second column of trees
  return(trees)
}
