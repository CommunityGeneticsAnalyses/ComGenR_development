gpmTrees <-
function(pheno='phenotype of genotypes',reps='replication'){
  if (any(pheno=='phenotype of genotypes')){
    pheno <- c(11.00,12.50,13.75,16.00,14.00,15.25,17.50,16.50,18.75,21.00)
  }
  if (reps=='replication'){reps <- 5}
  trees <- list()
  for (i in 1:length(pheno)){
    trees[[i]] <- rep(pheno[i],reps)
  }
  return(cbind(geno=gl(length(pheno),reps),pheno=unlist(trees)))
}
