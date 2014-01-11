co.net <-
function(x='species in cols',diag.zero=TRUE,std=TRUE){
  x[x!=0] <- 1
  out <- matrix(NA,nrow=ncol(x),ncol=ncol(x))
  rownames(out) <- colnames(out) <- colnames(x)
  for (j in 1:ncol(x)){
    for (k in 1:ncol(x)){
      out[j,k] <- null.prune(x[,j],x[,k],std=std)
    }
  }
  if (diag.zero){diag(out) <- 0}else{}
  return(out)
}
