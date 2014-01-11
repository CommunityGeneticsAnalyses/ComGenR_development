dep.net <-
function(x='species in cols',zero.na=TRUE,prune=TRUE,diag.zero=TRUE,pos=TRUE){
  out <- matrix(NA,nrow=ncol(x),ncol=ncol(x))
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      if (pos){
        out[i,j] <- calcDepend(x[,i],x[,j])
      }else{
        out[i,j] <- negDepend(x[,i],x[,j])
      }
    }
  }
  if (prune){
    out.rm <- co.net(x,diag.zero=diag.zero)
    out[out.rm==0] <- 0
  }else{}
  if (diag.zero){diag(out) <- 0}
  rownames(out) <- colnames(out) <- colnames(x)
  if (zero.na){out[is.na(out)] <- 0}
  return(out)
}
