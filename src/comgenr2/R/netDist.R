netDist <-
function(dn.t){
  net.d <- matrix(0,nrow=length(dn.t),ncol=length(dn.t))
  rownames(net.d) <- colnames(net.d) <- names(dn.t)
  for (i in 1:nrow(net.d)){
    for (j in 1:ncol(net.d)){
      net.d[i,j] <- sum(abs(dn.t[[i]]-dn.t[[j]])^2)
    }
  }
  net.d <- as.dist(net.d)
  return(net.d)
}
