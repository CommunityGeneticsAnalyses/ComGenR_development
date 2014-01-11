cscore <-
function(x,cu.mat=FALSE){
  x[x!=0] <- 1 #force binary
  cu <- matrix(0,nrow=ncol(x),ncol=ncol(x))
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      ri <- sum(x[,i])
      rj <- sum(x[,j])
      S <- x[,i]*0
      S[x[,i]==1&x[,j]==1] <- 1
      S <- sum(S)
      cu[i,j] <- (ri-S)*(rj-S)
    }
  }
  if (cu.mat){return(cu)}else{return(mean(cu))}
}
