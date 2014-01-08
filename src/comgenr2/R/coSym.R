coSym <-
function(x='dependency network',zero.na=TRUE){
  out <- x * 0
  for (i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      if (max(c(x[i,j],x[j,i]))==0){}else{
        out[i,j] <- abs(x[i,j]-x[j,i])/max(c(x[i,j],x[j,i]))
      }
    }
  }
  return(out)
}
