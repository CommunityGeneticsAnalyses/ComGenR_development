
net.dif <- function(dat,dat2=NULL,g1=NULL,g2=NULL){
  sum((dat[,,g1]-dat[,,g2])^2)
}

reduce.net <- function(x){
  if (length(x)!=2){warning('Error: matrix list not of length 2!')}
  x1 <- x[[1]]
  x2 <- x[[2]]
  if ((dim(x1)[1]!=dim(x1)[2])|(dim(x2)[1]!=dim(x2)[2])){warning('Matrices are not square!')}
  if (any(dim(x1)!=dim(x2))){warning('Dimensions are not equal!')}
  if (any(rownames(x1)!=rownames(x2))){warning('Node names are mismatched!')}
  out <- list()
  out[[1]] <- x1[(apply(x1,1,sum)!=0|apply(x2,1,sum)!=0),(apply(x1,2,sum)!=0|apply(x2,2,sum)!=0)]
  out[[2]] <- x2[(apply(x1,1,sum)!=0|apply(x2,1,sum)!=0),(apply(x1,2,sum)!=0|apply(x2,2,sum)!=0)]
  names(out) <- names(x)
  return(x)
}

cgQAP <- function(x,FUN=net.dif,nits=999){
  x <- reduce.net(x)  
  q <- array(0,dim=c(dim(x[[1]]),2))
  q[,,1] <- x[[1]]
  q[,,2] <- x[[2]]
  out <- qaptest(q,FUN,g1=1,g2=2,reps=nits)
  return(out)
}

