### CODED BY Matthew K. Lau
### 12Mar2012, updated 21May2012
### Method from Araujo et al. 2011
### M.B. Arau ́jo, A. Rozenfeld, C. Rahbek,
### and P.A. Marquet. Using species co- 
### occurrence networks to assess the 
### impacts of climate change. Ecography, 
### 34:897– 908, 2011.

### Function produces a co-occurrence network of dissimilarities
### using the Araujo et al. method.

araujoNet <- function(x,method='bray',min.abundance=1){
                                        #get the absent species names
  y <- colnames(x)[apply(x,2,sum) >= min.abundance]
                                        #remove low abundance species
  x <- x[,apply(x,2,sum) >= min.abundance]
                                        #assure presence-absence matrix
  x[x!=0] <- 1
                                        #warn if matrix is empty
  if (ncol(x) <= 1){warning('Community matrix is empty',quote=FALSE)}else{}
                                        #caluclate the Bray-Curtis dis.
  d <- as.matrix(vegdist(t(x),method=method))
                                        #calculate individual prob.
  Pa <- apply(x,2,function(x) length(x[x!=0])/length(x))
                                        #calculate null joint probability 
  Pab <- array(NA,dim=c(length(Pa),length(Pa)))
  rownames(Pab) <- colnames(Pab) <- names(Pa)
  for (i in 1:nrow(Pab)){
    for (j in 1:ncol(Pab)){
      Pab[i,j] <- Pa[i]*Pa[j]
    }
  }
                                        #Calculate confidence limits
                                        #number of co-occurrences
  N1 <- nrow(x) * Pab
                                        #calculate the variance
  V1 <- nrow(x) * Pab * (1-Pab)
                                        #+/- 2 SD limits (~95% confidence)
  cl.u <- N1 + 2*sqrt(V1)
  cl.l <- N1 - 2*sqrt(V1)
                                        #observed number of co-occurrences
  Nab <- N1 * 0 
  for (i in 1:nrow(Nab)){
    for (j in 1:ncol(Nab)){
      Nab[i,j] <- length(x[x[,i] != 0 & x[,j] != 0,i])
    }
  }
                                        #prune within confidence limits
  dp <- d * 0
  dp[Nab > cl.u] <- d[Nab > cl.u]
  dp[Nab < cl.l] <- d[Nab < cl.l]
                                        #pack for export
  out <- list(x=x,Pab=Pab,cl.l=cl.l,cl.u=cl.u,d=d,dp=dp)
  return(out)
}

### Example. Exploring the co-occurrence patterns in the 
### dune dataset from the vegan package.

## library(vegan)
## library(sna)
## data(dune)
## test <- araujoNet(dune)
## deg <- degree(test$dp[apply(test$dp,1,sum)!=0,apply(test$dp,2,sum)!=0]
##               ,rescale=FALSE)
## gplot(test$dp[apply(test$dp,1,sum)!=0,apply(test$dp,2,sum)!=0],
##       displaylabels=TRUE,gmode='graph',label.cex=0.65,
##       vertex.sides=50,vertex.col='black',edge.col='darkgrey',
##       vertex.cex=deg,edge.lwd=(test$dp[apply(test$dp,1,sum)!=0,
##         apply(test$dp,2,sum)!=0]+1)^2)

library(rjags)

bcJP <- function(x,y,file='',nits=100,n.adapt=10,n.chains=2,diagnostics=FALSE){
                                        #repack data
  x <- cbind(x,y)
                                        #calculate co-occurrence
                                        #and number of observations
  x <- apply(x,1,sum)
  x[x!=2] <- 0
  x[x==2] <- 1
  n <- length(x)
  x <- sum(x)
                                        #create data file
  if (file == ''){file = 'bcn'}else{}
  write.table(x,file = paste(file,'.data',sep=''),row.names = FALSE,col.names = FALSE)
                                        #create model
  write.table('model{',file=paste(file,'.bug',sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table('     x~dbin(p,n)',file=paste(file,'.bug',sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
  write.table('     p~dunif(0,1)',file=paste(file,'.bug',sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
  write.table('}',file=paste(file,'.bug',sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
                                        #initialize model
  
  jags <- jags.model('bcn.bug',data = list('x' = x,'n' = n),n.chains = n.chains,n.adapt = n.adapt,quiet=TRUE)
                                        #update model
  update(jags, nits,quite=TRUE)
                                        #sample from posterior
  samples <- coda.samples(jags,c('p'),nits,quiet=TRUE)
  if (diagnostics==TRUE){plot(samples)}else{}
                                        #pack for export
  return(c(summary(samples)$statistics[1],summary(samples)$quantiles[c(1,5)]))
}


bcNet <- function(x,method='bray',min.abundance=3){
                                        #get the absent species names
  y <- colnames(x)[apply(x,2,sum) >= min.abundance]
                                        #remove low abundance species
  x <- x[,apply(x,2,sum) >= min.abundance]
                                        #assure presence-absence matrix
  x[x!=0] <- 1
                                        #warn if matrix is empty
  if (ncol(x) <= 1){warning('Community matrix is empty',quote=FALSE)}else{}
                                        #caluclate the Bray-Curtis dis.
  d <- as.matrix(vegdist(t(x),method=method))
                                        #calculate individual prob.
  Pa <- apply(x,2,function(x) length(x[x!=0])/length(x))
                                        #calculate null joint probability 
  Pab <- array(NA,dim=c(length(Pa),length(Pa)))
  rownames(Pab) <- colnames(Pab) <- names(Pa)
  for (i in 1:nrow(Pab)){
    for (j in 1:ncol(Pab)){
      if (i >= j){Pab[i,j] <- Pa[i]*Pa[j]}else{}
    }
  }
  Pab[upper.tri(Pab)] <- t(Pab)[upper.tri(Pab)]
                                        #Calculate confidence limits
                                        #95% Posterior Credible Interval
  ci.l <- Pab * 0
  ci.u <- Pab * 0
  postMu <- Pab * 0
  count <- 0
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      count <- count + 1
      if (i >= j){
        postJP <- bcJP(x[,i],x[,j])
        ci.l[i,j] <- postJP[2]
        ci.u[i,j] <- postJP[3]
        postMu[i,j] <- postJP[1]
        print(paste(round((count/(ncol(x)^2)*100),0),'% complete',sep=''),quote=FALSE)
      }else{}
    }
  }
                                        #fill upper matrix (matrices are summetric)
  ci.l[upper.tri(ci.l)] <- t(ci.l)[upper.tri(ci.l)]
  ci.u[upper.tri(ci.u)] <- t(ci.u)[upper.tri(ci.u)]
  postMu[upper.tri(postMu)] <- t(postMu)[upper.tri(postMu)]
                                        #prune within confidence limits
                                        #distances/dissimilarities
  dp <- d * 0
  dp[Pab > ci.u] <- d[Pab > ci.u]
  dp[Pab < ci.l] <- d[Pab < ci.l]
                                        #posterior joint probabilities
  postMup <- postMu * 0
  postMup[Pab > ci.u] <- postMu[Pab > ci.u]
  postMup[Pab < ci.l] <- postMu[Pab < ci.l]
                                        #pack for export
  out <- list(x=x,d=d,dp=dp,postJP=postMu,postJPp=postMup)
  return(out)
}

