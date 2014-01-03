## Dependencies
library(pbapply)
library(audio)
library(bipartite)
library(gam)
library(glmnet)

###my gplots
mgp <- function(x,vp=2,vcol='black',empty.v=TRUE,ep=2,escale=2){
  x <- abs(x)
  gam.attach <- any(search()=='package:gam')
  if (gam.attach){detach(package:gam)}
  if (any(search()=='package:sna')==FALSE){require('sna',quiet=TRUE)}
  v.size <- ((abs(degree(x))/max(abs(degree(x))))+1)^vp
  if (empty.v){
    vcol <- rep(vcol,nrow(x))
    vcol[apply(x,1,function(x) sum(abs(x)))==FALSE] <- 'white'
  }else{}
  e.size <- ((abs(x)/max(abs(x)))^ep)*escale
  gplot(abs(x),displaylabels=FALSE,mode='circle',gmode='graph',
        vertex.cex=v.size,vertex.sides=25,vertex.col=vcol,
        edge.lwd=e.size,vertex.border='grey',edge.col='grey')
  if (gam.attach){require('gam',quiet=TRUE)}
}


###GAMS
model.gam <- function(x,y,prob=0.95,aic=TRUE,coef=TRUE){
  test.null <- gam(y~1)
  test. <- gam(y~x)
  test.. <- gam(y~x+I(x^2))
  tests <- list(test.null,test.,test..)
                                        #model probability
  pm <- unlist(lapply(tests,AIC))
  pm <- exp((min(pm)-pm)/2)
  pm[is.na(pm)] <- 0
                                        #decide which model is optimal
  if (all(pm<prob)){m.out <- 0}else{
    m.out <- ((1:3)[pm==max(pm)])[1]
    if (m.out==1){
      m.out <- 0
    }else if (m.out==2){
      m.out <- 1
    }else if (m.out==3){
      m.out <- 2
    }
  }
                                        #return aic or coefficient values
  if (aic){
    if (m.out==1){
      if (coef){m.out <- coefficients(tests[[m.out]])[2]}else{m.out <- abs(pm[1]-pm[2])}
    }else if (m.out==2){
      if (coef){m.out <- coefficients(tests[[m.out]])[2]}else{m.out <- abs(pm[1]-pm[3])}
    }
  }else{}
  if (is.na(m.out)){m.out <- 0}
                                        #finish
  return(m.out)
}

gamNet <- function(x,prob=0.95,aic=TRUE,coef=TRUE){
  out <- array(NA,dim=rep(ncol(x),2))
  count <- 1
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      if (i<=j){}else{
        out[i,j] <- model.gam(x[,i],x[,j],prob=prob,aic=aic,coef=coef)
        out[j,i] <- out[i,j]
        print(paste(round((count/(ncol(x)^2)*100),0),'% complete',sep=''),quote=FALSE)
      }
      count <- count + 1
    }
  }
  rownames(out) <- colnames(out) <- colnames(x)
  diag(out) <- 0
  return(out)
}

###LASSO
print('lasso')

##Spatial Co-occurrence networks (Araujo method)
##Ara\'ujo, M.B., Rozenfeld A., Rahbek C. and Marquet P.A. 2011 Using species co-occurrence networks to assess the impacts of climate change. Ecography 000:000-000.

##Building spaital networks of co-occurrence
#1. measure the number of sites where species co-occur
#calculate athe Bray-Curtis dissimilarity for all species

#2. generate the network based on the distributional overlap between species distributions

#3. prune the network of spurious distributions
                                        #notes from paper
## A <- #site adjacency matrix (e.g. lat X lon grid)
## N <- nrow(x) #total number of sites
## Spa <- x[,1] #species a observations
## Na <- length(Spa[Spa != 0]) #number of sites containing A
## Pa <- Na/N #probability of observing species A

                                        #Pruned network function
library(vegan)
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

pruneNetwork <- function(x,method='bray'){
  y <- colnames(x)[apply(x,2,sum) == 0]
  x <- x[,apply(x,2,sum) != 0]
  if (ncol(x) <= 1){break}else{}
  d <- as.matrix(vegdist(t(x),method=method))
                                        #calculate independent probabilities
Pa <- apply(x,2,function(x) length(x[x!=0])/length(x))
                                        #calculate joint probability under independence
Pab <- array(NA,dim=c(length(Pa),length(Pa)))
rownames(Pab) <- colnames(Pab) <- names(Pa)
for (i in 1:nrow(Pab)){
  for (j in 1:ncol(Pab)){
    Pab[i,j] <- Pa[i]*Pa[j]
  }
}
diag(Pab) <- Pa #set the diagonal to the indpendent probabilty of each species
                                        #probability of either species alone under independence
Paorb <- Pab*0 
for (i in 1:nrow(Paorb)){
  for (j in 1:ncol(Paorb)){
    Paorb[i,j] <- Pa[i]+Pa[j]-(2*Pa[i]*Pa[j])
}
}
diag(Paorb) <- 0 #set the diagonal to zero
                                        #probability of neither species under independence
Pnotab <- Paorb*0
for (i in 1:nrow(Pnotab)){
  for (j in 1:ncol(Pnotab)){
    Pnotab[i,j] <- 1-Pa[i]-Pa[j]+Pa[i]*Pa[j]
}
}
diag(Pnotab) <- 1-Pa #set the diagonal to the complement of Pa
                                        #check, Pab + Paorb + Pnotab = 0
if (all(round(Pab + Paorb + Pnotab,7) == 1)){}else{break}
                                        #Calculate confidence limits
                                        #number of sites in each state:
                                        #1. containing both species a and b
N1 <- nrow(x) * Pab
                                        #2. containing one of species a or b
                                        #N2 <- nrow(x) * Paorb
                                        #3. neither species a or b
                                        #N3 <- nrow(x) * Pnotab
                                        #Calculate the variances of each state
V1 <- nrow(x) * Pab * (1-Pab)
                                        #V2 <- nrow(x) * Paorb * (1-Paorb) 
                                        #V3 <- nrow(x) * Pnotab * (1-Pnotab)
                                        #+/- 2 SD limits (~95% confidence)
cl.u <- N1 + 2*sqrt(V1)
cl.l <- N1 - 2*sqrt(V1)
                                        #pruning
Nab <- N1 * 0 #observed number of co-occurrences
for (i in 1:nrow(Nab)){
  for (j in 1:ncol(Nab)){
    Nab[i,j] <- length(x[x[,i] != 0 & x[,j] != 0,i])
  }
}
                                        #prune within confidence limits
dp <- d * 0
dp[Nab > cl.u] <- d[Nab > cl.u]
dp[Nab < cl.l] <- d[Nab < cl.l]
                                        #add back absent species
x.names <- c(colnames(x),y)
x <- cbind(x,array(0,dim=c(nrow(x),length(y))))
colnames(x) <- x.names
d <- cbind(d,array(0,dim=c(nrow(d),length(y))))
d <- rbind(d,array(0,dim=c(length(y),ncol(d))))
dp <- cbind(dp,array(0,dim=c(nrow(dp),length(y))))
dp <- rbind(dp,array(0,dim=c(length(y),ncol(dp))))
out <- list(x=x,d=d,dp=dp)
return(out)
}

#4. use a percolation threshold to analyze and compare network properties
## namely used to look for the point at which the network breaks apart into isolated clusters
##?????

##Mutual information network

#Joint probability
#Given a matrix of presence absence observations for a set of species, determine the observed joint probabilities of all pairs of species

jpMat <- function(x,independent=FALSE){
  x[x != 0] <- 1
  y <- array(0,dim=c(ncol(x),ncol(x)))
  rownames(y) <- colnames(y) <- colnames(x)
  if (independent == FALSE){
    for (i in 1:nrow(y)){
      for (j in 1:ncol(y)){
        y[i,j] <- length((x[,i] + x[,j])[(x[,i] + x[,j]) == 2]) / nrow(x)
      }
    }
  }else{
    for (i in 1:nrow(y)){
      for (j in 1:ncol(y)){
        y[i,j] <- (sum(x[,i])/length(x[,i])) * (sum(x[,j])/length(x[,j])) 
      }
    }
  }
  return(y)
}

##Conditional probability
#Given a matrix of presence absence observations for a set of species, determine the observed conditional probabilities of all pairs of species

cpMat <- function(x,na.zero=TRUE){
  x[x != 0] <- 1
  y <- array(0,dim=c(ncol(x),ncol(x)))
  rownames(y) <- colnames(y) <- colnames(x)
     for (i in 1:nrow(y)){
      for (j in 1:ncol(y)){
        y[i,j] <- length((x[,i] + x[,j])[(x[,i] + x[,j]) == 2]) / nrow(x)
        if (i >= j){y[i,j] <- y[i,j] / (sum(x[,i])/length(x[,i]))}else{y[i,j] <- y[i,j] / (sum(x[,j])/length(x[,j]))}
      }
    }
  if (na.zero){y[is.na(y)] <- 0}else{}
  return(y)
}

##overlapDegree: degree of overlap between pairs of species
overlapDegree <- function(x,zero.na=TRUE){
  y <- array(0,dim=c(ncol(x),ncol(x)))
  rownames(y) <- colnames(y) <- colnames(x)
  for (i in 1:nrow(y)){
    for (j in 1:ncol(y)){
      y[i,j] <- length(x[x[,i] == 1 & x[,j] == 1,i])/sum(x[,i])
    }
  }
  if (zero.na){y[is.na(y)] <- 0}
  return(y)
}

##symLink: a measure of the degree of link between two species
symLink <- function(x,zero.na=TRUE){
  y <- array(0,dim=c(ncol(x),ncol(x)))
  rownames(y) <- colnames(y) <- colnames(x)
  for (i in 1:nrow(y)){
    for (j in 1:ncol(y)){
      O.ab <- length(x[x[,i] == 1 & x[,j] == 1,i])/sum(x[,i])
      O.ba <- length(x[x[,i] == 1 & x[,j] == 1,i])/sum(x[,j])
      y[i,j] <- abs(O.ab-O.ba)/max(c(O.ab,O.ba))
    }
  }
  if (zero.na){y[is.na(y)] <- 0}
  return(y)
}

## Automated analysis of communities for co-occurrence patterns
co.nets <- function(x,min.dens=3,decode.names=TRUE){
                                        #separate data files
  y <- x[,-1]
  y <- y[,apply(y,2,sum)>=min.dens]
  if (decode.names){
    y.names <- colnames(y)
    colnames(y) <- 1:ncol(y)
  }
  x <- x[,1]
                                        #store output into "out"
  out <- list()
                                        #composition analysis
  y. <- apply(y,2,function(x) x/sum(x))
  y.[is.na(y.)] <- 0
  out[[1]] <- adonis(y.~x)
                                        #indicator species analysis
  library(labdsv)
  out[[2]] <- indval(y,x)
  detach(package:labdsv)
                                        #network models
  x. <- array((rep(x,length(unique(x)))),dim=c(length(x),length(unique(x))))
  for (i in 1:length(unique(x))){
    x.[x.[,i] == unique(x.[,i])[i],i] <- 1
    x.[x.[,i] != unique(x.[,i])[i],i] <- 0
    x.[,i] <- (x.[,i])
  }
  x. <- array(as.numeric(x.),dim=dim(x.))
  colnames(x.) <- unique(x)
  comat <- cbind(x.,y)
  comat[comat!=0] <- 1
  out[[3]] <- comat 
  if (decode.names){y.names <- c(colnames(comat)[1:length(unique(x))],y.names)}
                                        #pruned networks
  out[[4]] <- pruneNetwork(comat)$dp
                                        #conditional probability network
  out[[5]] <- cpMat(comat)
  out[[5]][out[[4]]==0] <- 0
                                        #define the components of out
  if (decode.names){
    out[[6]] <- y.names
    names(out) <- c('perm.anova','ind.spp','x','dist.net','cp.net','spp.names')
  }else{
    names(out) <- c('perm.anova','ind.spp','x','dist.net','cp.net')
  }
  return(out)
}

#Binary sum
bin.sum=function(x){
	x[x!=0]<-1
	sum(x)
	}


#Partial Correlation

pcor=function(x,method='pearson'){
	R=cor(x,method=method)
	p=1/(1-R^2)
	
	}
x=array(rnorm(100),c(10,10))

#produce a correlation network
cor.net=function(x,alpha,lim){
  
  lim.n=colnames(x)[apply(x,2,sum)<lim]
  x=x[,apply(x,2,sum)>=lim]
  col.n=append(colnames(x),lim.n)
  
  I=nrow(x)
  x=dchisq(eed(cor(x),I),1)
  diag(x)=0
  x[x<qchisq((1-alpha),1)]<-0
  
  x=cbind(x,array(0,c(nrow(x),length(lim.n))))
  x=rbind(x,array(0,c(length(lim.n),ncol(x))))
  rownames(x)=colnames(x)=col.n
  
  return(x)
  	
}

cor.net=function(x,alpha,lim,method='pearson',p.adj=TRUE){
 
 lim.n=colnames(x)[apply(x,2,sum)<lim]
 x=x[,apply(x,2,sum)>=lim]
 col.n=append(colnames(x),lim.n)
 
 I=nrow(x)
 y=cor(x,method=method)
 x=dchisq(eed(y,I),1)
 diag(x)=0
 
 if (p.adj==TRUE){
 	
 	}
 
 y=cbind(y,array(0,c(nrow(y),length(lim.n))))
 y=rbind(y,array(0,c(length(lim.n),ncol(y))))
 rownames(y)=colnames(y)=col.n
 
 return(y)
 	
 	}

#Correlation Pairwise Test
kendall.pairs=function(x,alpha=0.05,p.adj=TRUE,adj.method='holm',raw=FALSE){

  p=array(0,c(ncol(x),ncol(x)))
  colnames(p)=rownames(p)=colnames(x)
  r=p
  
  for (i in seq(along=rownames(p))){
    for (j in seq(along=colnames(p))){
      if (i==j|i<j){}
      else{
        if (sum(x[,i])==0|sum(x[,j])==0){
          p[i,j]=1
          r[i,j]=0
        }else{
          out=unlist(suppressWarnings(cor.test(x[,i],x[,j],method='kendall')))
          p[i,j]=round(as.numeric(out[2]),7)
          r[i,j]=round(as.numeric(out[3]),7)
        }
      }	
    }
  }
  
  if (p.adj==TRUE){
    p=p[lower.tri(p)]
    p=p.adjust(p,method=adj.method)
  }else{p=p[lower.tri(p)]}
  
  out=r*0
  out[lower.tri(out)]=p
  p=out
  p=as.matrix(as.dist(p))
  r=as.matrix(as.dist(r))
  r.=r
  r.[p>=alpha]=0	
  if (raw==TRUE){
    out=list(r,p)	
  }else{
    out=r.
  }
  
  return(out)
  
}

#Distance based calculation of correlations


my.solve=function (a, b, tol = 1e-07) 
{
    if (!is.qr(a)) 
        a <- qr(a, tol = tol)
    nc <- ncol(a$qr)
    nr <- nrow(a$qr)
    if (a$rank != min(nc, nr)){
    	out="singular matrix 'a' in solve"
    	} 
    if (missing(b)) {
        if (nc != nr) 
            out="only square matrices can be inverted"
        b <- diag(1, nc)
    }else{
    out <- qr.coef(a, b)
    out[is.na(out)] <- 0	
    	}
    out
}

dcor.net=function(x,alpha=0.05,by.lim=0.01,max.lim=1,rm.spp=TRUE){
si.=x
si.p='si.p'
si.c=0
lim=0

lim=lim+by.lim
si.=si.[,apply(si.,2,sum)>lim]
si.d=vegdist(t(si.))
si.d=as.matrix(si.d)
si.c=gower(si.d)
si.p=my.solve(si.c)

while(any(si.p=="singular matrix 'a' in solve")){
lim=lim+by.lim
si.=si.[,apply(si.,2,sum)>lim]
si.d=vegdist(t(si.))
si.d=as.matrix(si.d)
si.c=gower(si.d)
si.p=my.solve(si.c)
if (lim>=max.lim){stop('ERROR: Singular matrix')}
}

si.p=solve(si.c)

print(lim)
print(si.p)

si.r=si.p*0

for (i in seq(1:nrow(si.p))){
	for (j in seq(1:ncol(si.p))){
	si.r[i,j]=(-si.p[i,j])/sqrt((si.p[i,i]*si.p[j,j]))
		}
	}

si.R=si.r*-1
diag(si.R)=1

si.eed=eed(si.R,nrow(x))
si.R[si.eed<qchisq((1-alpha),1)]<-0

if (rm.spp==TRUE){
rm.spp=colnames(x)[apply(x,2,sum)<=lim]
si.R=cbind(si.R,array(0.0000000,c(nrow(si.R),length(rm.spp))))
si.R=rbind(si.R,array(0.0000000,c(length(rm.spp),ncol(si.R))))
colnames(si.R)=rownames(si.R)=c(colnames(si.R)[1:(ncol(si.R)-length(rm.spp))],rm.spp)
	}else{}

return(si.R)
	}
	
#Given a list of adjacency matrixes, return the average of all edge weights
edge.mean=function(x){
	y=x[[1]]*0
	for (i in (1:length(x))){
		y=y+x[[i]]
		}
	y=y/length(y)
	return(y)
	}

#Kullback's chi-square test statistic for the difference between two correlation matrixes

chi.cor=function(x1,x2,n1,n2){
n=c(n1,n2)
Nk=n-1
N=sum(Nk)
Rk=list(x1,x2)
R=Rk[[1]]*0
for (i in (1:length(Rk))){
	R=R+(Nk[i]*Rk[[i]])
	}
R=R/N
Chi2.k=0
for (i in (1:length(Rk))){
	Chi2.k=Chi.k+(Nk[i]*log((det(R)/det(Rk[[i]]))))
	}
Chi2.df=((length(Rk)-1)*ncol(R)*(ncol(R)-1))/2	
p.val=pchisq(Chi2.k,Chi2.df,lower.tail=FALSE)
return(c(Chi2.df,Chi2.k,p.val))
}


#PermCor: permutes two input matrixes from a list
PermCor=function(x){
	x1=x[[1]]
	x2=x[[2]]
	x.=c(x1[lower.tri(x1)],x2[lower.tri(x2)])
	perm=sample(1:length(x.))
	x1.=array(dim=dim(x1))
	x2.=array(dim=dim(x2))
	x1.[lower.tri(x1.)]=x.[perm[(1:(length(x.)/2))]]
	x2.[lower.tri(x2.)]=x.[perm[(((length(x.)/2)+1):length(x.))]]
	x1.=as.matrix(as.dist(x1.))
	x2.=as.matrix(as.dist(x2.))
	colnames(x1.)=colnames(x2.)=rownames(x1.)=rownames(x2.)=colnames(x1)
	perm=list(x1.,x2.)
	names(perm)=names(x)[1:2]
	return(perm)
}

#Shared number of connections between networks, all graphs assumed to be undirected
shared=function(x,diag=FALSE){
if (diag==FALSE){
	for (i in (1:length(x))){diag(x[[i]])=0}}
for (i in (1:length(x))){x[[i]][x[[i]]!=0]=1}
for (i in (2:length(x))){x[[1]]=x[[1]]+x[[i]]}
x=x[[1]]
return(table(x[upper.tri(x)]))
}

#rm.isolated

rm.isolated=function(x){
	y=x
	x[x!=0]=1
	y=y[apply(x,1,sum)>=1,]
	y=y[,apply(x,2,sum)>=1]
	return(y)
}

#fragmentation (Begotti's )
fragmentation=function(x){
	require('igraph',quietly=TRUE)
	N=nrow(x)
	x=graph.adjacency(x,weighted=TRUE,mode='undirected')
	Nj=clusters(x)$csize
	F=1-(sum(Nj*(Nj-1))/(N*(N-1)))
	detach(package:igraph)
	return(F)
}


##Network Distances
#Given a list of input networks of equal size returns a square matrix of the pairwise differences

getFlow <- function(x='network object list'){
                                        #retrieve the flow matrices along the length of x
  for (i in (1:length(x))){
    x[[i]] = x[[i]]%n%'flow'
  }

  return(x)

}

mDist <- function(M1,M2){
  sqrt(sum((M1 - M2)^2)) #compute the Euclidean distance between two matrices (square-root of the sum of squared differences between all pairs of cells)
}


netDist <- function(x='matrix list',dist.func=mDist){

  D = matrix(NA,length(x),length(x)) #initialize the output distance matrix

                                        #pairwise distance calculations
  for (i in (1:length(x))){
    for (j in (1:length(x))){
      D[i,j] = dist.func(x[[i]],x[[j]])
    }
  }  
                                        #row and column naming
  if (length(names(x)) == length(x)){rownames(D) = colnames(D) = names(x)}
  
  return(D)
}	

	
##Graph size
netSize <- function(x='graph matrix'){
x[x != 0] <- 1
x <- apply(x,1,sum)
length(x[x > 0])
}


##Graph degree
netDegree <- function(x='graph matrix'){
x[x != 0] <- 1
sum(x)
}


##Local network shifts
netShift <- function(x='graph matrix',y='graph matrix',order=TRUE){
x <- abs(x - y)
x <- apply(x,1,sum)
if (order == TRUE){x <- sort(x,decreasing=TRUE)}else{}
return(x)
}

##miniGraph: produces a minimal set of edges for two graphs having removed all nodes that have no connections in either graph.
miniGraph <- function(x){
  if (mode(x) != 'list' | length(x) != 2){
    warning('Error: x must be a list of length 2.')
  }else{
    x <- conform(x[[1]],x[[2]])
    y <- abs(x[[1]]) + abs(x[[2]])
    x[[1]] <- x[[1]][apply(y,1,bin.sum) != 0,apply(y,2,bin.sum) != 0]
    x[[2]] <- x[[2]][apply(y,1,bin.sum) != 0,apply(y,2,bin.sum) != 0]
    return(x)
  }
}

##getBiModules: returns the graph and the module memberships from a moduleWeb object
getBiModules <- function(mod,moduleweb=TRUE){
  if (moduleweb == TRUE){
    out <- list(rownames(slot(mod,'moduleWeb')),colnames(slot(mod,'moduleWeb')))
  }else{
    out <- list(rownames(slot(mod,'originalWeb')),colnames(slot(mod,'originalWeb')))
  }
  x <- listModuleInformation(mod)[[2]]
  for (i in 1:length(x)){
    for (j in 1:length(x[[i]][[1]])){
      out[[1]][out[[1]] == x[[i]][[1]][j]] <- i
    }
    for (j in 1:length(x[[i]][[2]])){
      out[[2]][out[[2]] == x[[i]][[2]][j]] <- i
    }
  }
  if (moduleweb == TRUE){
    out <- list(slot(mod,'moduleWeb'),out)
  }else{
    out <- list(slot(mod,'originalWeb'),out)
  }
  return(out)
}

##Co-occurrence Analysis
Cscore.pairs <- function(x){
  x[x!=0] <- 1
  out <- matrix(0,nrow=ncol(x),ncol=ncol(x))
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      out[i,j] <- C.score(x[,c(i,j)],normalize=FALSE)
    }
  }
  rownames(out) <- colnames(out) <- colnames(x)
  return(out)
}
                                        #random co-occurrence patterns
nullSim <- function(x,mc=1000,method='r1',alert=FALSE){
  x[x!=0] <- 1
  null <- list()
  for (i in 1:mc){
    null[[i]] <- commsimulator(x,method=method)
    print(paste(round(i/mc*100,0),'%',' Complete',sep=''))
  }
  if (alert == TRUE){play(c(sin(1:10000/20),rep(0,1000),sin(1:10000/10)))}else{}
  return(null)
}

                                        #standardized co-occurrence analysis output
CA.results <- function(x,nmc=5000,sim.out=TRUE,zero.na=TRUE){
  sim <- nullSim(x,mc=nmc)
  cs <- unlist(pblapply(sim,C.score))
  if (zero.na){cs[is.na(cs)] <- 0}else{}
  SES <- (C.score(x)-mean(cs))/sd(cs)
  if (SES>=0){p.val <- length(cs[cs>=C.score(x)])}else{p.val <- length(cs[cs<=C.score(x)])}
  cat <- c(obs=C.score(x),mu=mean(cs),
           sd=sd(cs),SES=SES,
           p.val=p.val/length(cs)) 
  cat <- list(sim=cs,cat=cat)
  return(cat)
}

##Conform matrices
conform <- function(x,y){
  mismatch <- list(colnames(y)[is.na(match(colnames(y),colnames(x)))],
                   colnames(x)[is.na(match(colnames(x),colnames(y)))])
                                        #adjust x
  add.x <- array(0,dim=c(length(mismatch[[1]]),ncol(x)))
  rownames(add.x) <- mismatch[[1]]
  x <- rbind(x,add.x)
  add.x <- array(0,dim=c(nrow(x),length(mismatch[[1]])))
  colnames(add.x) <- mismatch[[1]]
  x <- cbind(x,add.x)
                                        #adjust y
  add.y <- array(0,dim=c(length(mismatch[[2]]),ncol(y)))
  rownames(add.y) <- mismatch[[2]]
  y <- rbind(y,add.y)
  add.y <- array(0,dim=c(nrow(y),length(mismatch[[2]])))
  colnames(add.y) <- mismatch[[2]]
  y <- cbind(y,add.y)
                                        #order
  x <- x[order(colnames(x)),order(colnames(x))]
  y <- y[order(colnames(y)),order(colnames(y))]
                                        #
  return(list(x=x,y=y))
}
