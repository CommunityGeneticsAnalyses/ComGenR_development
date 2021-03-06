###enaInformation
##Calculates the information content of an ecosystem model or community dataset. 
##For community datasets, an optional graph of significant co-occurrence patterns is returned for all species pairs
##The information content calculations are from Ullanowicz 2011 Information
##The co-occurrence Arau'jo et al. 2011 Ecography

                                        #probability of x
indProb <- function(x){
  length(x[x == 1])/length(x)
}

                                        #probability of x and y
jointProb <- function(x,y){
  length((x+y)[(x+y)==2])/length((x+y))
}

                                        #conditional probability
condProb <- function(x,y){
  jointProb(x,y) / indProb(y)
}

                                        #probability of x not y
exclProb <- function(x,y){
  indProb(x) - jointProb(x,y)
}

                                        
probNet <- function(x,d=2,sig.test=TRUE){
                                        #convert to presence absence
  x[x!=0] <- 1
                                        #calculate probabilities
  Pa <- apply(x,2,sum) / nrow(x) #independent probability vector
  Pab <- array(NA,dim=c(ncol(x),ncol(x))) #joint probability matrix
  Paorb <- array(NA,dim=c(ncol(x),ncol(x))) #directional non-independent probability matrix
  Pagb <-  array(NA,dim=c(ncol(x),ncol(x))) #condiational probability matrix
  Nab <-  array(NA,dim=c(ncol(x),ncol(x))) #number of co-occurrences
                                        #name rows and columns
  rownames(Pab) <- colnames(Pab) <- colnames(x)
  rownames(Paorb) <- colnames(Paorb) <- colnames(x)
  rownames(Pagb) <- colnames(Pagb) <- colnames(x)
  rownames(Nab) <- colnames(Nab) <- colnames(x)
                                        #loop for calculations
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      Pab[i,j] <- jointProb(x[,i],x[,j])
      Paorb[i,j] <- exclProb(x[,i],x[,j])
      Pagb[i,j] <- condProb(x[,i],x[,j])
      Nab[i,j] <- length((x[,i]+x[,j])[(x[,i]+x[,j]) == 2])
    }
  }
                                        #significance test?
  if (sig.test){
    Pabexp <- Pa %*% t(Pa) #expected probabilities (null)
    A <- nrow(x) #total number of observations
    Nexp <- A * Pabexp #number of expected co-occurrences
    V <- A * Pabexp * (1 - Pabexp) #variance of co-occurrences
    cl.l <- Nexp - d * sqrt(V) #lower confidence limit
    cl.u <- Nexp + d * sqrt(V) #upper confidence limit
    sig <- Nab < cl.l | Nab > cl.u #test whether the observed co-occurrences are below or above the confidence limits
                                        #zero non-significant co-occurrence patterns
    Pab[sig == FALSE] <- 0
    Paorb[sig == FALSE] <- 0 
    Pagb[sig == FALSE] <- 0
                                        #pack for export
    out <- list(Pa,Pab,Paorb,Pagb,sig)
    names(out) <- c('Pa','Pab','Paorb','Pagb','sig')
  }else{
                                        #pack for export
    out <- list(Pa,Pab,Paorb,Pagb)
    names(out) <- c('Pa','Pab','Paorb','Pagb')
  }
                                        #return
  return(out)
}
                                        #calculate information theoretic statistics
enaInformation <- function(x){}

                                        #test the community organization
library(vegan)
library(sna)
x <- read.csv('~/Documents/Active_Projects/ONC_Lichen/ONC_LCO_data/ONCLichenCooc_12May2011.csv')
                                        #remove upper quadrat
x <- x[x$Quadrat == 'n45.55',]
com <- x[,7:ncol(x)]
tree <- paste(x[,1],x[,2],sep='_')
com.l <- split(com, tree)
geno <- factor(sapply(names(com.l),function(x) strsplit(x,split='_')[[1]][2]))
net.Pab <- lapply(com.l,function(x) probNet(x)$Pab)
net.Paorb <- lapply(com.l,function(x) probNet(x)$Paorb)
net.Pagb <- lapply(com.l,function(x) probNet(x)$Pagb)
net.sig <- lapply(com.l,function(x) probNet(x)$sig)
gplot(net.Pagb[[5]],edge.lwd=(net.Paorb[[2]]+1)^3,vertex.cex=2,vertex.sides=25,displaylabels=TRUE)
                                        #network stats
net <- net.Paorb
deg <- unlist(lapply(net,function(x) (length(x[x != 0]) - length(diag(x)[diag(x) != 0])))) #degree
cen <- unlist(lapply(net,function(x) centralization(x,degree))) #centralization

                                        #analysis
rough <- na.omit(read.csv('~/lco/data/ONC_raw_roughness.csv')[,c(1,5)])
rough[,1] <- sub('-','.',rough[,1])
rough[,1] <- sub('\\.0','.',rough[,1])

tree. <- sapply(tree,function(x) strsplit(x,split='_')[[1]][1])
for (i in 1:length(rough[,1])){
  tree.[tree. == rough[i,1]] <- rough[i,2]
}
rough <- tapply(as.numeric(tree.),tree,mean)
                                        #analysis
summary(deg.glm <- glm(deg~rough*geno,family=poisson))
summary(cen.glm <- glm(cen~rough*geno,family=))
hist(residuals(deg.glm))
hist(residuals(cen.glm))
shapiro.test(residuals(deg.glm))
shapiro.test(residuals(cen.glm))
                                        #multivariate
                                        #make diagonals zero
net.Pab. <- net.Pab
net.Paorb. <- net.Paorb
net.Pagb. <- net.Pagb
for (i in 1:length(net.Pab)){
  diag(net.Pab.[[i]]) <- 0
  diag(net.Paorb.[[i]]) <- 0
  diag(net.Pagb.[[i]]) <- 0
}
                                        #calculate distances
d.Pab <- as.dist(netDist(net.Pab.))
d.Paorb <- as.dist(netDist(net.Paorb.))
d.Pagb <- as.dist(netDist(net.Pagb.))
adonis(d.Pab~geno*rough)
adonis(d.Paorb~geno*rough)
adonis(d.Pagb~geno*rough)
                                        #PCA
pca.Pab <- princomp(d.Pab)
pca.Paorb <- princomp(d.Paorb)
pca.Pagb <- princomp(d.Pagb)
par(mfrow=c(1,3))
plot(pca.Pab$scores,col=rainbow(nlevels(geno))[as.numeric(geno)],pch=19)
plot(pca.Paorb$scores,col=rainbow(nlevels(geno))[as.numeric(geno)],pch=19)
plot(pca.Pagb$scores,col=rainbow(nlevels(geno))[as.numeric(geno)],pch=19)

                                        #look into the following motif based approaches
dya <- do.call(rbind,lapply(net.Pab,function(x) dyad.census(x)))
(dya)
tri <- do.call(rbind,lapply(net.Pab,function(x) triad.census(x,mode='graph')))
tri
‘triad.classify’, ‘dyad.census’, ‘kcycle.census’, ‘kpath.census’,
     ‘gtrans’
