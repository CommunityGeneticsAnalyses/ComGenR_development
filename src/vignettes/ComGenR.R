###ComGenR - A community genetics analysis package

## package.skeleton(name='ComGenR',path='~/projects/ComGenR')

## Community genetics seeks to identify the community level impacts of
## evolutionary processes. Thus, this field merges analyses across
## several fields including population genetics and community ecology. In
## this package, I present a core set of functions that provide
## analytical tools for community genetics analyses that can be grouped
## broadly into three categories:

## I. Compositional Analyses
## II. Quantifying Heritability 
## III. Network and Co-occurrence Modeling

## Create an analytical package that can be downloaded from github
## Use your standard project folder structure
## Instruct the user that there is an analytical template in the src folder

## Package takes the input of a community matrix and a vector of genotype codes
## The package is designed for two levels of analysis: 1) Stand level and 2) Tree level
## stand level = replicates on a tree by tree basis
## stand level is expected to be much more common
## data are handled as matrices
## tree level = replicates within an individual tree
## data input as matrices with a tree indicator and converted internally to lists

## 0. Dependencies

require(vegan)
require(ecodist)
require(sna)
require(bipartite)
require(pbapply)

## I. Compositional Analyses

## PerMANOVA (including dispersion tests)
## %%Vegan - adonis


## Pairwise permanovas
## %%internal - use code from trihybrid analyses and web analyses
pair.permanova <- function(x,f,nits=999){
  require(vegan)
  f. <- sort(unique(f))
  out <- list()
  p.out <- array(NA,dim=c(length(f.),length(f.)))
  rownames(p.out) <- colnames(p.out) <- f.
  h <- 1
  for (i in 1:length(f.)){
    k <- i + 1
    for (j in k:length(f.)){
      if (i!=j & j<=length(f.)){
        print(paste(f.[i],f.[j],sep=' vs '))
        y <- x[f == f.[i] | f == f.[j],]
        yf <- factor(f[f == f.[i] | f == f.[j]])
        out[[h]] <- as.matrix(adonis(y~yf)$aov.tab,permutations=nits)
        p.out[i,j] <- out[[h]][1,dim(out[[h]])[2]]
        names(out)[h] <- paste(f.[i],f.[j],sep=' vs ')
        h <- h + 1
      }else{}
    }
  }
  out <- list(f.tables=out,p.mat=p.out)
  return(out)
}

## NMDS (include test using scores)
## %%Ecodist - nmds

## Cross-hair plots
## %%Internal/Vegan

###cross hair plots
## Usage
## xtype <- as.character(x.[,colnames(x.)=="Genotype..FINAL."])
## xtype[xtype=='BC1-DxAA '] <- 'BC1-DxAA'
## xtype[xtype=='trihybride'] <- 'Trihybrid'
## ch.plot(nmds.min(nms.),xtype,cex=2,plot.legend=TRUE,loc='topright')

ch.plot <- function(x='ordination matrix',g='groupings',cex=1,plot.legend=FALSE,loc='topleft'){
  mu <- apply(x,2,function(x,g) tapply(x,g,mean),g=g)
  se <- apply(x,2,function(x,g) tapply(x,g,function(x) sd(x)/sqrt(length(x))),g=g)
  mu <- na.omit(mu)
  se <- na.omit(se)
                                        #coloring
  mu.col <- 'black'
  mu.pch <- 19
                                        #error bars
  cl.xu <- mu[,1] + se[,1]
  cl.xl <- mu[,1] - se[,1]
  cl.yu <- mu[,2] + se[,2]
  cl.yl <- mu[,2] - se[,2]
  plot(mu,pch=mu.pch,cex=cex,xlim=c(min(cl.xl),max(cl.xu)),ylim=c(min(cl.yl),max(cl.yu)),col=mu.col)
  for (i in 1:nrow(mu)){
    lines(x=c(cl.xl[i],cl.xu[i]),y=c(mu[i,2],mu[i,2]))
    lines(x=c(mu[i,1],mu[i,1]),y=c(cl.yl[i],cl.yu[i]))
  }
  if (plot.legend){legend(loc,legend=rownames(se),cex=cex*0.5,pch=mu.pch,col=mu.col,border='grey')}else{}
}

## II. Quantifying Heritability

###Genotype stats function

mean.g <- function(x='community matrix',g='grouping'){
  return(apply(x,2,function(x,g) tapply(x,g,mean),g=g))
}

se.g <- function(x='community matrix',g='grouping'){
  return(apply(x,2,function(x,g) tapply(x,g,function(x) sd(x)/sqrt(length(x))),g=g))
}

## Shuster's calculator
## %%Internal - h2c


getH2C <- function(x='NMDS axis',g='grouping vector',sibs='sibling proportion (1=clones)'){

g <- as.character(g)
aov.tab <- matrix(unlist(summary(aov(x~g))),nrow=2)
if (sibs=='sibling proportion (1=clones)'){sibs <- 1}

#S = number of clones, families, sires, etc.
S <- length(unique(as.character(g))) 
#ni = number of individuals in ith clone
ni <- table(g)
#n. = total number of inidividuals in the analysis
n. <- sum(ni)
#k = ni in expected means squares 
#= ni if design is balanced 
#= kl if design is unbalanced
if (all(ni==ni[1])){
  k <- mean(ni)
}else{
  ##this is kl
  k <- 1/(S-1)*(n.-(sum(ni^2)/n.))
}
#s = sigma here
#variance among = s2s = (MSs-MSw)/k
#variance within = s2w = MSw
#total variance = Vp = s2total = s2s + s2w
#H2c = community heritability
s2s <- (aov.tab[1,3]-aov.tab[2,3])/k
s2w <- aov.tab[2,3]
s2total <- s2s+s2w
H2C <- s2s / s2total

###Confidence limits for H2C
t <- H2C * (1/sibs)
if (all(ni==ni[1])){
  SE <- ((2*(((1-t)^2)*(1+(k-1)*t)^2))/(k*(k-1)*(S-1)))^(1/2)
}else{
  ##Unbalanced design with kl
  SE <- ((2*(n.-1)*((1-t)^2)*(1+(k-1)*t)^2)/((k^2)*(n.-S)*(S-1)))^(1/2)
}
CI <- SE*1.96

return(c(H2C=H2C,CI=CI,SE=SE))

}

## Community Simulator from Shuster 2006
## %%Internal - cgSim

###Function for simulating community genetics patterns 
###based on the function from Shuster et al. 2006.
###Started 13Dec2013 but taken from work begun in 2012
### tree genotypic values - range from 11 to 21

gpmTrees <- function(pheno='phenotype of genotypes',reps='replication'){
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

###Initiate arthropods

gpmCom <- function(n='number of species',het.values=c(5,21),allelic.range=c(0,3)){
if (n=='number of species'){n <- 25}
                                        #generate heterozygous alleles for n species
com <- matrix(NA,nrow=n,ncol=2)
com[,1] <- runif(n,het.values[1],het.values[2]) #heterozygote value between 5 and 21
com[,2] <- runif(n,allelic.range[1],allelic.range[2]) #range between 0 and 3
                                        #map genotype to phenotype
com. <- com
com.[,1] <- (com[,1]-0.5*com[,2])/2 #C allelic value = (HET - 0.5*range)/2
com.[,2] <- (com[,1]+0.5*com[,2])/2 #D allelic value = (HET + 0.5*range)/2

return(com.)

}

##Initiate tree

cgSim <- function(tree.pheno='tree phenotype matrix',insect='community phenotype matrix',reps=10,GG=8,YY=5,VeT=8,Ve=0.1,VeN=15,K=100){

if (any(tree.pheno=='tree phenotype matrix')){tree.pheno <- gpmTrees()}
if (any(insect=='community phenotype matrix')){insect <- gpmCom()}

###Experimental design
# reps <- 10 #number of times to run simulation
# GG <- 8 #number of selection scenarios
# YY <- 5 #number of environmental scenarios
###Parameter values
# VeT <- 8 #environmental variance in tree trait influences tree heritability (2 for high H2 and 8 for low H2)
# Ve <- 0.1 #environmental variance for insect trait
# VeN <- 15 #step size for environmental variance in interactions (0 15 30 45 60)
# K <- 100 #insect pop carrying capacity

###Initiate output objects
T <- nrow(tree.pheno) #number of trees
I <- nrow(insect) #number of insects
art_g <- matrix(NA,nrow=T,ncol=I) #insect gene frequency
art_z <- matrix(NA,nrow=T,ncol=I) #insect trait mean
art_Vg <- matrix(NA,nrow=T,ncol=I) #insect gene variability
art_Vz <- matrix(NA,nrow=T,ncol=I) #insect trait variability
gen_load <- matrix(NA,nrow=T,ncol=I) #insect genetic load
dem_load <- matrix(NA,nrow=T,ncol=I) #insect demographic load
art_pop <- matrix(NA,nrow=T,ncol=I) #insect population community matrix
                                        #STORAGE LIST
out <- gg.list <- yy.list <- list()
###Simulation
tic <- Sys.time() #start stopwatch
for (RR in 1:reps){
###generating trees to use
  scores_XX <- matrix(NA,nrow=nrow(tree.pheno),ncol=GG)
  scores_XX[,1] <- tree.pheno[,2] + runif(T,0,1) * VeT - VeT/2 #scores are in second column of trees
###filling phenotypic values for trees (same set of trees for all scenarios
  for (z in 2:GG){
    scores_XX[1:T,z] <- scores_XX[1:T,1]
  }
  trees <- scores_XX #just renaming phenotypes to be trees
  for (y in 1:YY){
                                        #YY variation scenarios of other ecological interactions
    for (z in 1:GG){
                                        #GG selection intensity scenarios
      for (i in 1:T){
                                        #insects on trees
                                        #for each tree i
        for (j in 1:I){
                                        #for each insect j
###Equation 6 from MS - Arthropod (art) Gene frequency
          if (trees[i,z] < 2*insect[j,1]){
            art_g[i,j] <- 0
          }else if (trees[i,z] > 2*insect[j,2]){
            art_g[i,j] <- 1
          }else{
            art_g[i,j] <- (trees[i,z] - 2*insect[j,1]) / (2*insect[j,2] - 2*insect[j,1])
          }
###Equation 6 from MS - calculating mean trait Z from art gene frequency
          art_z[i,j] <- 2*insect[j,2]*art_g[i,j]^2 + 2*art_g[i,j]*(1-art_g[i,j])*(insect[j,1]+insect[j,2])+2*insect[j,1]*(1-art_g[i,j])^2 + runif(1)*Ve - Ve/2
###Art genetic (Vg) and trait variance (Vz)
          art_Vg[i,j] <- 2*art_g[i,j]*(1-art_g[i,j])
          art_Vz[i,j] <- art_Vg[i,j]*(insect[j,2] - insect[j,1])^2+Ve
###Evolutionary (gen) and demographic (dem) loads from selection
          gen_load[i,j] <- 0.5*(0.00007924*2.511886^(z-1))*(art_z[i,j]- trees[i,z])^2
          dem_load[i,j] <- 0.5*(0.00007924*2.511886^(z-1))*(art_Vz[i,j])
###Equation 7 from MS - art predicted population size as a function of loads and ecological variance
          art_pop[i,j] <- K * (1 - gen_load[i,j] - dem_load[i,j])+runif(1)*VeN*(y-1)-VeN*(y-1)/2
###preventing art pops from going negative or zero slightly above zero
          if (art_pop[i,j] < 0){
            art_pop[i,j] <- runif(1)*3;
          }else{
            art_pop[i,j] #not sure if this is correct
          }
        } #end insect loop
      } #end tree loop
###keeping track of which iteration the simulation is on reps*GG*(y-1)+reps*(z-1)+RR
###dlmwrite(int2str(reps*GG*(y-1)+reps*(z-1)+RR),art_pop,'\t');  %%% write arthropod sampling to file %%%
                                        #print(paste(reps*GG*(y-1),reps*(z-1),RR,sep='_'))
      print(paste(RR,y,z,sep=' '))
      gg.list[[z]] <- art_pop
                                        # names(gg.list)[z] <- paste(reps*GG*(y-1),reps*(z-1),RR,sep='')
      names(gg.list)[z] <- paste(RR,y,z,sep='_')
    } #end selection loop
    yy.list[[y]] <- gg.list
  } #end environment loop
  toc <- Sys.time() #stop stopwatch
  out[[RR]] <- yy.list
} #end REP loop and END simulation

return(out)
}



## III. Network and Co-occurrence Modeling

## Araujo's Method
## %%internal - coNet

null.prune <- function(a,b,std=TRUE){

###Method for producing network models from co-occurrence data.
###Developed in Araujo et al. 2011 Ecography Using species co-occurrence networks to assess the impacts of climate change
###Coded 18 Sep 2013 MKLau

## E is an environmental lattice with rows = N and cols = M and has size A = N*M.
## A = total number of sites
## Na = number of sites with species a
                                        #determine number of sites
A <- length(a)
                                        #NULL probabilities
## Assuming no interactions among species
## Probability of observing species a
## Pa <- Na/A
## Probability of observing both species a and b simultaneously
## Pa&b = Pa * Pb
## Probability of observing either species a and b
## Pa+b = Pa(1-Pb) + (1-Pa)Pb = Pa + Pb - 2PaPb 
## Probability of observing neither species a nor b
## P!a&!b = (1-Pb)(1-Pa) = 1 - Pa - Pb + PaPb
## and it follows that 
## Pa&b + Pa+b + P!a&!b = 1

pa <- sum(a)/A
pb <- sum(b)/A
paANDb <- pa * pb
paORb <- pa + pb - 2*pa*pb
pNOTaORb <- 1- pa - pb + pa*pb

## Three states, i = ab, ii a+b = (a&!b + !a&b), iii = !a!b
                                        #expected frequencies
## E(i) = APa&b
## E(ii) = APa+b
## E(iii) = AP!a&!b

Ei <- A*paANDb
Eii <- A*paORb
Eiii <- A*pNOTaORb
                                        #variance of frequencies
## Var(i) = APa&b(1-Pa&b)
## Var(ii) = APa+b(1-APa+b)
## Var(iii) = AP!a&!b(1-AP!a&!b)

Vi <- A*paANDb*(1-paANDb)
Vii <- A*paORb*(1-paORb)
Viii <- A*pNOTaORb*(1-pNOTaORb)
                                        
## Critical interval
## Critical threshold (+ = attraction and - = repulstion)
## CI.u = E(i) + 2*Var(i)^2 = attractive overlap
## CI.l = E(i) - 2*Var(i)^2 = repulsive overlap

ci.u <- Ei + 2*Vi^(1/2)
ci.l <- Ei - 2*Vi^(1/2)
                                        #Oi = observed a AND b
Oi <- length(a[(a+b)==2])
                                        #Determine output (Oi = significant, 0 = not)
if (Oi > ci.u | Oi < ci.l){
  if (std){
    out <- (Oi - Ei) / sqrt(Vi)
  }else{
    out <- Oi
  }

}else{
  out <- 0
}
return(out)
}

co.net <- function(x='species in cols',diag.zero=TRUE,std=TRUE){
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

##Species Dependence
calcDepend <- function(a,b){
  return(length(a[(a+b)==2])/sum(a)) #intersection of a with b divided by the total of a
}

###Dependency Network
dep.net <- function(x='species in cols',zero.na=TRUE,prune=TRUE,diag.zero=TRUE){
  out <- matrix(NA,nrow=ncol(x),ncol=ncol(x))
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
        out[i,j] <- calcDepend(x[,i],x[,j])
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

percThreshold <- function(x='network matrix',step.size=0.01){
  no.c <- no.clusters(graph.adjacency(x,weighted=TRUE))
  step <- 1
  while (no.c==1){
    x[x<=(step*step.size)] <- 0
    no.c <- no.clusters(graph.adjacency(x,weighted=TRUE))
    step <- step + 1
  }
  out <- list(threshold=((step-1)*step.size),isolated.nodes=colnames(x)[clusters(graph.adjacency(x,weighted=TRUE))$membership==2])
  return(out)
}


CoNetwork <- function(x='community matrix'){
###Runs all steps of the process for modeling
###Co-occurrence networks described by Araujo et al. 2011.
###It depends on the seenetR.R script which contains both the
###Araujo functions and related co-occurrence null modeling
###functions.

###Inputs: 
#x = matrix of co-occurrence with species in columns

#Step 1. Calculate a Bray-Curtis distance matrix

bc.d <- as.matrix(vegdist(t(x)))

#Step 2. Prune distance matrix based on co-occurrence probabilities

prune <- co.net(x)
bc.d[prune==0] <- 0

#Step 3. Reduce to percolation threshold
thresh <- percThreshold(bc.d)$threshold
pruned.net <- bc.d
pruned.net[bc.d<thresh] <- 0

return(pruned.net)

}

## Network plots
min.net <- function(net='network',com='community matrix'){
  out <- list(net=net,com=com)
  if (all(net!='network')){out[[1]] <- net[apply(net,1,sum)>0,apply(net,2,sum)>0]}
  if (all(com!='community matrix')){out[[2]] <- com[,apply(net,2,sum)>0]}
  return(out)
}

mgp <- function(scn,com,loc=TRUE,my.coord=''){
e.col <- sign(scn)
e.col[e.col==1] <- 'grey'
e.col[e.col==-1] <- 'red'
v.cex <- apply(com,2,sum)
v.cex <- log(v.cex,10)
v.cex <- v.cex * (1/min(v.cex))
v.cex <- v.cex/2
if (length(my.coord)==1){
  coord <- gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
                 edge.col=e.col,edge.lwd=abs(scn),
                 vertex.cex=v.cex,vertex.col='darkgrey',vertex.border='darkgrey')
}else{
  coord <- gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
                 edge.col=e.col,edge.lwd=abs(scn),
                 vertex.cex=v.cex,vertex.col='darkgrey',vertex.border='darkgrey',
                 coord=my.coord)
}
if (loc){return(coord)}else{}
}


## Co-ocurrence using appropriate null model
## %%vegan/internal - nullSim

###Calculate Stone and Roberts C-score
cscore <- function(x,cu.mat=FALSE){
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



nullCom <- function(com,method='r1',nits=5000,burn=500,thin=10){
                                        #force binary
  com[com!=0] <- 1
  for (i in 1:burn){post.burn <- commsimulator(x=com,method=method,thin=thin)}
  out <- list()
  for (i in 1:nits){out[[i]] <- commsimulator(x=post.burn,method=method,thin=thin)}
  return(out)
}

cnm.test <- function(com,method='r1',nits=5000,burn=500,thin=10){
  ###Co-occurrence Null Modeling test
  com.nul <- nullCom(com,method=method,nits=nits,burn=burn,thin=thin)
  cs.nul <- unlist(pblapply(com.nul,cscore))
  cs.obs <- cscore(com)
  ses <- (cs.obs-mean(cs.nul))/sd(cs.nul)
  return(c(SES=ses,lower.p=length(cs.nul[cs.nul<=cs.obs])/length(cs.nul),
           upper.p=length(cs.nul[cs.nul>=cs.obs])/length(cs.nul)))
}


## %%internal - cgNet
## %%stand level option
## Model stand level network using all data
## Model stand level network using compressed data
## %%tree level option
## Model tree level networks using repeated measurements within
##   trees
## Average tree level networks
## %%bipartite
## Model tree level co-occurrence and return C-score, SES and p-values
## Model bipartite networks between genotypes and species

## %%sna - qap.test
## Test for correlations between networks

## %%internal - netDist
## Calculate network distances
netDist <- function(dn.t){
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

## %%vegan or bipartite

## Network Structural Analyses

## Test for nestedness of bipartite networks
## %%igraph?
## Test for modularity 

###Araujo network statistics

#species strength
coStrength <- function(x='dependency network',direction='in'){
  if (direction=='in'){
    return(apply(x,2,sum))    
  }else{
    return(apply(x,1,sum))
  }
}

#symmetry
coSym <- function(x='dependency network',zero.na=TRUE){
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


