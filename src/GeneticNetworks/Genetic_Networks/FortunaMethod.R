#Fortuna et al. 2009's Method for Resolving Independence Structure in Spatial Genetic Networks

#Using the crust community data from Matt's Spain dataset: /Users/artemis/Documents/Active_Projects/SpainCrustNetworks/SpainCrustAnalysis.R

#Using Bray-Curtis Distance
library(vegan)
library(corpcor)

#My microsatellite distance function (see Fortuna et al. 2009 using Smouse and Peakall 1999)

#1 Translate the genotype of an individual into a codification vector Y of length K, where K is the number of k alleles in the population. The y values can be 0, 1 or 2 depending on whether the individual has zero one or two copies of the k allele

#a. separate the loci
loci=list()
z=1:ncol(x)
z=z[which(z %% 2 != 0)]
for (i in z){
	loci[[i]]=paste(colnames(x)[i],x[,i],sep='.')
	loci[[(i+1)]]=paste(colnames(x)[(i)],x[,(i+1)],sep='.')
	}
names(loci)=colnames(x)
#b. get a vector of all of the alleles in the population
K=unique(as.vector(unlist(loci)))

#c. Determine if an individual has zero, one or two copies of each allele for the codification vector Y
x.indiv=data.frame(loci)
Y=array(NA,c(nrow(x.indiv),length(K)))

for (i in 1:nrow(x.indiv)){
	for (j in 1:length(K)){
	Y[i,j]=length(x.indiv[i,x.indiv[i,]==K[j]])
	}
	}

#Independence network function (also from Fortuna et al. 2009)

d=as.matrix(vegdist(t(x)))
	
source('/Users/Aeolus/Documents/Active_Projects/GeneticNetworks/Genetic_Networks/ind_net.R')

par(mfrow=c(4,6))

for (i in 1:length(com.l)){
  ind.net(com.l[[i]])
}

