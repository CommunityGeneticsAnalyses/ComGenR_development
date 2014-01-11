#Fortuna et al. 2009's Method for Resolving Independence Structure in Spatial Genetic Networks

#Using the crust community data from Matt's Spain dataset: /Users/artemis/Documents/Active_Projects/SpainCrustNetworks/SpainCrustAnalysis.R

#Using Bray-Curtis Distance
library(vegan)
x=com.l[[1]]
d=as.matrix(vegdist(t(x)))

c=d*0

c_ij=function(d='distance matrix'){
	
	
	
	}
