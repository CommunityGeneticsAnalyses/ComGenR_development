gpmCom <-
function(n='number of species',het.values=c(5,21),allelic.range=c(0,3)){
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
