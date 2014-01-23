###Source simTrees and simSpp


###Make the trees
tree.gpm <- gpmTrees()
trees <- simTrees(tree.gpm=tree.gpm,VeT=2)

###Make the insects
insects <- simSpp()

###Increase selection intensity
sim.z <- list()
for (i in 1:8){
  sim.z[[i]] <- cgSim(trees=trees,insects=insects,z=i,VeN=45)
}

##Look at gamma values
gamma.shuster06 <- read.csv('../data/shuster_2006_points_c0.csv',header=TRUE)
gamma.shuster06
gamma.cgsim <- round(unlist(lapply(sim.z,function(x) x$gamma)),5)
cor(gamma.cgsim[1:nrow(gamma.shuster06)],gamma.shuster06[,1])

##Get heritabilities
#h2c <- unlist(lapply(sim.z,function(x,g) getH2C(nmds.min(nmds(vegdist(x[[1]]),2,2))[,1],g)[2],g=tree.gpm[,1]))
#h2c.out <- cbind(h2c.0,h2c.15,h2c.30,h2c.45,h2c.60)

                                        #this is for no tree environmental variation
h2c <- read.csv(file='../results/h2c_zvsEn.csv')[,-1]
h2c[h2c<0] <- 0
plot(gamma.cgsim,h2c[,1],log='x',pch='',xaxt='none',xlim=c(0.00001,0.1),ylim=c(0,1),xlab='Gamma',ylab='H2C')
axis(1,at=c(0.00001,0.0001,0.001,0.01,0.1),labels=c(0.00001,0.0001,0.001,0.01,0.1))
h2c.pch <- c(18,22,17,1,8)
for (i in 1:ncol(h2c)){
  points(gamma.cgsim,h2c[,i],pch=h2c.pch[i],type='b')
}
pairs(h2c)

##Networks
sim.z <- list()
for (i in 1:8){
  sim.z[[i]] <- cgSim(trees=trees,insects=insects,z=i,VeN=15)
}
nets <- lapply(lapply(sim.z,function(x) x[[1]]),CoNetwork,threshold=0)
par(mfrow=c(2,4))
for (i in 1:length(nets)){mgp(nets[[i]],sim.z[[i]][[1]])}
nest.nets <- lapply(sim.z,function(x,g) nestedness(mean.g(x[[1]],g)),g=tree.gpm[,1])
nest.check <- oecosimu(mean.g(sim.z[[8]][[1]],tree.gpm[,1]),nestedtemp,method='r1',nsimul=1000,burnin=100)
nestedtemp(mean.g(sim.z[[7]][[1]],tree.gpm[,1]))
