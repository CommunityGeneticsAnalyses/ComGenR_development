library(ComGenR)
test <- cgSim(reps=1,GG=2,YY=1)
com <- test[[1]][[1]][[1]]
geno <- gpmTrees()[,1]
nms <- nmds(vegdist(com),maxdim=1)
getH2C(x=nmds.min(nms)[,1],g=geno,sibs=1)
myh2c <- getH2C(x=nmds.min(nms)[,1],g=geno,sibs=1)
c(myh2c[1] - myh2c[2],myh2c[1],myh2c[1] + myh2c[2])

min.stress <- numeric()
lower.ci <- numeric()
min.r2 <- numeric()
for (i in 1:100){
  nms <- nmds(vegdist(com),maxdim=1)
  min.stress[i] <- min(nms$stress)
  min.r2[i] <- nms$r2[nms$stress==min(nms$stress)]
  myh2c <- getH2C(x=nmds.min(nms)[,1],g=geno,sibs=1)
  lower.ci[i] <- myh2c[1] - myh2c[2]
}

pairs(cbind(min.stress,min.r2,lower.ci))
summary(lm(lower.ci~min.r2))
