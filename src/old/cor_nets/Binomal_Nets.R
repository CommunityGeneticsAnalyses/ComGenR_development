###Binomal Analysis of Species Co-occurrence Data

library(limma) #From the bioconductor repository source("http://www.bioconductor.org/biocLite.R")
library(gdata)

#library(venn) #installed from http://fisher.stats.uwo.ca/faculty/murdoch/repos/html/vennv1.5.htmlm
#This function did not work.

#Binomial distribution
my.pbinom <- function(k,n,p){choose(n,k)*p^k*(1-p)^(n-k)}

plot.binom <- function(n,p,by=1,pch=19,lty=1,xlim=c(0,n),p.col=1,l.col=1,ylim=c(0,0.2),empty=TRUE){p
if (empty == TRUE){p.type <- 'n'}else{p.type <- 'p'}
plot(my.pbinom(seq(0,n,by=by),n,p)~seq(0,n,by=by),type=p.type,pch=pch,col=p.col,ylim=ylim,xlab='',ylab='')
if (empty == TRUE){}else{lines(my.pbinom(seq(0,n,by=by),n,p)~seq(0,n,by=by),lty=lty,col=l.col)}
}

plot.bdist <- function(A,B,ylim=c(0,0.35)){
AB <- A+B

##For the AB species example
by <- 1
Na <- length(A[A == 1]) #number of times A occurs in total is the number of trials
Ka <- length(AB[AB == 2]) #number of times A and B occur together is a success
Nb <- length(B[B == 1]) #number of times B occurs in total is the number of trials
Kb <- length(AB[AB == 2]) #number of times A and B occur together is a success
p <- 0.5 #NULL HYPOTHESIS success and failure have equal probabilities

##Generate a plot of the null distribution and the observed statistics
par(mfrow=c(1,2))
plot.binom(max(Na,Nb),p,p.col=1,ylim=ylim)
lines(my.pbinom(seq(0,Na,by=by),Na,p)~seq(0,Na,by=by),lty=1,col=2)
points(Ka,my.pbinom(Ka,Na,p),col=2,pch=1,cex=2);points(Ka,my.pbinom(Ka,Na,p),col='black',pch=19,cex=0.5)
lines(my.pbinom(seq(0,Nb,by=by),Nb,p)~seq(0,Nb,by=by),lty=1,col=3)
points(Kb,my.pbinom(Kb,Nb,p),col=3,pch=1,cex=2);points(Kb,my.pbinom(Kb,Nb,p),col='black',pch=19,cex=0.5)
legend('topright',legend=c('A','B'),fill=c(2,3))
title(xlab='k = Number of Sucesses',ylab='Probability')

#Generate a Venn Diagram
ven.dat <- cbind(A==1,B==1)
colnames(ven.dat) <- c('A','B')
ven.count <- vennCounts(ven.dat)
vennDiagram(ven.count,circle.col=c(2,3),lwd=2)

}

mcbinom <- function(x,y,perm=999){
xy <- x + y
s. <- numeric()
s <- length(xy[xy==2])

for (i in 1:perm){
x. <- sample(x)
y. <- sample(y)
xy. <- x. + y.
s.[i] <- length(xy.[xy. == 2]) #calculate simulated success (i.e. where x and y co-occur)
}

out <- c(s,mean(s.),length(s.[s.<=s])/length(s.),length(s.[s.>=s])/length(s.))
names(out) <- c('obs','mu.sim','p.lower','p.upper')

return(out)
}


##General example from Wikipedia binomial distribution page
plot.binom(20,0.5,p.col=3,ylim=c(0,0.35))
n=20;p=0.5;by=1
lines(my.pbinom(seq(0,n,by=by),n,p)~seq(0,n,by=by),lty=1,col=3)
p=0.1
lines(my.pbinom(seq(0,n,by=by),n,p)~seq(0,n,by=by),lty=1,col=4)
p=0.8
lines(my.pbinom(seq(0,n,by=by),n,p)~seq(0,n,by=by),lty=1,col=2)
title(xlab='k = Number of Sucesses',ylab='Probability')

##Start
n <- 100 #number of observations
A <- rbinom(n,1,0.5) #sample for species A
B <- rbinom(n,1,0.5) #sample for species B
AB <- A+B

plot.binom(A,B)
binom.test(length(AB[AB==2]),length(A[A==1]))
binom.test(length(AB[AB==2]),length(B[B==1]))


##Setup scenarios for different levels of dependence
n <- 100 #number of observations
A <- rbinom(n,1,0.1) #sample for species A
B <- rbinom(n,1,0.9) #sample for species B
AB <- A+B

plot.bdist(A,B,ylim=c(0,0.35))
binom.test(length(AB[AB==2]),length(A[A==1]))
binom.test(length(AB[AB==2]),length(B[B==1]))

##The problem with this is that it cannot be applied when the population sizes are un-equal
##Try a Monte-Carlo approach.

mcbinom(A,B)
plot.bdist(A,B)
