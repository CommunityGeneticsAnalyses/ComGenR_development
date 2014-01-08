getH2C <-
function(x='NMDS axis',g='grouping vector',sibs='sibling proportion (1=clones)'){

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
