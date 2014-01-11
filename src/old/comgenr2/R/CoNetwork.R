CoNetwork <-
function(x,plot.net=TRUE,scalar=3,min.vsize=0.1){
###Runs all steps of the process for modeling
###Co-occurrence networks described by Araujo et al. 2011.
###It depends on the seenetR.R script which contains both the
###Araujo functions and related co-occurrence null modeling
###functions.

###Inputs: 
#x = matrix of co-occurrence with species in columns
#plot.net = logical. Should the network be plotted?
#scalar = scales the size of all vertices
#min.vsize = sets the minimum size for vertices

#Step 1. Calculate a Bray-Curtis distance matrix

bc.d <- as.matrix(vegdist(t(x)))

#Step 2. Prune distance matrix based on co-occurrence probabilities

prune <- co.net(x)
bc.d[prune==0] <- 0

#Step 3. Reduce to percolation threshold
thresh <- percThreshold(bc.d)$threshold
pruned.net <- bc.d
pruned.net[bc.d<thresh] <- 0

if (plot.net){
  v.cex <- apply(x,2,sum) #scaling node size by the log of species frequencies
  v.cex <- (((v.cex/sum(v.cex))/max((v.cex/sum(v.cex))))*scalar)+min.vsize
  gplot(abs(pruned.net),displaylabels=TRUE,gmode='graph',pad=1.5,
        edge.lwd=(abs(pruned.net)),vertex.cex=v.cex,vertex.col='grey')
}

return(pruned.net)

}
