if (all(ls()!='f.list')){
  require(codetools)
  library(enaR)
  called.by <- function(tarFunc, tarPack){
    flist <-   sapply(lsf.str(tarPack, all=TRUE), c)
    names(flist) <- NULL
    gotit <- sapply(flist, function(x) tarFunc %in% findGlobals(get(x, tarPack),FALSE)$functions)
    flist[gotit]
  }

  f.list <- as.character(sapply(lsf.str('package:ComGenR',all=TRUE),c))
  f.array <- array(0,dim=rep(length(f.list),2))
  rownames(f.array) <- colnames(f.array) <- f.list
  for (i in 1:length(f.list)){
    f.array[match(called.by(f.list[i],'package:ComGenR'),rownames(f.array)),i] <- 1
  }
  f.net <- network(t(f.array))
}

plot(f.net,displaylabels=TRUE,label.cex=0.85,arrowhead.cex=0.65,
     edge.lwd=0.75,vertex.col='lightblue',vertex.border='white',edge.col='darkgrey')

