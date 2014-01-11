ch.plot <-
function(x='ordination matrix',g='groupings',cex=1,buffer=0.1,plot.legend=TRUE,loc='topleft'){
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
