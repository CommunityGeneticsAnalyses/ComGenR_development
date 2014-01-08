mgp <-
function(scn,com,loc=TRUE,my.coord=''){
e.col <- sign(scn)
e.col[e.col==1] <- 'grey'
e.col[e.col==-1] <- 'red'
v.cex <- apply(com,2,sum)
v.cex <- log(v.cex,10)
v.cex <- v.cex * (1/min(v.cex))
v.cex <- v.cex/2
if (length(my.coord)==1){
  coord <- gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
                 edge.col=e.col,edge.lwd=abs(scn),
                 vertex.cex=v.cex,vertex.col='darkgrey',vertex.border='darkgrey')
}else{
  coord <- gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
                 edge.col=e.col,edge.lwd=abs(scn),
                 vertex.cex=v.cex,vertex.col='darkgrey',vertex.border='darkgrey',
                 coord=my.coord)
}
if (loc){return(coord)}else{}
}
