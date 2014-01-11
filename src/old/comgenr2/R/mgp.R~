mgp <- function(net='species network',com='community matrix',my.coord='',loc=TRUE,v.scale=3,v.min=0.1){

v.cex <- apply(com[,apply(com,2,sum)!=0],2,sum) #scaling node size by the log of species frequencies
v.cex <- (((v.cex/sum(v.cex))/max((v.cex/sum(v.cex))))*v.scale)+v.min
e.col <- net
e.col[net>0.5] <- 'red'
e.col[net<0.5] <- 'black'
e.col[net==0.5] <- 'grey'

if (length(my.coord)==1){
  coord <- gplot(abs(net),displaylabels=FALSE,gmode='graph',pad=1.5,
        edge.lwd=(abs(net)),vertex.cex=v.cex,vertex.col='grey',
        edge.col=e.col)

}else{
  gplot(abs(net),displaylabels=FALSE,gmode='graph',pad=1.5,
        edge.lwd=(abs(net)),vertex.cex=v.cex,vertex.col='grey',
        edge.col=e.col,coord=my.coord)
}
if (loc){return(coord)}else{}
}
