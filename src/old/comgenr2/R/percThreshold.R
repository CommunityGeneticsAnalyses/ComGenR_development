percThreshold <-
function(x='network matrix',step.size=0.01){
  no.c <- no.clusters(graph.adjacency(x,weighted=TRUE))
  step <- 1
  while (no.c==1){
    x[x<=(step*step.size)] <- 0
    no.c <- no.clusters(graph.adjacency(x,weighted=TRUE))
    step <- step + 1
  }
  out <- list(threshold=((step-1)*step.size),isolated.nodes=colnames(x)[clusters(graph.adjacency(x,weighted=TRUE))$membership==2])
  return(out)
}
