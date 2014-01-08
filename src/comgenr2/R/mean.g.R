mean.g <-
function(x='community matrix',g='grouping'){
  return(apply(x,2,function(x,g) tapply(x,g,mean),g=g))
}
