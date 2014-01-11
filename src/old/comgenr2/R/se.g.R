se.g <-
function(x='community matrix',g='grouping'){
  return(apply(x,2,function(x,g) tapply(x,g,function(x) sd(x)/sqrt(length(x))),g=g))
}
