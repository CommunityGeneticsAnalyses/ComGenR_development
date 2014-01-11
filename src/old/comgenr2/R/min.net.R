min.net <-
function(net='network',com='community matrix'){
  out <- list(net=net,com=com)
  if (all(net!='network')){out[[1]] <- net[apply(net,1,sum)>0,apply(net,2,sum)>0]}
  if (all(com!='community matrix')){out[[2]] <- com[,apply(net,2,sum)>0]}
  return(out)
}
