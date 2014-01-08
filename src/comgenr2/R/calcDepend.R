calcDepend <-
function(a,b){
  return(length(a[(a+b)==2])/sum(a)) #intersection of a with b divided by the total of a
}
