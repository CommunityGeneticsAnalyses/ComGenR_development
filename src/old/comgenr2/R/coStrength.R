coStrength <-
function(x='dependency network',direction='in'){
  if (direction=='in'){
    return(apply(x,2,sum))    
  }else{
    return(apply(x,1,sum))
  }
}
