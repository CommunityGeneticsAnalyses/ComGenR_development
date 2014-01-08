nullCom <-
function(com,method='r1',nits=5000,burn=500,thin=10){
                                        #force binary
  com[com!=0] <- 1
  for (i in 1:burn){post.burn <- commsimulator(x=com,method=method,thin=thin)}
  out <- list()
  for (i in 1:nits){out[[i]] <- commsimulator(x=post.burn,method=method,thin=thin)}
  return(out)
}
