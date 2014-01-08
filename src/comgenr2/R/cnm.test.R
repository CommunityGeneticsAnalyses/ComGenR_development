cnm.test <-
function(com,method='r1',nits=5000,burn=500,thin=10){
  ###Co-occurrence Null Modeling test
  com.nul <- nullCom(com,method=method,nits=nits,burn=burn,thin=thin)
  cs.nul <- unlist(pblapply(com.nul,cscore))
  cs.obs <- cscore(com)
  ses <- (cs.obs-mean(cs.nul))/sd(cs.nul)
  return(c(SES=ses,lower.p=length(cs.nul[cs.nul<=cs.obs])/length(cs.nul),
           upper.p=length(cs.nul[cs.nul>=cs.obs])/length(cs.nul)))
}
