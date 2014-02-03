###Easy single and two random effect REML
##3 Feb 2014

library(lme4)
cgREML <- function(x,g,cv='covariate'){
  if (any(cv=='covariate')){
    x.lmer <- lmer(x~(1|g))
    x.lm <- lm(x~1)
    chi2 <- -2*logLik(x.lm, REML=T) +2*logLik(x.lmer, REML=T)
    p.chi2 <- pchisq(chi2,df=1,lower.tail=FALSE)
  }else{
    x.lmer <- lmer(x~(1|g)+cv)
    x.lm <- lm(x~1)
    chi2 <- -2*logLik(x.lm, REML=T) +2*logLik(x.lmer, REML=T)
    p.chi2 <- pchisq(chi2,df=1,lower.tail=FALSE)
  }
  return(c(chi2=chi2,P.value=p.chi2))
}
