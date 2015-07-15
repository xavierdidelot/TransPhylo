#' Calculates the log-probability of a transmission tree
#' @param ttree Transmission tree
#' @param off.r First parameter of the negative binomial distribution for offspring number
#' @param off.p Second parameter of the negative binomial distribution for offspring number
#' @param pi probability of sampling an infected individual
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation length w
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation length w 
#' @return Probability of the transmission tree
probTTree = function(ttree,off.r,off.p,pi,w.shape,w.scale)  {
  prob <- 0 
  n <- nrow(ttree)
  f <- function(x) {x-(1-pi)*((1-off.p)/(1-off.p*x))^off.r}
  p0 <- uniroot(f,c(0,1))$root
  pd0 <- rep(NA, n+1)
  for (i in (1:n)) { 
    if (is.na(ttree[i,2])) prob<-prob+log(1-pi) else 
      prob<-prob+log(pi)+dgamma((ttree[i,2]-ttree[i,1]),shape=w.shape,scale=w.scale,log=TRUE) 
    offspring <- which( ttree[ ,3] == i ) 
    d0 <- length(offspring)
    if (is.na(pd0[d0+1])) {pd0[d0+1]=0
    notinf=d0+2*off.r*off.p/(1-off.p)
    for (k in d0:notinf) pd0[d0+1]=pd0[d0+1]+choose(k,d0)*dnbinom(k,off.r,off.p)*p0^{k-d0}
    pd0[d0+1]=log(pd0[d0+1])
    }
      #R*(p0-1)-d0*log(p0)+pgamma(R*p0,d0,log.p=T)#This is log of alpha_*(d)
    prob <- prob + pd0[d0+1]
    #prob <- prob + dpois(length(offspring),R,log=TRUE)    
    for (j in offspring) {
      prob <- prob + dgamma((ttree[j,1]-ttree[i,1]),shape=w.shape,scale=w.scale,log=TRUE) 
    } 
  } 
#  print(c(p0,exp(pd0[1])*(1-pi)))#Sanity check: should be the same
  return(prob)
} 
