#' Calculates the log-probability of a transmission tree
#' (assumes w is Gamma(2,1)) #TODO relax this
#' @param ttree Transmission tree
#' @param R Basic reproduction number
#' @param pi probability of sampling an infected individual
#' @return Probability of the transmission tree
probTTree = function(ttree,R,pi)  {
  prob <- 0 
  n <- nrow(ttree)
  p0 <- -lambert_W0(R*(pi-1)/exp(R))/R
  pd0 <- rep(NA, n+1)
  for (i in (1:n)) { 
    if (is.na(ttree[i,2])) {prob<-prob+log(1-pi)} else 
    {prob<-prob+log(pi);
    prob <- prob + dgamma((ttree[i,2]-ttree[i,1]), shape = 2, scale = 1, log=TRUE)} 
    offspring <- which( ttree[ ,3] == i ) 
    d0 <- length(offspring)
    if (is.na(pd0[d0+1])) {
      pd0[d0+1]=d0*log(R)+R*(p0-1)-d0*log(R*p0)+log(1-d0*gamma_inc(d0,R*p0)/gamma(1+d0))
    }
    prob <- prob + pd0[d0+1]
    #prob <- prob + dpois(length(offspring),R,log=TRUE)    
    for (j in offspring) {
      prob <- prob + dgamma((ttree[j,1]-ttree[i,1]), shape = 2, scale = 1,log=TRUE) 
    } 
  } 
  return(prob)
} 
