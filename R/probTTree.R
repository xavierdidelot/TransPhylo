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
  omegaStar <- uniroot(function(x) {x-(1-pi)*((1-off.p)/(1-off.p*x))^off.r},c(0,1))$root #This is Equation (1)
  alphaStar <- rep(NA, n+1)
  for (i in (1:n)) { 
    if (is.na(ttree[i,2])) prob<-prob+log(1-pi) #This is the first term in the product in Equation (5)
    else prob<-prob+log(pi)+dgamma((ttree[i,2]-ttree[i,1]),shape=w.shape,scale=w.scale,log=TRUE) #This is the second term in the product in Equation (5)
    offspring <- which( ttree[ ,3] == i ) 
    d0 <- length(offspring)
    if (is.na(alphaStar[d0+1])) {
      alphaStar[d0+1]=0
      notinf=d0+2*off.r*off.p/(1-off.p)
      for (k in d0:notinf) alphaStar[d0+1]=alphaStar[d0+1]+choose(k,d0)*dnbinom(k,off.r,off.p)*omegaStar^{k-d0} #This is in Equation (2)
      alphaStar[d0+1]=log(alphaStar[d0+1])
    }
    prob <- prob + alphaStar[d0+1] #This is the third term in the product in Equation (5)
    for (j in offspring) {
      prob <- prob + dgamma((ttree[j,1]-ttree[i,1]),shape=w.shape,scale=w.scale,log=TRUE) #This is the fourth term in the product in Equation (5)
    } 
  } 
  return(prob)
} 
