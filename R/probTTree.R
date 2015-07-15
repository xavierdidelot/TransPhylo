#' Calculates the log-probability of a transmission tree
#' @param ttree Transmission tree
#' @param R Basic reproduction number
#' @param pi probability of sampling an infected individual
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation length w
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation length w 
#' @return Probability of the transmission tree
probTTree = function(ttree,R,pi,w.shape,w.scale)  {
  prob <- 0 
  n <- nrow(ttree)
  p0 <- -lambert_W0(R*(pi-1)/exp(R))/R
  pd0 <- rep(NA, n+1)
  for (i in (1:n)) { 
    if (is.na(ttree[i,2])) prob<-prob+log(1-pi) else 
      prob<-prob+log(pi)+dgamma((ttree[i,2]-ttree[i,1]),shape=w.shape,scale=w.scale,log=TRUE) 
    offspring <- which( ttree[ ,3] == i ) 
    d0 <- length(offspring)
    if (is.na(pd0[d0+1])) pd0[d0+1]=R*(p0-1)-d0*log(p0)+pgamma(R*p0,d0,log.p=T)    
    prob <- prob + pd0[d0+1]
    #prob <- prob + dpois(length(offspring),R,log=TRUE)    
    for (j in offspring) {
      prob <- prob + dgamma((ttree[j,1]-ttree[i,1]),shape=w.shape,scale=w.scale,log=TRUE) 
    } 
  } 
#  print(c(p0,exp(pd0[1])*(1-pi)))#Sanity check: should be the same
  return(prob)
} 
