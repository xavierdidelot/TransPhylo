#' Calculates the log-probability of a transmission tree
#' (assumes w is Gamma(2,1)) #TODO relax this
#' @param ttree Transmission tree
#' @param R Basic reproduction number
#' @param pi probability of sampling an infected individual
#' @return Probability of the transmission tree
probTTree = function(ttree,R,pi)  {
  #TODO need to account for unobserved branches, cf writelatex file
  prob <- 0 
  n <- nrow(ttree)
  for (i in (1:n)) { 
    if (is.na(ttree[i,2])) {prob<-prob+log(1-pi)} else 
    {prob<-prob+log(pi);
    prob <- prob + dgamma((ttree[i,2]-ttree[i,1]), shape = 2, scale = 1, log=TRUE)} 
    offspring <- which( cbind(ttree[ ,3]) == i ) 
    prob <- prob + dpois(length(offspring),R,log=TRUE)
    for (j in (offspring)) {
      prob <- prob + dgamma((ttree[j,1]-ttree[i,1]), shape = 2, scale = 1,log=TRUE) 
    } 
  } 
  return(prob)
} 
