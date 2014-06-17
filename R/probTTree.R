#' Calculates the log-probability of a transmission tree
#' (assumes w is Gamma(2,1)) #TODO relax this
#' @param ttree Transmission tree
#' @param R Basic reproduction number
#' @param pi probability of sampling an infected individual
#' @return Probability of the transmission tree
probTTree = function(ttree,R,pi)  {
  #TODO need to account for unobserved branches, cf writelatex file
  prob <- 0 
  for (i in (1:nrow(ttree))) { 
    if (is.na(ttree[i,2])) {prob<-prob+log(1-pi)} else 
    {prob<-prob+log(pi);
    prob <- prob + log(dgamma((ttree[i,2]-ttree[i,1]), shape = 2, scale = 1))} 
    offspring <- which( cbind(ttree[ ,3]) == i ) 
    prob <- prob + log(dpois(length(offspring),R)) 
    for (j in (offspring)) {
      prob <- prob + log(dgamma((ttree[j,1]-ttree[i,1]), shape = 2, scale = 1)) 
    } 
  } 
  return(prob)
} 
