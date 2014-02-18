#' Calculates the log-probability of a transmission tree
#' (assumes w is Gamma(2,1)) 
#' @param ttree Transmission tree
#' @param R Basic reproduction number
#' @return Probability of the transmission tree
probTTree = function(ttree,R)  {
  prob <- 0 
  for (i in (1:nrow(ttree))) { 
    prob <- prob + log(dgamma((ttree[i,2]-ttree[i,1]), shape = 2, scale = 1)) 
    offspring <- which( cbind(ttree[ ,3]) == i ) 
    prob <- prob + log(dpois(length(offspring),R)) 
    for (j in (offspring)) {
      prob <- prob + log(dgamma((ttree[j,1]-ttree[i,1]), shape = 2, scale = 1)) 
    } 
  } 
  return(prob)
} 
