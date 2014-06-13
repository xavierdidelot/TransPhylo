#' Extracts transmission tree from a combined phylogenetic/transmission tree 
#' @param fulltree Combined tree
#' @return transmission tree
ttreeFromFullTree = function(fulltree)  {
  host <- fulltree[ ,4]
  ttree <- fulltree[fulltree[ ,2] == 0&fulltree[ ,3] == 0,1] 
  nsam <- length(ttree) 
  nh <- nrow(fulltree)-3*nsam+1
  ntot <- nsam+nh
  ttree <- c(ttree,rep(NA,nh))
  ttree <- cbind(matrix(0, ntot, 1),ttree,matrix(0, ntot, 1),deparse.level=0) 
  for (i in (1:ntot)) { 
    j <- min(which(host==i)) 
    while (host[j] == i)  { 
      j <- which( fulltree[ ,2] == j | fulltree[ ,3] == j ) 
    } 
    ttree[i,1] <- fulltree[j,1] 
    ttree[i,3] <- host[j] 
  } 
  return(ttree)
} 
