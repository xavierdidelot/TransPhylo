#' Extracts transmission tree from a combined phylogenetic/transmission tree 
#' @param fulltree Combined tree
#' @return transmission tree
ttreeFromFullTree = function(fulltree)  {
  host <- fulltree[ ,4]
  ttree <- fulltree[fulltree[ ,2] == 0&fulltree[ ,3] == 0,1] 
  n <- length(ttree) 
  ttree <- cbind(matrix(0, length(ttree), 1),ttree,matrix(0, length(ttree), 1),deparse.level=0) 
  for (i in (1:n)) { 
    j <- i 
    while (host[j] == i)  { 
      j <- which( fulltree[ ,2] == j | fulltree[ ,3] == j ) 
    } 
    ttree[i,1] <- fulltree[j,1] 
    ttree[i,3] <- host[j] 
  } 
  return(ttree)
} 
