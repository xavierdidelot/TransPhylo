#' Extracts transmission tree from a combined phylogenetic/transmission tree 
#' @param ctree Combined tree
#' @return transmission tree
#' @export
extractTTree = function(ctree)  {
  nam=ctree$nam
  ctree=ctree$ctree
  host <- ctree[ ,4]
  ttree <- ctree[ctree[ ,2] == 0&ctree[ ,3] == 0,1] 
  nsam <- length(ttree) 
  nh <- nrow(ctree)-3*nsam+1
  ntot <- nsam+nh
  ttree <- c(ttree,rep(NA,nh))
  ttree <- cbind(matrix(0, ntot, 1),ttree,matrix(0, ntot, 1),deparse.level=0) 
  for (i in (1:ntot)) { 
    j <- min(which(host==i)) 
    while (host[j] == i)  { 
      j <- which( ctree[ ,2] == j | ctree[ ,3] == j ) 
    } 
    ttree[i,1] <- ctree[j,1] 
    ttree[i,3] <- host[j] 
  } 
  return(list(ttree=ttree,nam=nam))
} 
