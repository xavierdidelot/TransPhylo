#' Extracts transmission tree from a combined phylogenetic/transmission tree 
#' @param ctree Combined tree
#' @return transmission tree
#' @examples 
#' extractTTree(simulateOutbreak())
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
  parents <- rep(NA, nrow(ctree))
  parents[ctree[ ,2:3] + 1] <- 1:nrow(ctree)
  parents=parents[-1]
  maxs=rep(0,ntot)
  for (i in 1:length(host)) maxs[host[i]]=i
  for (i in (1:ntot)) { 
    j<-maxs[i]#(which(host==i))
    j<-parents[j]
    ttree[i,1] <- ctree[j,1] 
    ttree[i,3] <- host[j] 
  } 
  l=list(ttree=ttree,nam=nam)
  class(l)<-'ttree'
  return(l)
} 
