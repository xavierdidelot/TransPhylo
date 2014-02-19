#' Converts a phylogenetic tree into an ape phylo object
#' @param ptree phylogenetic tree
#' @return phylo object
phyloFromPtree <- function(ptree) {
  n<-ceiling(nrow(ptree)/2)
  tr<-list()
  tr$Nnode<-n-1
  tr$tip.label<-as.character(1:n)
  tr$edge<-matrix(0,n*2-2,2)
  tr$edge.length<-rep(0,n*2-2)
  iedge<-1
  tra<-c(1:n,(2*n-1):(n+1))
  for (i in nrow(ptree):(n+1)) {
    tr$edge[iedge,]<-c(tra[i],tra[ptree[i,2]])
    tr$edge.length[iedge]<-ptree[ptree[i,2],1]-ptree[i,1]
    iedge<-iedge+1
    tr$edge[iedge,]<-c(tra[i],tra[ptree[i,3]])
    tr$edge.length[iedge]<-ptree[ptree[i,3],1]-ptree[i,1]
    iedge<-iedge+1
  }
  class(tr)<-'phylo'
  return(tr)
}