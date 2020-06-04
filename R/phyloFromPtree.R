#' Converts a phylogenetic tree into an ape phylo object
#' @param ptree phylogenetic tree
#' @return phylo object
#' @examples 
#' phyloFromPTree(extractPTree(simulateOutbreak()))
#' @export
phyloFromPTree <- function(ptree) {
  nam=ptree$nam
  ptree=ptree$ptree
  n<-ceiling(nrow(ptree)/2)
  if (n==1) return(ape::read.tree(text='(1);'))
  tr<-list()
  tr$Nnode<-n-1
  tr$tip.label<-nam
  tr$edge<-matrix(0,n*2-2,2)
  tr$edge.length<-rep(0,n*2-2)
  iedge<-1
  root<-which(ptree[,1]==min(ptree[,1]))
  tra<-c(1:n,root,setdiff((n+1):(2*n-1),root))
  tra2<-1:length(tra)
  tra[tra]<-tra2
  for (i in (n+1):(2*n-1)) {
    tr$edge[iedge,]<-c(tra[i],tra[ptree[i,3]])
    tr$edge.length[iedge]<-ptree[ptree[i,3],1]-ptree[i,1]
    iedge<-iedge+1
    tr$edge[iedge,]<-c(tra[i],tra[ptree[i,2]])
    tr$edge.length[iedge]<-ptree[ptree[i,2],1]-ptree[i,1]
    iedge<-iedge+1
  }
  tr$root.time=min(ptree[,1])
  class(tr)<-'phylo'
  tr=ape::reorder.phylo(tr, order = "cladewise")
  return(tr)
}
