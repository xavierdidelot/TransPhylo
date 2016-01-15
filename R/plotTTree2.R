#' Plot a transmission tree without showing all missing links
#' @param ttree Transmission tree
#' @param showLabels Whether or not to show the labels 
plotTTree2 = function(ttree,showLabels=TRUE) {
  tr=list()
  tr$Nnode=length(comb)-n
  tr$tip.label=as.character(1:n)
  edge=cbind(parents[which(!is.na(parents))],which(!is.na(parents)))
  tr$edge=edge;tr$edge[edge==n+1]=root;tr$edge[edge==root]=n+1#Make the root the (n+1)-th
  tr$edge.length=rep(NA,nrow(tr$edge))
  for (i in 1:nrow(tr$edge)) tr$edge.length[i]=median(comb[[tr$edge[i,2]]]$w)
  class(tr)<-'phylo'
  plot(tr)
}