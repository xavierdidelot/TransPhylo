#' Converts an ape phylo object into a phylogenetic tree
#' @param tr phylo object
#' @param dateLastSample date of the last sample
#' @return phylogenetic tree
#' @export
ptreeFromPhylo <- function(tr,dateLastSample) {
  if (dateLastSample>1900 && dateLastSample<2100 && sum(tr$edge.length)<1)
    warning('Warning: input tree has small branch lengths. This needs to be a timed tree in the same unit as dateLastSample (eg years).\n')
  n<-length(tr$tip.label)
  ed<-tr$edge
  le<-tr$edge.length
  
  #Preserve leaf names if they are already named from 1 to n
  #na<-tr$tip.label
  #co<-as.character(1:n)
  #if (length(setdiff(na,co))==0 && length(setdiff(co,na))==0) {
  #  na<-as.numeric(na)
  #  for (i in 1:nrow(ed)) if(ed[i,2]<=n) ed[i,2]=na[ed[i,2]]
  #}
  
  tra<-c(1:n,(2*n-1):(n+1))
  ptree<-matrix(0,2*n-1,3)
  if (n==1) {ptree[1,1]=dateLastSample;return(ptree)}
  for (i in 1:nrow(ed)) {
    father<-tra[ed[i,1]]
    son<-tra[ed[i,2]]
    if (ptree[father,2]==0) {ptree[father,2]=son} else {ptree[father,3]=son}
    ptree[son,1]<-le[i]
  }
  todo<-2*n-1
  while (length(todo)>0) {
    t1=todo[1]
    if (ptree[t1,2]==0) {todo=todo[-1];next}
    ptree[ptree[t1,2],1]<-ptree[ptree[t1,2],1]+ptree[t1,1]
    ptree[ptree[t1,3],1]<-ptree[ptree[t1,3],1]+ptree[t1,1]
    todo=c(todo[-1],ptree[t1,2],ptree[t1,3])
  }
  ptree[,1]=ptree[,1]-max(ptree[,1])+dateLastSample
  l=list(ptree=ptree,nam=tr$tip.label)
  class(l)<-'ptree'
  return(l)
}
