#' Simulate an outbreak
#' @param R the basic reproduction number
#' @param neg the within-host effective population size (Ne) times  generation duration (g)
#' @param ninf number of sampled infected individuals, or NA for any
#' @param pi probability of sampling an infected individual
#' @return Combined phylogenetic and transmission tree
#' @examples
#' plotBothTree(simulateOutbreak())
simulateOutbreak = function(R=1,neg=0.25,ninf=NA,pi=0.5,dateLastSample=2014) {
  #Create a transmission tree with ninf infected sampled individuals
  nsam<-0
  nh<-0
  while (is.na(ninf)||nsam!=ninf) {
    ttree<-makeTTree(R,pi)[[1]]
    if (is.null(ttree)) {nsam<-0;nh=0} else {
      nsam<-length(which(!is.na(ttree[,2])))
      nh=nrow(ttree)-nsam
      if (is.na(ninf)) ninf=nsam
    }
  }
  n<-nsam+nh
  
  #Create a within-host phylogenetic tree for each infected host
  wtree<-vector('list',n)
  for (i in (1:n)) {
    if (is.na(ttree[i,2])) {times<-c(           ttree[which(ttree[,3]==i),1])-ttree[i,1]}
                      else {times<-c(ttree[i,2],ttree[which(ttree[,3]==i),1])-ttree[i,1]}
    wtree[[i]]<-.withinhost(times,neg)[[1]];
  }
  
  #Glue these trees together
  truth<-.glueTrees(ttree,wtree)
  truth[,1]<-truth[,1]-max(truth[,1])+dateLastSample
  return(truth)
}  
