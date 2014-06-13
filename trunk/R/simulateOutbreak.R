#' Simulate an outbreak
#' @param R the basic reproduction number
#' @param neg the within-host effective population size (Ne) times  generation duration (g)
#' @param ninf number of infected individuals
#' @param pi proportion of sampled infected individuals
#' @return Combined phylogenetic and transmission tree
#' @examples
#' plotBothTree(simulateOutbreak())
simulateOutbreak = function(R=1,neg=0.25,ninf=10,pi=0.5) {
  #Create a transmission tree with ten individuals
  n<-1
  while (n!=ninf) {
    ttree<-makeTTree(R)[[1]]
    n<-nrow(ttree)
    if (is.null(ttree)) {n<-0}
  }
  
  #Create a within-host phylogenetic tree for each infected host
  wtree<-vector('list',n)
  for (i in (1:n)) {
    times<-c(ttree[i,2],ttree[which(ttree[,3]==i),1])-ttree[i,1]
    wtree[[i]]<-.withinhost(times,neg)[[1]];
  }
  
  #Glue these trees together
  truth<-.glueTrees(ttree,wtree)
  truth[,1]<-truth[,1]+2005#Epidemic started in 2005
  return(truth)
}  
