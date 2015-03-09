#' Simulate an outbreak
#' @param R the basic reproduction number
#' @param neg the within-host effective population size (Ne) times  generation duration (g)
#' @param ninf number of sampled infected individuals, or NA for any
#' @param pi probability of sampling an infected individual
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation length w
#' @param w.scape Scale parameter of the Gamma probability density function representing the generation length w 
#' @param dateStartOutbreak Date when index case becomes infected
#' @param dataPresent Date when process stops (this can be Inf for fully simulated outbreaks)
#' @return Combined phylogenetic and transmission tree
#' @examples
#' plotBothTree(simulateOutbreak())
simulateOutbreak = function(R=1,neg=0.25,ninf=NA,pi=0.5,w.shape=2,w.scale=1,dateStartOutbreak=2000,datePresent=Inf) {
  #Create a transmission tree with ninf infected sampled individuals
  nsam<-0
  nh<-0
  while (is.na(ninf)||nsam!=ninf) {
    ttree=NULL
    while (is.null(ttree)) ttree<-makeTTree(R,pi,w.shape,w.scale,datePresent-dateStartOutbreak)[[1]]
    nsam<-length(which(!is.na(ttree[,2])))
    nh=nrow(ttree)-nsam
    if (is.na(ninf)) ninf=nsam
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
  truth[,1]<-truth[,1]+dateStartOutbreak
  return(truth)
}  
