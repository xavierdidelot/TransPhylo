#' Simulate an outbreak
#' @param off.r First parameter of the negative binomial distribution for offspring number
#' @param off.p Second parameter of the negative binomial distribution for offspring number
#' @param neg the within-host effective population size (Ne) timesgeneration duration (g)
#' @param nSampled number of sampled infected individuals, or NA for any
#' @param pi probability of sampling an infected individual
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time 
#' @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time
#' @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time 
#' @param dateStartOutbreak Date when index case becomes infected
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @return Combined phylogenetic and transmission tree
#' @examples
#' plotCTree(simulateOutbreak())
#' @export
simulateOutbreak = function(off.r=1,off.p=0.5,neg=0.25,nSampled=NA,pi=0.5,w.shape=2,w.scale=1,ws.shape=w.shape,ws.scale=w.scale,dateStartOutbreak=2000,dateT=Inf) {
  #Create a transmission tree with nSampled infected sampled individuals
  nsam<-0
  nh<-0
  rejected=-1
  while (is.na(nSampled)||nsam!=nSampled) {
    ttree=NULL
    while (is.null(ttree)) {
      mtt<-makeTTree(off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT-dateStartOutbreak,nSampled)
      rejected=rejected+1
      ttree<-mtt$ttree
      if (mtt$pruned>0) {
        dateStartOutbreak=dateStartOutbreak+mtt$pruned
        cat(sprintf('Note that simulated outbreak was pruned: in order to have %d sampled by present date %f, the start date was set to %f\n',nSampled,dateT,dateStartOutbreak))
      }
      }
    nsam<-length(which(!is.na(ttree[,2])))
    nh=nrow(ttree)-nsam
    if (is.na(nSampled)) nSampled=nsam
  }
  if (rejected>0) cat(sprintf('Note that rejection sampling was used %d times to simulate outbreak with %d sampled individuals\n',rejected,nSampled))
  n<-nsam+nh
  
  #Create a within-host phylogenetic tree for each infected host
  wtree<-vector('list',n)
  for (i in (1:n)) {
    if (is.na(ttree[i,2])) {times<-c(           ttree[which(ttree[,3]==i),1])-ttree[i,1]}
                      else {times<-c(ttree[i,2],ttree[which(ttree[,3]==i),1])-ttree[i,1]}
    wtree[[i]]<-.withinhost(times,neg)[[1]]
  }
  
  #Glue these trees together
  truth<-.glueTrees(ttree,wtree)
  truth[,1]<-truth[,1]+dateStartOutbreak
  return(list(ctree=truth,nam=mtt$nam))
}  
