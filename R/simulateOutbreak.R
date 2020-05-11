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
#' @param w.mean Mean of the Gamma distribution representing the generation time
#' @param w.std Std of the Gamma distribution representing the generation time 
#' @param ws.mean Mean of the Gamma distribution representing the sampling time
#' @param ws.std Std of the Gamma distribution representing the sampling time 
#' @param dateStartOutbreak Date when index case becomes infected
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @return Combined phylogenetic and transmission tree
#' @examples
#' simulateOutbreak()
#' simulateOutbreak(off.r=2,dateStartOutbreak=2010,dateT=2015)
#' @export
simulateOutbreak = function(off.r=1,off.p=0.5,neg=0.25,nSampled=NA,pi=0.5,w.shape=2,w.scale=1,ws.shape=NA,ws.scale=NA,w.mean=NA,w.std=NA,ws.mean=NA,ws.std=NA,dateStartOutbreak=2000,dateT=Inf) {
  if (!is.na( w.mean)&&!is.na( w.std)) { w.shape= w.mean^2/ w.std^2; w.scale= w.std^2/ w.mean}
  if (!is.na(ws.mean)&&!is.na(ws.std)) {ws.shape=ws.mean^2/ws.std^2;ws.scale=ws.std^2/ws.mean}
  if (is.na(ws.shape)) ws.shape=w.shape
  if (is.na(ws.scale)) ws.scale=w.scale
  
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
      probttree<-mtt$prob
      if (mtt$pruned>0) {
        dateStartOutbreak=dateStartOutbreak+mtt$pruned
        message(sprintf('Note that simulated outbreak was pruned: in order to have %d sampled by present date %f, the start date was set to %f',nSampled,dateT,dateStartOutbreak))
      }
      }
    nsam<-length(which(!is.na(ttree[,2])))
    nh=nrow(ttree)-nsam
    if (is.na(nSampled)) nSampled=nsam
  }
  if (rejected>0) message(sprintf('Note that rejection sampling was used %d times',rejected))
  n<-nsam+nh
  
  #Create a within-host phylogenetic tree for each infected host
  wtree<-vector('list',n)
  probwithin=0
  for (i in (1:n)) {
    if (is.na(ttree[i,2])) {times<-c(           ttree[which(ttree[,3]==i),1])-ttree[i,1]}
                      else {times<-c(ttree[i,2],ttree[which(ttree[,3]==i),1])-ttree[i,1]}
    a<-withinhost(times,neg)
    wtree[[i]]=a$nodes
    probwithin=probwithin+a$prob
  }
  
  #Glue these trees together
  truth<-.glueTrees(ttree,wtree)
  truth[,1]<-truth[,1]+dateStartOutbreak
  l=list(ctree=truth,nam=mtt$nam,probttree=probttree,probwithin=probwithin)
  class(l)<-'ctree'
  return(l)
}  
