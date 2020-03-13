#' Simulate a transmission tree
#' @param off.r First parameter of the negative binomial distribution for offspring number
#' @param off.p Second parameter of the negative binomial distribution for offspring number
#' @param pi probability of sampling an infected individual
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time 
#' @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time
#' @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time 
#' @param maxTime Duration of simulation (can be Inf)
#' @param nSampled Number of sampled individuals (can be NA for any)
#' @return A N*3 matrix in the following format with one row per infected host, first column is time of infection, second column is time of sampling, third column is infector
#' @export
makeTTree <-function(off.r,off.p,pi,w.shape,w.scale,ws.shape=w.shape,ws.scale=w.scale,maxTime=Inf,nSampled=NA) { 
  ttree<-matrix(0,1,3)
  prob<-0
  todo<-1
  while (length(todo)>0) {
    rand<-runif(1)
    if (rand<pi) {
      #This individual is sampled
      prob<-prob+log(pi)
      draw<-rgamma(1,shape=ws.shape,scale=ws.scale)
      if (ttree[todo[1],1]+draw<maxTime)
      ttree[todo[1],2]<-ttree[todo[1],1]+draw
      else ttree[todo[1],2]<-NA
      prob<-prob+log(dgamma(draw,shape=ws.shape,scale=ws.scale))}
    else {
      #This individual is not sampled
      prob<-prob+log(1-pi)
      ttree[todo[1],2]<-NA}
    
    offspring<-rnbinom(1,off.r,1-off.p)
    prob<-prob+log(dnbinom(offspring,off.r,1-off.p))
    if (offspring>0) {
      for (i in 1:offspring) {
        draw<-rgamma(1,shape=w.shape,scale=w.scale)
        prob<-prob+log(dgamma(draw,shape=w.shape,scale=w.scale))
        if (ttree[todo[1],1]+draw>maxTime) next
        ttree<-rbind(ttree,c(ttree[todo[1],1]+draw,0,todo[1]))
        todo<-c(todo,nrow(ttree))
      }
    }
    todo<-todo[-1] 
  }
  
  #Prune down number of sampled individuals if needed (only for ongoing outbreaks)
  pruned=0
  samTimes=ttree[which(!is.na(ttree[,2])),2]
  if (!is.na(nSampled)&&maxTime<Inf&&nSampled<length(samTimes)) {
    samTimes=sort(samTimes)
    cutoff=(samTimes[nSampled]+samTimes[nSampled+1])/2
    pruned=maxTime-cutoff
    ttree[ttree[,2]>cutoff,2]=NA
  }
  
  #Remove infected individuals who are not sampled and have no children
  while (TRUE) {
    if (nrow(ttree)==1 && is.na(ttree[1,2])) {return(list(ttree=NULL,prob=NULL,pruned=0))} #Nothing left
    torem=c()
    for (i in 1:nrow(ttree)) {if (is.na(ttree[i,2])&&length(which(ttree[,3]==i))==0) {torem=c(torem,i)}}
    if (length(torem)==0) {break}
    ttree<-ttree[setdiff(1:nrow(ttree),torem),,drop=FALSE]
    for (i in 1:nrow(ttree)) {ttree[i,3]=ttree[i,3]-length(which(torem<ttree[i,3]))}
  }  
  
  #Remove root if not sampled and only has one child
  while (TRUE) {
    if (nrow(ttree)==1 && is.na(ttree[1,2])) {return(list(ttree=NULL,prob=NULL,pruned=0))} #Nothing left
    if (is.na(ttree[1,2])&&length(which(ttree[,3]==1))==1) {
      ttree[,2]<-ttree[,2]-ttree[2,1]
      ttree[,1]<-ttree[,1]-ttree[2,1]
      ttree<-ttree[-1,,drop=FALSE]
      ttree[,3]=ttree[,3]-1
    } else {break}
  }
  
  #Reorder so that sampled hosts are first
  order<-c(which(!is.na(ttree[,2])),which(is.na(ttree[,2])))
  invorder=1:length(order);invorder[order]=1:length(order)
  ttree<-ttree[order,,drop=FALSE]
  ttree[ttree[,3]>0,3]=invorder[ttree[ttree[,3]>0,3]]
  l=list(ttree=ttree,nam=sprintf('%d',seq(1:length(which(!is.na(ttree[,2]))))),prob=prob,pruned=pruned)
  class(l)<-'ttree'
  return(l)
}
