#' Simulate a transmission tree
#' (assumes generation distribution w is Gamma(2,1)) TODO relax this
#' @param R the basic reproduction number
#' @param pi probability of sampling an infected individual
#' @return A N*3 matrix in the following format with one row per infected host, first column is time of infection, second column is time of sampling, third column is infector
makeTTree <-function(R,pi) { 
  ttree<-matrix(0,1,3)
  prob<-0
  todo<-1
  while (length(todo)>0) {
    rand<-runif(1)
    if (rand<pi) {
      #This individual is sampled
      prob<-prob+log(pi)
      draw<-rgamma(1,2,1)
      ttree[todo[1],2]<-ttree[todo[1],1]+draw
      prob<-prob+log(dgamma(draw,2,1))}
    else {
      #This individual is not sampled
      prob<-prob+log(1-pi)
      ttree[todo[1],2]<-NA}
    offspring<-rpois(1,R)
    prob<-prob+log(dpois(offspring,R))
    if (offspring>0) {
      for (i in 1:offspring) {
        draw<-rgamma(1,2,1)
        prob<-prob+log(dgamma(draw,2,1))
        ttree<-rbind(ttree,c(ttree[todo[1],1]+draw,0,todo[1]))
        todo<-c(todo,nrow(ttree))
        if (nrow(ttree)>100) {return(list(ttree=NULL,prob=NULL))}
      }
    }
    todo<-todo[-1] 
  }
  
  #Remove infected individuals who are not sampled and have no children
  while (TRUE) {
    if (nrow(ttree)==1 && is.na(ttree[1,2])) {return(list(ttree=NULL,prob=NULL))} #Nothing left
    torem=c()
    for (i in 1:nrow(ttree)) {if (is.na(ttree[i,2])&&length(which(ttree[,3]==i))==0) {torem=c(torem,i)}}
    if (length(torem)==0) {break}
    ttree<-ttree[setdiff(1:nrow(ttree),torem),,drop=FALSE]
    for (i in 1:nrow(ttree)) {ttree[i,3]=ttree[i,3]-length(which(torem<ttree[i,3]))}
  }  
  
  #Remove root if not sampled and only has one child
  while (TRUE) {
    if (nrow(ttree)==1 && is.na(ttree[1,2])) {return(list(ttree=NULL,prob=NULL))} #Nothing left
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
  return(list(ttree=ttree,prob=prob))
}