#' Select the most representative transmission tree from a MCMC output
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return The index of the selected transmission tree
#' @export 
selectTTree = function(record,burnin=0.5)
{
  #Remove burnin
  burned=round(length(record)*burnin)
  record=record[(burned+1):length(record)]
  m=length(record)
  n=sum(record[[1]]$ctree$ctree[,2]==0&record[[1]]$ctree$ctree[,3]==0) #Number of sampled individuals

  #Record partitions in sampled transmission tree
  hash=vector('list',n*m*10)
  tinf=matrix(0,m,n)
  a=floor(n*10*runif(n))+1
  for (i in 1:m)
  {
    ttree=extractTTree(record[[i]]$ctree)$ttree
    tinf[i,]=ttree[1:n,1]
    #Find children of nodes
    todo=1:n
    children=matrix(NA,nrow(ttree),n)
    while (length(todo)>0) {
      w=todo[1]
      todo=todo[-1]
      direct=which(ttree[,3]==w)
      children[w,]=colSums(children[direct,,drop=F])>0
      if (is.na(children[w,1])) next
      if (w<=n) children[w,w]=1
      if (ttree[w,3]>0) todo=c(todo,ttree[w,3])
    }
    #Add nodes one by one
    for (j in 1:nrow(ttree)) {
      c=children[j,]
      hh=1+(sum(a[which(c==1)])%%(length(hash)))
      found=0
      while (1) {
        if (is.null(hash[[hh]])) break
        if (all(c==hash[[hh]]$c)) break
        hh=hh+1
        if (hh>length(hash)) hh=1
      }
      if (is.null(hash[[hh]])) {hash[[hh]]$c=c;hash[[hh]]$n=0;hash[[hh]]$w=c();hash[[hh]]$last=0;hash[[hh]]$t=c()}
      if (hash[[hh]]$last==i) {hash[[hh]]$w[1]=hash[[hh]]$w[1]+1}
      else {
        hash[[hh]]$last=i
        hash[[hh]]$n=hash[[hh]]$n+1
        hash[[hh]]$w=c(1,hash[[hh]]$w)
        hash[[hh]]$t=c(i,hash[[hh]]$t)
      }
    }
    
    #Add bonus leaves
    for (j in 1:n) {
      if ((sum(children[j,]))==1) next
      c=rep(0,n);c[j]=1
      hh=1+(sum(a[which(c==1)])%%(length(hash)))
      found=0
      while (1) {
        if (is.null(hash[[hh]])) break
        if (all(c==hash[[hh]]$c)) break
        hh=hh+1
        if (hh>length(hash)) hh=1
      }
      if (is.null(hash[[hh]])) {hash[[hh]]$c=c;hash[[hh]]$n=0;hash[[hh]]$w=c();hash[[hh]]$last=0;hash[[hh]]$t=c()}
      hash[[hh]]$n=hash[[hh]]$n+1
      hash[[hh]]$w=c(0,hash[[hh]]$w)
      hash[[hh]]$t=c(i,hash[[hh]]$t)
    }
  }
  
  comb=hash[!sapply(hash,is.null)]
  dist=matrix(0,m,m)
  for (i in 1:length(comb)) {
    a=comb[[i]]$t
    b=setdiff(1:m,a)
    s=mean(comb[[i]]$w)^2
    dist[a,b]=dist[a,b]+s
    dist[b,a]=dist[b,a]+s
  }
  dist=rowSums(dist)
  
  dist2=matrix(0,m,m)
  for (i in 2:m) for (j in 1:i) {dist2[i,j]=sum((tinf[i,]-tinf[j,])^2);dist2[j,i]=dist2[i,j]}
  dist2=rowSums(dist2)
  ws=which(dist==min(dist))
  dist2=dist2[ws]
  w2=which(dist2==min(dist2))
  elected=ws[w2]
  elected=elected[1]
  
  return(burned+elected)
}
