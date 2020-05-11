#' Build a consensus transmission tree from a MCMC output
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @param minimum Minimum probability for inclusion of a partition in the consensus
#' @param debug Used for debugging
#' @return The consensus transmission tree
#' @export
consTTree = function(record,burnin=0.5,minimum=0.2,debug=F)
{
  #Remove burnin
  if (burnin>0) record=record[round(length(record)*burnin):length(record)]
  m=length(record)
  n=sum(record[[1]]$ctree$ctree[,2]==0&record[[1]]$ctree$ctree[,3]==0) #Number of sampled individuals

  #Record partitions in sampled transmission tree
  hash=vector('list',n*m*10)
  a=floor(n*10*runif(n))+1
  for (i in 1:length(record))
  {
    ttree=extractTTree(record[[i]]$ctree)$ttree
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
      if (hash[[hh]]$last==i) {hash[[hh]]$w[1]=hash[[hh]]$w[1]+1;hash[[hh]]$t[1]=max(hash[[hh]]$t[1],ttree[j,1])}
      else {
        hash[[hh]]$last=i
        hash[[hh]]$n=hash[[hh]]$n+1
        hash[[hh]]$w=c(1,hash[[hh]]$w)
        hash[[hh]]$t=c(ttree[j,1],hash[[hh]]$t)
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
      hash[[hh]]$t=c(ttree[j,1],hash[[hh]]$t)
    }
  }
  
  #Choose partitions to include in consensus transmission tree
  comb=hash[!sapply(hash,is.null)]
  if (minimum>0) comb=comb[sapply(comb,'[[',2)>(minimum*m)]
  comb=comb[order(sapply(comb,'[[',2),decreasing=T)]
  keep=c()
  for (i in 1:n) {
    c=rep(0,n);c[i]=1
    for (j in 1:length(comb)) if (all(c==comb[[j]]$c)) {keep=c(keep,j);break}
  }
  for (i in 1:length(comb)) {
    if (is.element(i,keep)) next
    c1=which(comb[[i]]$c==1)
    exclude=F
    for (j in keep) {
      c2=which(comb[[j]]$c==1)
      if (length(intersect(c1,c2))>0 && length(setdiff(c1,c2))>0 && length(setdiff(c2,c1))>0) {exclude=T;break}
    }
    if (exclude==F) keep=c(keep,i)
  }
  comb=comb[keep]
  inds=c()
  for (i in 1:n) for (j in 1:length(comb)) if (sum(comb[[j]]$c)==1&&comb[[j]]$c[i]) inds=c(inds,j)
  for (j in 1:length(comb)) if (sum(comb[[j]]$c)==n) inds=c(inds,j)
  comb=comb[c(inds,setdiff(1:length(comb),inds))]
    
  #Build list of parents from selected partitions
  parents=rep(NA,length(comb))
  bralen =rep(NA,length(comb))
  inftim =rep(NA,length(comb))
  for (i in 1:length(comb)) {
    bralen[i]=round(mean(comb[[i]]$w)*comb[[i]]$n/m)
    inftim[i]=sum(comb[[i]]$t)/comb[[i]]$n
    ci=which(comb[[i]]$c==1)
    bestscore=Inf
    for (j in 1:length(comb)) {
      if (i==j) next
      cj=which(comb[[j]]$c==1)
      if (length(setdiff(ci,cj))>0) next
      if (length(setdiff(cj,ci))<bestscore) {bestscore=length(setdiff(cj,ci));parents[i]=j}
    }
  }
  
  #Plot transmission tree as a phylogenetic tree; this is only for the purpose of debugging
  if (debug) {
    for (i in 1:length(comb)) message(comb[[i]]$c,', w=',mean(comb[[i]]$w),', n=',comb[[i]]$n,', bralen=',round(mean(comb[[i]]$w)*comb[[i]]$n/m))
    tr=list()
    tr$Nnode=length(comb)-n
    tr$edge=cbind(parents[which(!is.na(parents))],which(!is.na(parents)))
    tr$edge.length=bralen[which(!is.na(parents))]
    tr$tip.label=as.character(1:n)
    #tr$tip.label=record[[1]]$ctree$nam
    class(tr)<-'phylo'
    plot(tr)
  }
  
  #Update vectors parents and bralen and inftim so that branches of length zero are removed and branches of length>1 are deduplicated
  i=1
  while (i<=length(parents)) {
    if (bralen[i]==0 && parents[i]>n) {
      torem=parents[i]
      parents[i]=parents[torem]
      parents[which(parents==torem)]=i
      parents[which(parents>torem)]=parents[which(parents>torem)]-1
      parents=parents[-torem]
      bralen[i]=bralen[torem]
      bralen=bralen[-torem]
      inftim[i]=inftim[torem]
      inftim=inftim[-torem]
      i=1
      next
    } 
    if (bralen[i]==0 && i>n) {
      torem=i
      parents[which(parents==torem)]=parents[torem]
      parents[which(parents>torem)]=parents[which(parents>torem)]-1
      parents=parents[-torem]
      bralen[i]=bralen[torem]
      bralen=bralen[-torem]
      inftim[i]=inftim[torem]
      inftim=inftim[-torem]
      i=1
      next
    } 
    if (bralen[i]>1) {
      dt=abs(inftim[i]-inftim[parents[i]])
      inftim=c(inftim,inftim[i]-dt/bralen[i])
      bralen=c(bralen,bralen[i]-1)
      bralen[i]=1
      parents=c(parents,parents[i])
      parents[i]=length(bralen)
      i=1
      next
    }
    i=i+1
  }

  #Build transmission tree
  cons=matrix(NA,length(parents),3)
  cons[,1]=inftim
  parents[which(is.na(parents))]=0
  cons[,3]=parents
  cons[1:n,2]=ttree[1:n,2]#copy sampling dates from any ttree
  l=list(ttree=cons,nam=record[[1]]$ctree$nam)
  class(l)<-'ttree'
  return(l)
}
