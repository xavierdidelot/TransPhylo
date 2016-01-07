#' Build a consensus transmission tree from a MCMC output
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
consTTree = function(record,burnin=0.5)
{
  #Remove burnin
  record=record[round(length(record)*burnin):length(record)]
  m=length(record)
  n=sum(record[[1]]$tree[,2]==0&record[[1]]$tree[,3]==0) 

  #Record partitions in sampled transmission tree
  hash=vector('list',n*m*10)
  a=floor(n*10*runif(n))+1
  for (i in 1) #DEBUGGING: WHEN ONLY ONE TTREE, OUTPUT SHOULD BE SAME
  {
    ttree=ttreeFromFullTree(record[[i]]$tree)
    ttree=t(matrix(c(0,0,0,0,0,1,0,0,1),3,3));n=3
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
      hh=1+(sum(a[c])%%(length(hash)))
      found=0
      while (1) {
        if (is.null(hash[[hh]])) break
        if (all(c==hash[[hh]]$c)) break
        hh=hh+1
        if (hh>length(hash)) hh=1
      }
      if (is.null(hash[[hh]])) {hash[[hh]]$c=c;hash[[hh]]$n=0;hash[[hh]]$w=c();hash[[hh]]$last=0}
      if (hash[[hh]]$last==i) hash[[hh]]$w[1]=hash[[hh]]$w[1]+1
      else {
        hash[[hh]]$last=i
        hash[[hh]]$n=hash[[hh]]$n+1
        hash[[hh]]$w=c(1,hash[[hh]]$w)
      }
    }
    
    #Add bonus leaves
    for (j in 1:n) {
      if ((sum(children[j,]))==1) next
      c=rep(0,n);c[j]=1
      hh=1+(sum(a[c])%%(length(hash)))
      found=0
      while (1) {
        if (is.null(hash[[hh]])) break
        if (all(c==hash[[hh]]$c)) break
        hh=hh+1
        if (hh>length(hash)) hh=1
      }
      if (is.null(hash[[hh]])) {hash[[hh]]$c=c;hash[[hh]]$n=0;hash[[hh]]$w=c();hash[[hh]]$last=0}
      hash[[hh]]$n=hash[[hh]]$n+1
      hash[[hh]]$w=c(0,hash[[hh]]$w)
    }
  }
  
  #Choose partitions to include in consensus transmission tree
  comb=hash[!sapply(hash,is.null)]
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
      if (length(intersect(c1,c2)>0) && length(setdiff(c1,c2))>0 && length(setdiff(c2,c1))>0) {exclude=T;print('excluded');break}
      #if (length(intersect(c1,c2)>0) && length(union(c1,c2))==n) {exclude=T;print('excluded');break}
    }
    if (exclude==F) keep=c(keep,i)
  }
  comb=comb[keep]
  
  parents=rep(NA,length(comb))
  for (i in 1:length(comb)) {
    ci=which(comb[[i]]$c==1)
    bestscore=Inf
    for (j in 1:length(comb)) {
      cj=which(comb[[j]]$c==1)
      if (i==j) next
      if (length(setdiff(ci,cj))>0) next
      if (length(setdiff(cj,ci))<bestscore) {bestscore=length(setdiff(cj,ci));parents[i]=j}
    }
  }
  
  tr=list()
  tr$Nnode=length(comb)-n
  tr$tip.label=as.character(1:n)
  tr$edge=cbind(parents[which(!is.na(parents))],which(!is.na(parents)))
  tr$edge.length=rep(NA,nrow(tr$edge))
  for (i in 1:nrow(tr$edge)) tr$edge.length[i]=median(comb[[tr$edge[i,2]]]$w)
  class(tr)<-'phylo'
  return(list(comb=comb,tr=tr))
}