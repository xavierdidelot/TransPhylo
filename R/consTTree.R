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
  hash=vector('list',n*m)
  a=floor(n*runif(n))+1
  for (i in 1:m)
  {
    ttree=ttreeFromFullTree(record[[i]]$tree)
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
      c=children[[j]]
      hh=1+mod(sum(a[c]),m*n)
      found=0
      for (k in 1:length(hash[[hh]]))
        if (children[j,]==hash[hh]$c) {found=k;break}
      
    }
  }
  
  #Choose partitions to include in consensus transmission tree
  cons=ttree
  return(cons)
}