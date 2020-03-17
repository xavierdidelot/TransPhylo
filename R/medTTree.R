#' Return the medoid from a MCMC output
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return The index of the medoid
#' @export 
medTTree = function(record,burnin=0.5)
{
  #Remove burnin
  burned=round(length(record)*burnin)
  record=record[(burned+1):length(record)]
  m=length(record)
  if (m==1) return(burned+1)
  n=sum(record[[1]]$ctree$ctree[,2]==0&record[[1]]$ctree$ctree[,3]==0) #Number of sampled individuals

  #Compute tinf and nhosts
  tinf=matrix(0,m,n)
  nhosts=rep(0,m)
  for (i in 1:m)
  {
    ttree=extractTTree(record[[i]]$ctree)$ttree
    tinf[i,]=ttree[1:n,1]
    nhosts[i]=nrow(ttree)
  }
  
  #Select trees with nhosts in IQR
  ws=which(nhosts>=quantile(nhosts,1/4)&nhosts<=quantile(nhosts,3/4))
  m=length(ws)

  #Compute distances based on tinf
  d=matrix(NA,m,m)
  for (i in 1:(m-1)) for (j in (i+1):m) {
    dist=sum((tinf[ws[i],]-tinf[ws[j],])^2)
    d[i,j]=dist
    d[j,i]=dist
  }
  
  #Find medoid
  c=colSums(d,na.rm = T)
  w=ws[which(c==min(c))]
  w=w[1]
  return(record[[w]]$ctree)
}
