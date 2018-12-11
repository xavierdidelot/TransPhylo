#' Build a matrix of probability of who infected whom from a MCMC output
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return Matrix of probability of who infected whom
#' @export
computeMatWIW = function(record,burnin=0.5)
{
  #Remove burnin
  if (burnin>0) record=record[round(length(record)*burnin):length(record)]
  m=length(record)
  n=sum(record[[1]]$ctree$ctree[,2]==0&record[[1]]$ctree$ctree[,3]==0) #Number of sampled individuals
  mat=matrix(0,n,n)
  colnames(mat)<-record[[1]]$ctree$nam
  rownames(mat)<-record[[1]]$ctree$nam
  
  for (i in 1:length(record))
  {
    ttree=extractTTree(record[[i]]$ctree)$ttree
    infectors=ttree[1:n,3]
    infecteds=1:n
    w=which(infectors==0|infectors>n)
    infectors=infectors[setdiff(1:n,w)]
    infecteds=infecteds[setdiff(1:n,w)]
    mat[cbind(infectors,infecteds)]=mat[cbind(infectors,infecteds)]+1/length(record)
  }
  return(mat)
}

#' Build a matrix indicating for each pairs of individuals how many intermediates there are in the transmission chain
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return Matrix of intermediates in transmission chains between pairs of hosts
#' @export
computeMatTDist = function(record,burnin=0.5)
{
  #Remove burnin
  if (burnin>0) record=record[round(length(record)*burnin):length(record)]
  m=length(record)
  n=sum(record[[1]]$ctree$ctree[,2]==0&record[[1]]$ctree$ctree[,3]==0) #Number of sampled individuals
  mat=matrix(0,n,n)
  colnames(mat)<-record[[1]]$ctree$nam
  rownames(mat)<-record[[1]]$ctree$nam
  
  for (i in 1:length(record))
  {
    ttree=extractTTree(record[[i]]$ctree)$ttree
    for (a in 2:n) for (b in 1:a) {
      aa=a
      bb=b
      count=0
      while (aa!=bb) {
        if (ttree[aa,1]>ttree[bb,1]) aa=ttree[aa,3] else bb=ttree[bb,3]
        count=count+1
      }
      mat[a,b]=mat[a,b]+count/length(record)
      mat[b,a]=mat[a,b]
    }
  }
  return(mat)
}

