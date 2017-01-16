#' Build a matrix of probability of who infected whom from a MCMC output
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
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