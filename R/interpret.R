#' Returns and/or plot numbers of sampled and unsampled cases over time
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @param numBins Number of time bins to compute and display incident cases
#' @param show.plot Show a plot of incident cases over time with stacked bars
#' @return List with four entries. Time is a vector of the time points. allCases is the average number of cases at each time in the posterior. sampledCases: average number of sampled cases. unsampCases: average number of unsampled cases.
#' @export
getIncidentCases <- function(record,burnin=0.5,numBins=10,show.plot=FALSE) {
  record=record[max(1,round(length(record)*burnin)):length(record)]
  minTime=Inf
  maxTime=-Inf
  for (n in 1:length(record))  {
    #minTime=min(minTime,min(record[[n]]$ctree$ctree[,1]))
    #maxTime=max(maxTime,max(record[[n]]$ctree$ctree[,1]))
    thisTT <- extractTTree(record[[n]]$ctree)$ttree
    minTime=min(minTime,min(thisTT[,1]))
    maxTime=max(maxTime,max(thisTT[,1]))
  }
  breaks=(0:numBins)/numBins*(maxTime-minTime)+minTime
  allcounts=list()
  sampcounts=list()
  unsampcounts=list()
  for (n in 1:length(record))  {
    thisTT <- extractTTree(record[[n]]$ctree)$ttree
    allcounts[[n]]=hist(thisTT[,1],breaks=breaks,plot=FALSE)
    IND=!is.na(thisTT[,2]) # SAMPLED
    sampcounts[[n]]=hist(thisTT[IND,1],breaks=breaks,plot=FALSE) # sampled
    unsampcounts[[n]]=hist(thisTT[!IND,1],breaks=breaks,plot=FALSE) # not sampled
  }
  # find average numbers of sampled, unsampled and total cases in each time bin
  runsum=0*allcounts[[1]]$counts
  runsamp=0*sampcounts[[1]]$counts
  rununsamp=0*unsampcounts[[1]]$counts
  for (n in 1:length(allcounts)) {
    runsum=runsum+allcounts[[n]]$counts
    runsamp=runsamp+sampcounts[[n]]$counts
    rununsamp=rununsamp+unsampcounts[[n]]$counts
  }
  runsum=runsum/length(allcounts)
  runsamp=runsamp/length(sampcounts)
  rununsamp=rununsamp/length(unsampcounts)
  
  res=list(Time=allcounts[[1]]$mids,allCases=runsum, sampledCases=runsamp, unsampCases=rununsamp);
  
  if (show.plot) {
    mydata=t(cbind(res$sampledCases,res$unsampCases))
    barplot(mydata,col=c('blue','red'))
    labs=pretty(c(minTime,maxTime),6)
    axis(1,at = (labs-minTime)/(maxTime-minTime)*1.2*numBins,labels=labs)
    legend("topleft", legend=c('Unsampled','Sampled'), pch=15, col=c('red','blue'))
  }
  return(res)
}

#' Extract and return distribution of infection time of a given sampled case
#' @param record MCMC output produced by inferTTree
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @param k Case whose posterior infection times are to be extracted. Either an integer or a string matching one of the case names in the data
#' @param show.plot Show a barplot of the distribution
#' @return A vector of posterior infection times for case k 
#' @export
getInfectionTimeDist <- function(record,burnin=0.5,k,show.plot=F) {
  record=record[max(1,round(length(record)*burnin)):length(record)]  
  if (is.numeric(k)) 
    mytimes= vapply(1:length(record), function(x) { tt=extractTTree(record[[x]]$ctree); return(tt$ttree[k,1])},FUN.VALUE=1)
  else 
    mytimes= vapply(1:length(record), function(x) { tt=extractTTree(record[[x]]$ctree); ii=which(tt$nam==k); return(tt$ttree[ii,1])},FUN.VALUE=1)
  if (show.plot) hist(mytimes,main='',xlab='',ylab='')
  return(mytimes)
}

#' Extract and return realised generation time distribution
#' @param record MCMC output produced by inferTTree
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @param maxi Maximum generation time to consider
#' @param numBins Number of time bins to compute and display distribution
#' @param show.plot Show a barplot of the distribution
#' @return Vector of times between becoming infected and infecting others (generation times) in the posterior
#' @export
getGenerationTimeDist <-  function(record,burnin=0.5,maxi=2,numBins=20,show.plot=F) 
{ 
  record=record[max(1,round(length(record)*burnin)):length(record)]  
  breaks=c((0:numBins)/numBins*maxi,Inf)
  dist=rep(0,numBins+1)
  for (i in 1:length(record)) {
  tt=extractTTree(record[[i]]$ctree)$ttree
  t1=tt[tt[tt[,3]!=0,3],1]
  t2=tt[tt[,3]!=0,1]
  dist=dist+hist(t2-t1,breaks=breaks,plot=F)$counts
  }
  dist=dist/sum(dist)
  dist=dist[1:numBins]
  meds=(breaks[2:(length(breaks)-1)]+breaks[1:(length(breaks)-2)])/2
  res=cbind(meds,dist)
  colnames(res)<-c('Time','Proba')
  
  if (show.plot) {
    barplot(res[,2])
    labs=pretty(c(0,maxi),6)
    axis(1,at = labs/maxi*1.2*numBins,labels=labs)
  }
  
  return(res)
}

#' Extract and return realised sampling time distribution
#' @param record MCMC output produced by inferTTree
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @param maxi Maximum generation time to consider
#' @param numBins Number of time bins to compute and display distribution
#' @param show.plot Show a barplot of the distribution
#' @return Vector of times between becoming infected and becoming sampled in the posterior
#' @export
getSamplingTimeDist <-  function(record,burnin=0.5,maxi=2,numBins=20,show.plot=F) 
{ 
  record=record[max(1,round(length(record)*burnin)):length(record)]  
  breaks=c((0:numBins)/numBins*maxi,Inf)
  dist=rep(0,numBins+1)
  for (i in 1:length(record)) {
    tt=extractTTree(record[[i]]$ctree)$ttree
    t1=tt[!is.na(tt[,2]),1]
    t2=tt[!is.na(tt[,2]),2]
    dist=dist+hist(t2-t1,breaks=breaks,plot=F)$counts
  }
  dist=dist/sum(dist)
  dist=dist[1:numBins]
  meds=(breaks[2:(length(breaks)-1)]+breaks[1:(length(breaks)-2)])/2
  res=cbind(meds,dist)
  colnames(res)<-c('Time','Proba')
  
  if (show.plot) {
    barplot(res[,2])
    labs=pretty(c(0,maxi),6)
    axis(1,at = labs/maxi*1.2*numBins,labels=labs)
  }
  
  return(res)
}

