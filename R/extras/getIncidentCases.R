#' @param record MCMC output produced by inferTTree
#' @param skipnum Number of record entries to skip in obtaining posterior incident cases, ie skipnum=10 uses the 1st, 11th, 21st ... elements of the record
#' @param numBins Number of time bins to compute and display incident cases
#' @param show.plot Show a plot of incident cases over time using ggplot with stacked bars
#' @return List with four entries. Time is a vector of the time points. allCases is the average number of cases at each time in the posterior. sampledCases: average number of sampled cases. unsampCases: average number of unsampled cases.
getIncidentCases <- function(record,skipnum=1,numBins=10,show.plot=FALSE) {
require(ggplot2)
    thisTT <- extractTTree(record[[length(record)]]$ctree)$ttree;
    breaks=hist(thisTT[,1],numBins,plot=FALSE)$breaks; breaks=c(min(breaks)-3,breaks);
    # avoids error in hist where some events are not counted
    # note using 3 here because finding the minimum time in all the record entries takes too long. First bin is therefore not optimal. 

  allcounts=list(); sampcounts=list(); unsampcounts=list()
    whichones=seq(from=1,to=length(record),by=skipnum); # entries to use
    for (n in 1:length(whichones))  {
    thisTT <- extractTTree(record[[n]]$ctree)$ttree;
    allcounts[[n]]=hist(thisTT[,1],breaks=breaks,plot=FALSE)
    IND=!is.na(thisTT[,2]) # SAMPLED
    sampcounts[[n]]=hist(thisTT[IND,1],breaks=breaks,plot=FALSE) # sampled
    unsampcounts[[n]]=hist(thisTT[!IND,1],breaks=breaks,plot=FALSE) # not sampled
}
# find average numbers of sampled, unsampled and total cases in each time bin
runsum=0*allcounts[[1]]$counts; runsamp=0*sampcounts[[1]]$counts; rununsamp=0*unsampcounts[[1]]$counts;  # initialize
   for (n in 1:length(allcounts)) {
runsum=runsum+allcounts[[n]]$counts
runsamp=runsamp+sampcounts[[n]]$counts
rununsamp=rununsamp+unsampcounts[[n]]$counts
}
    runsum=runsum/length(allcounts)
    runsamp=runsamp/length(sampcounts)
    rununsamp=rununsamp/length(unsampcounts)

test=list(Time=allcounts[[1]]$mids,allCases=runsum, sampledCases=runsamp, unsampCases=rununsamp);
#ggplot call 
   if (show.plot) {
    mydata=data.frame(Time=c(test$Time,test$Time),Cases=c(test$sampledCases,test$unsampCases), IsSampled=c(rep("Yes",length(test$sampledCases)),rep("No",length(test$unsampCases))))
ggplot(data=mydata)+geom_bar(aes(Time,Cases,fill=IsSampled),stat="identity")
}
    return(test)
}
