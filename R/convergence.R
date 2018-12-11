#' Plot MCMC traces
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @export
plotTraces = function(record,burnin=0) {
  par(mfrow=c(2,2))
  record=record[max(1,round(length(record)*burnin)):length(record)]
  plot(sapply(record,function(x) x$pTTree+x$pPTree),ylab='Posterior probability',
       xlab='MCMC iterations',type='l')
  plot(sapply(record,function(x) x$pi),ylab='Sampling proportion pi',
       xlab='MCMC iterations',type='l')
  plot(sapply(record,function(x) x$neg),ylab='Within-host coalescent rate Ne*g',
       xlab='MCMC iterations',type='l')
  plot(sapply(record,function(x) x$off.r),ylab='Basic reproduction number R',
       xlab='MCMC iterations',type='l')
}

#' Convert to coda mcmc format
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @export
convertToCoda = function(record,burnin=0.5) {
  record=record[max(1,round(length(record)*burnin)):length(record)]
  mat=cbind(
  sapply(record,function(x) x$pi),
  sapply(record,function(x) x$neg),
  sapply(record,function(x) x$off.r))
  colnames(mat)<-c('pi','neg','R')
  return(coda::as.mcmc(mat))
}

