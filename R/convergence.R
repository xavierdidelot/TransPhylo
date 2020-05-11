#' Plot MCMC traces
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @param extend Whether to also show traces of off.r and off.p
#' @return Returns invisibly the first parameter
#' @export
plotTraces = function(record,burnin=0,extend=F) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(2,ifelse(extend,3,2)))
  record=record[max(1,round(length(record)*burnin)):length(record)]
  plot(sapply(record,function(x) x$pTTree+x$pPTree),ylab='Posterior probability',
       xlab='MCMC iterations',type='l')
  plot(sapply(record,function(x) x$pi),ylab='Sampling proportion pi',
       xlab='MCMC iterations',type='l')
  plot(sapply(record,function(x) x$neg),ylab='Within-host coalescent rate Ne*g',
       xlab='MCMC iterations',type='l')
  plot(sapply(record,function(x) x$off.r*x$off.p/(1-x$off.p)),ylab='Basic reproduction number R',
       xlab='MCMC iterations',type='l')
  if (extend) {
    plot(sapply(record,function(x) x$off.r),ylab='off.r',
         xlab='MCMC iterations',type='l')    
    plot(sapply(record,function(x) x$off.p),ylab='off.p',
         xlab='MCMC iterations',type='l')    
  }
  return(invisible(record))
}

#' Convert to coda mcmc format
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return Object of class mcmc from coda package
#' @export
convertToCoda = function(record,burnin=0.5) {
  record=record[max(1,round(length(record)*burnin)):length(record)]
  mat=cbind(
  sapply(record,function(x) x$pi),
  sapply(record,function(x) x$neg),
  sapply(record,function(x) x$off.r),
  sapply(record,function(x) x$off.p))
colnames(mat)<-c('pi','neg','off.r','off.p')
  return(coda::as.mcmc(mat))
}

