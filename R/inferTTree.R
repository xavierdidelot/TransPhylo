#' Infer transmission tree given a phylogenetic tree
#' @param ptree Phylogenetic tree
#' @param w.shape Shape parameter of the Gamma distribution representing the generation time
#' @param w.scale Scale parameter of the Gamma distribution representing the generation time 
#' @param ws.shape Shape parameter of the Gamma distribution representing the sampling time
#' @param ws.scale Scale parameter of the Gamma distribution representing the sampling time 
#' @param w.mean Mean of the Gamma distribution representing the generation time
#' @param w.std Std of the Gamma distribution representing the generation time 
#' @param ws.mean Mean of the Gamma distribution representing the sampling time
#' @param ws.std Std of the Gamma distribution representing the sampling time 
#' @param mcmcIterations Number of MCMC iterations to run the algorithm for
#' @param thinning MCMC thinning interval between two sampled iterations
#' @param startNeg Starting value of within-host coalescent parameter Ne*g
#' @param startOff.r Starting value of parameter off.r
#' @param startOff.p Starting value of parameter off.p
#' @param startPi Starting value of sampling proportion pi
#' @param updateNeg Whether of not to update the parameter Ne*g
#' @param updateOff.r Whether or not to update the parameter off.r
#' @param updateOff.p Whether or not to update the parameter off.p
#' @param updatePi Whether or not to update the parameter pi
#' @param startCTree Optional combined tree to start from
#' @param updateTTree Whether or not to update the transmission tree
#' @param optiStart Type of optimisation to apply to MCMC start point (0=none, 1=slow, 2=fast)
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @param verbose Whether or not to use verbose mode (default is false)
#' @return posterior sample set of transmission trees
#' @examples 
#' inferTTree(ptreeFromPhylo(ape::rtree(5),2020),mcmcIterations=100)
#' @export
inferTTree = function(ptree, w.shape=2, w.scale=1, ws.shape=NA, ws.scale=NA, 
                      w.mean=NA,w.std=NA,ws.mean=NA,ws.std=NA,mcmcIterations=1000,
                      thinning=1, startNeg=100/365, startOff.r=1, startOff.p=0.5, startPi=0.5, updateNeg=TRUE,
                      updateOff.r=TRUE, updateOff.p=FALSE, updatePi=TRUE, startCTree=NA, updateTTree=TRUE,
                      optiStart=2, dateT=Inf,verbose=F) {

  ptree$ptree[,1]=ptree$ptree[,1]+runif(nrow(ptree$ptree))*1e-10#Ensure that all leaves have unique times
  if (dateT<dateLastSample(ptree)) stop('The parameter dateT cannot be smaller than the date of last sample')
  for (i in (ceiling(nrow(ptree$ptree)/2)+1):nrow(ptree$ptree)) for (j in 2:3) 
    if (ptree$ptree[ptree$ptree[i,j],1]-ptree$ptree[i,1]<0) 
      stop("The phylogenetic tree contains negative branch lengths!")
  
  if (!is.na( w.mean)&&!is.na( w.std)) { w.shape= w.mean^2/ w.std^2; w.scale= w.std^2/ w.mean}
  if (!is.na(ws.mean)&&!is.na(ws.std)) {ws.shape=ws.mean^2/ws.std^2;ws.scale=ws.std^2/ws.mean}
  if (is.na(ws.shape)) ws.shape=w.shape
  if (is.na(ws.scale)) ws.scale=w.scale
  
  #MCMC algorithm
  neg <- startNeg
  off.r <- startOff.r
  off.p <- startOff.p
  pi <- startPi
  if (is.na(sum(startCTree))) ctree <- makeCTreeFromPTree(ptree,off.r,off.p,neg,pi,w.shape,w.scale,ws.shape,ws.scale,dateT,optiStart)#Starting point 
  else ctree<-startCTree
  ttree <- extractTTree(ctree)
  record <- vector('list',mcmcIterations/thinning)
  pTTree <- probTTree(ttree$ttree,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT) 
  pPTree <- probPTreeGivenTTree(ctree$ctree,neg) 
  if (verbose==F) pb <- utils::txtProgressBar(min=0,max=mcmcIterations,style = 3)
  
  for (i in 1:mcmcIterations) {#Main MCMC loop
    if (i%%thinning == 0) {
      #Record things 
      if (verbose==F) utils::setTxtProgressBar(pb, i)
      if (verbose==T) message(sprintf('it=%d,neg=%f,off.r=%f,off.p=%f,pi=%f,Prior=%e,Likelihood=%e,nind=%d',i,neg,off.r,off.p,pi,pTTree,pPTree,nrow(ttree$ttree)))
      record[[i/thinning]]$ctree <- ctree
      record[[i/thinning]]$pTTree <- pTTree 
      record[[i/thinning]]$pPTree <- pPTree 
      record[[i/thinning]]$neg <- neg 
      record[[i/thinning]]$off.r <- off.r
      record[[i/thinning]]$off.p <- off.p
      record[[i/thinning]]$pi <- pi
      record[[i/thinning]]$w.shape <- w.shape
      record[[i/thinning]]$w.scale <- w.scale
      record[[i/thinning]]$ws.shape <- ws.shape
      record[[i/thinning]]$ws.scale <- ws.scale
      record[[i/thinning]]$source <- ctree$ctree[ctree$ctree[which(ctree$ctree[,4]==0),2],4]
      if (record[[i/thinning]]$source<=length(ctree$nam)) record[[i/thinning]]$source=ctree$nam[record[[i/thinning]]$source] else record[[i/thinning]]$source='Unsampled'
    }
    
    if (updateTTree) {
    #Metropolis update for transmission tree 
    if (verbose) message("Proposing ttree update")
    prop <- proposal(ctree$ctree) 
    ctree2 <- list(ctree=prop$tree,nam=ctree$nam)
    class(ctree2)<-'ctree'
    ttree2 <- extractTTree(ctree2)
    pTTree2 <- probTTree(ttree2$ttree,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT) 
    #pPTree2 <- probPTreeGivenTTree(ctree2$ctree,neg) 
    pPTreeDiff <- probPTreeGivenTTree(ctree2$ctree,neg,prop$new)-probPTreeGivenTTree(ctree$ctree,neg,prop$old)
    if (log(runif(1)) < log(prop$qr)+pTTree2+pPTreeDiff-pTTree)  { 
      ctree <- ctree2 
      ttree <- ttree2
      pTTree <- pTTree2 
      pPTree <- pPTree+pPTreeDiff
    } 
    }
    
    if (updateNeg) {
      #Metropolis update for Ne*g, assuming Exp(1) prior 
      neg2 <- abs(neg + (runif(1)-0.5)*0.5)
      if (verbose) message(sprintf("Proposing Ne*g update %f->%f",neg,neg2))
      pPTree2 <- probPTreeGivenTTree(ctree$ctree,neg2) 
      if (log(runif(1)) < pPTree2-pPTree-neg2+neg)  {neg <- neg2;pPTree <- pPTree2} 
    }
    
    if (updateOff.r) {
      #Metropolis update for off.r, assuming Exp(1) prior 
      off.r2 <- abs(off.r + (runif(1)-0.5)*0.5)
      if (verbose) message(sprintf("Proposing off.r update %f->%f",off.r,off.r2))
      pTTree2 <- probTTree(ttree$ttree,off.r2,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT) 
      if (log(runif(1)) < pTTree2-pTTree-off.r2+off.r)  {off.r <- off.r2;pTTree <- pTTree2}
    }
    
    if (updateOff.p) {
      #Metropolis update for off.p, assuming Unif(0,1) prior 
      off.p2 <- abs(off.p + (runif(1)-0.5)*0.1)
      if (off.p2>1) off.p2=2-off.p2
      if (verbose) message(sprintf("Proposing off.p update %f->%f",off.p,off.p2))
      pTTree2 <- probTTree(ttree$ttree,off.r,off.p2,pi,w.shape,w.scale,ws.shape,ws.scale,dateT) 
      if (log(runif(1)) < pTTree2-pTTree)  {off.p <- off.p2;pTTree <- pTTree2}
    }

        if (updatePi) {
      #Metropolis update for pi, assuming Unif(0.01,1) prior 
      pi2 <- pi + (runif(1)-0.5)*0.1
      if (pi2<0.01) pi2=0.02-pi2
      if (pi2>1) pi2=2-pi2
      if (verbose) message(sprintf("Proposing pi update %f->%f",pi,pi2))
      pTTree2 <- probTTree(ttree$ttree,off.r,off.p,pi2,w.shape,w.scale,ws.shape,ws.scale,dateT) 
      if (log(runif(1)) < pTTree2-pTTree)  {pi <- pi2;pTTree <- pTTree2}       
    }
    
  }#End of main MCMC loop
  
  class(record)<-'resTransPhylo'
  return(record)
}
