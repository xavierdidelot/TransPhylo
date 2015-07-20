#' Infer transmission tree given a phylogenetic tree
#' @param ptree Phylogenetic tree
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation length w
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation length w 
#' @param mcmcIterations Number of MCMC iterations to run the algorithm for
#' @param startNeg Starting value of within-host coalescent parameter Ne*g
#' @param startOff.r Starting value of parameter off.r
#' @param startOff.p Starting value of parameter off.p
#' @param startPi Starting value of sampling proportion pi
#' @param updateNeg Whether of not to update the parameter Ne*g
#' @param updateOff.r Whether or not to update the parameter off.r
#' @param updateOff.p Whether or not to update the parameter off.p
#' @param updatePi Whether or not to update the parameter pi
#' @param datePresent Date when process stops (this can be Inf for fully simulated outbreaks)
#' @return posterior sample set of transmission trees
inferTTree = function(ptree,w.shape=2,w.scale=1,mcmcIterations=1000,startNeg=100/365,startOff.r=1,startOff.p=0.5,startPi=0.5,updateNeg=TRUE,updateOff.r=TRUE,updateOff.p=TRUE,updatePi=TRUE,datePresent=Inf) {
  #print(is.memoised(.getOmegabar))
  memoise::forget(.getOmegabar)
  #print(is.memoised(.probSubtree))
  memoise::forget(.probSubtree)
  #if (testing) {
  #  v=ceiling(nrow(ptree)/2+1):nrow(ptree)
  #  totbralen=-sum(ptree[v,1]-ptree[ptree[v,2],1])-sum(ptree[v,1]-ptree[ptree[v,3],1])
  #  probPTreeGivenTTree = function(fulltree,neg) {return(0)} 
  #  probTTree = function(ttree,pi.r,pi.p,pi,w.shape,w.scale) {#In test mode, the prob of a ttree is just prob of number of sampled cases given pi
  #    nsam=length(which(!is.na(ttree[,2])))
  #    n=nrow(ttree)
  #    return(dbinom(nsam,n,pi,log=TRUE)-log(totbralen)*(n-1)+lfactorial(n-nsam))}
  #}
  #MCMC algorithm
  neg <- startNeg
  off.r <- startOff.r
  off.p <- startOff.p
  pi <- startPi
  fulltree <- makeFullTreeFromPTree(ptree)#Starting point 
  ttree <- ttreeFromFullTree(fulltree)
  record <- vector('list',mcmcIterations)
  pTTree <- probTTree(ttree,off.r,off.p,pi,w.shape,w.scale,datePresent) 
  pPTree <- probPTreeGivenTTree(fulltree,neg) 
  for (i in 1:mcmcIterations) {#Main MCMC loop
    if (i%%100 == 0) message(i) 
    #Record things 
    record[[i]]$tree <- fulltree
    record[[i]]$pTTree <- pTTree 
    record[[i]]$pPTree <- pPTree 
    record[[i]]$neg <- neg 
    record[[i]]$off.r <- off.r
    record[[i]]$off.p <- off.p
    record[[i]]$pi <- pi
    record[[i]]$w.shape <- w.shape
    record[[i]]$w.scale <- w.scale
    record[[i]]$source <- fulltree[fulltree[which(fulltree[,1]==0),2],4] 
    
    #Metropolis update for transmission tree 
    prop <- .proposal(fulltree) 
    fulltree2 <- prop$tree
    ttree2 <- ttreeFromFullTree(fulltree2)
    pTTree2 <- probTTree(ttree2,off.r,off.p,pi,w.shape,w.scale,datePresent) 
    pPTree2 <- probPTreeGivenTTree(fulltree2,neg) 
    if (log(runif(1)) < log(prop$qr)+pTTree2 + pPTree2-pTTree-pPTree)  { 
      fulltree <- fulltree2 
      ttree <- ttree2
      pTTree <- pTTree2 
      pPTree <- pPTree2 
    } 
    
    if (updateNeg) {
      #Metropolis update for Ne*g, assuming Exp(1) prior 
      neg2 <- abs(neg + (runif(1)-0.5)*0.1)
      pPTree2 <- probPTreeGivenTTree(fulltree,neg2) 
      if (log(runif(1)) < pPTree2-pPTree-neg2+neg)  {neg <- neg2;pPTree <- pPTree2} 
    }
    
    if (updateOff.r) {
      #Metropolis update for off.r, assuming Exp(1) prior 
      off.r2 <- abs(off.r + (runif(1)-0.5)*0.5)
      pTTree2 <- probTTree(ttree,off.r2,off.p,pi,w.shape,w.scale,datePresent) 
      if (log(runif(1)) < pTTree2-pTTree-off.r2+off.r)  {off.r <- off.r2;pTTree <- pTTree2}
    }
    
    if (updateOff.p) {
      #Metropolis update for off.p, assuming Unif(0,1) prior 
      off.p2 <- abs(off.p + (runif(1)-0.5)*0.1)
      if (off.p2>1) off.p2=2-off.p2
      pTTree2 <- probTTree(ttree,off.r,off.p2,pi,w.shape,w.scale,datePresent) 
      if (log(runif(1)) < pTTree2-pTTree)  {off.p <- off.p2;pTTree <- pTTree2}
    }

        if (updatePi) {
      #Metropolis update for pi, assuming Unif(0,1) prior 
      pi2 <- abs(pi + (runif(1)-0.5)*0.1)
      if (pi2>1) pi2=2-pi2
      pTTree2 <- probTTree(ttree,off.r,off.p,pi2,w.shape,w.scale,datePresent) 
      if (log(runif(1)) < pTTree2-pTTree)  {pi <- pi2;pTTree <- pTTree2}       
    }
    
  }#End of main MCMC loop
  
  return(record)
}
