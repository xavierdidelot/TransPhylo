#' Infer transmission tree given a phylogenetic tree
#' @param ptree Phylogenetic tree
#' @param mcmcIterations Number of MCMC iterations to run the algorithm for
#' @param startNeg Starting value of within-host coalescent parameter Ne*g
#' @param startR Starting value of basic reproduction number R
#' @param startPi Starting value of sampling proportion pi
#' @param updateNeg Whether of not to update the parameter Ne*g
#' @param updateR Whether or not to update the parameter R
#' @param updatePi Whether or not to update the parameter pi
#' @return posterior sample set of transmission trees
inferTTree = function(ptree,mcmcIterations=1000,startNeg=100/365,startR=1,startPi=0.5,updateNeg=TRUE,updateR=TRUE,updatePi=TRUE) {
  #MCMC algorithm
  neg <- startNeg
  R <- startR
  pi <- startPi
  fulltree <- makeFullTreeFromPTree(ptree);#Starting point 
  record <- vector('list',mcmcIterations)
  pTTree <- probTTree(ttreeFromFullTree(fulltree),R,pi) 
  pPTree <- probPTreeGivenTTree(fulltree,neg) 
  for (i in 1:mcmcIterations) {#Main MCMC loop
    if (i%%100 == 0) print(i) 
    #Record things 
    record[[i]]$tree <- fulltree
    record[[i]]$pTTree <- pTTree 
    record[[i]]$pPTree <- pPTree 
    record[[i]]$neg <- neg 
    record[[i]]$R <- R
    record[[i]]$pi <- pi
    record[[i]]$source <- fulltree[fulltree[which(fulltree[,1]==0),2],4] 
    
    #Metropolis update for transmission tree 
    prop <- .proposal(fulltree) 
    fulltree2 <- prop$tree
    pTTree2 <- probTTree(ttreeFromFullTree(fulltree2),R,pi) 
    pPTree2 <- probPTreeGivenTTree(fulltree2,neg) 
    if (log(runif(1)) < log(prop$qr)+pTTree2 + pPTree2-pTTree-pPTree)  { 
      fulltree <- fulltree2 
      pTTree <- pTTree2 
      pPTree <- pPTree2 
    } 
    
    if (updateNeg) {
      #Metropolis update for Ne*g, assuming Exp(1) prior 
      neg2 <- abs(neg + (runif(1)-0.5)*0.1)
      pPTree2 <- probPTreeGivenTTree(fulltree,neg2) 
      if (log(runif(1)) < pPTree2-pPTree-neg2+neg)  {neg <- neg2;pPTree <- pPTree2} 
    }
    
    if (updateR) {
      #Metropolis update for R, assuming Exp(1) prior 
      R2 <- abs(R + (runif(1)-0.5)*0.1)
      pTTree2 <- probTTree(ttreeFromFullTree(fulltree),R2,pi) 
      if (log(runif(1)) < pTTree2-pTTree-R2+R)  {R <- R2;pTTree <- pTTree2}
    }
    
    if (updatePi) {
      #Metropolis update for pi, assuming Unif(0,1) prior 
      pi2 <- abs(pi + (runif(1)-0.5)*0.1)
      if (pi2>1) pi2=2-pi2
      pTTree2 <- probTTree(ttreeFromFullTree(fulltree),R,pi2) 
      if (log(runif(1)) < pTTree2-pTTree)  {pi <- pi2;pTTree <- pTTree2}       
    }
    
  }#End of main MCMC loop
  
  return(record)
}
