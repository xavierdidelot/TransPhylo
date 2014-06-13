#' Infer transmission tree given a phylogenetic tree
#' @param ptree Phylogenetic tree
#' @param mcmcIterations Number of MCMC iterations to run the algorithm for
#' @return posterior sample set of transmission trees
inferTTree = function(ptree,mcmcIterations=1000) {
  #MCMC algorithm
  neg <- 100/365
  R <- 1
  fulltree <- makeFullTreeFromPTree(ptree);#Starting point 
  record <- vector('list',mcmcIterations)
  pTTree <- probTTree(ttreeFromFullTree(fulltree),R) 
  pPTree <- probPTreeGivenTTree(fulltree,neg) 
  for (i in 1:mcmcIterations) { 
    if (i%%100 == 0)  { 
      print(i) 
    } 
    #Record things 
    record[[i]]$tree <- fulltree
    record[[i]]$pTTree <- pTTree 
    record[[i]]$pPTree <- pPTree 
    record[[i]]$neg <- neg 
    record[[i]]$source <- fulltree[nrow(fulltree)-1,4] 
    #Metropolis update for transmission tree 
    fulltree2 <- .proposal(fulltree) 
    pTTree2 <- probTTree(ttreeFromFullTree(fulltree2),R) 
    pPTree2 <- probPTreeGivenTTree(fulltree2,neg) 
    if (log(runif(1)) < pTTree2 + pPTree2-pTTree-pPTree)  { 
      fulltree <- fulltree2 
      pTTree <- pTTree2 
      pPTree <- pPTree2 
    } 
    #Metropolis update for Ne*g,assuming Exp(1) prior 
    neg2 <- neg + (runif(1)-0.5)*record[[1]]$neg 
    if (neg2 < 0)  { 
      neg2 <- -neg2 
    } 
    pPTree2 <- probPTreeGivenTTree(fulltree,neg2) 
    if (log(runif(1)) < pPTree2-pPTree-neg2 + neg)  { 
      neg <- neg2 
      pPTree <- pPTree2 
    } 
  } 
  return(record)
}
