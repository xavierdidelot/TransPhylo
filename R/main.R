#' Simulate an outbreak
#' @param neg the within-host effective population size (Ne) times  generation duration (g)
#' @param R the basic reproduction number
#' @return Combined phylogenetic and transmission tree
simulateOutbreak = function(R,neg) {
  #Create a transmission tree with ten individuals
  n<-1
  while (n!=10) {
    ttree<-makeTTree(R)[[1]]
    n<-nrow(ttree)
    if (is.null(ttree)) {n<-0}
  }
  
  #Create a within-host phylogenetic tree for each infected host
  wtree<-vector('list',n)
  for (i in (1:n)) {
    times<-c(ttree[i,2],ttree[which(ttree[,3]==i),1])-ttree[i,1]
    wtree[[i]]<-.withinhost(times,neg)[[1]];
  }
  
  #Glue these trees together
  truth<-.glueTrees(ttree,wtree)
  truth[,1]<-truth[,1]+2005#Epidemic started in 2005
  return(truth)
}  

#' Infer transmission tree given a phylogenetic tree
#' @param ptree Phylogenetic tree
#' @return posterior sample set of transmission trees
mcmc = function(ptree) {
  #MCMC algorithm
  neg <- 100/365
  R <- 1
  fulltree <- makeFullTreeFromPTree(ptree);#Starting point 
  mcmc <- 1000 
  record <- vector('list',mcmc)
  pTTree <- probTTree(ttreeFromFullTree(fulltree),R) 
  pPTree <- probPTreeGivenTTree(fulltree,neg) 
  for (i in 1:mcmc) { 
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
