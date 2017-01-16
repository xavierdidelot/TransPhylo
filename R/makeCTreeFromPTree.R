#' Create a transmission tree compatible with the provided phylogenetic tree
#' @param ptree Phylogenetic tree
#' @param off.r First parameter of the negative binomial distribution for offspring number
#' @param off.p Second parameter of the negative binomial distribution for offspring number
#' @param neg the within-host effective population size (Ne) timesgeneration duration (g)
#' @param pi probability of sampling an infected individual
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time 
#' @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time
#' @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time 
#' @param T Date when process stops (this can be Inf for fully simulated outbreaks)
#' @param allowTransPostSamp Whether or not to allow transmission after sampling of a host
#' @return A minimal non-zero probability phylogenetic+transmission tree, or an optimised version if parameters are provided
#' @export
makeCtreeFromPTree = function(ptree,off.r=NA,off.p=NA,neg=NA,pi=NA,w.shape=NA,w.scale=NA,ws.shape=NA,ws.scale=NA,T=NA,allowTransPostSamp=NA)  {
  nam=ptree$nam
  tree=ptree$ptree
  if (is.na(off.r)) {
    #Simple minimal method
    n <- ceiling( nrow(tree)/2 ) 
    tree <- rbind(tree,matrix(0, n, 3)) 
    tree[nrow(tree),1] <- min(tree[1:( 2*n-1 ),1])-1
    tree[nrow(tree),2] <- 2*n-1 
    source <- which.max(tree[1:n,1])
    notsource <- setdiff(1:n,source) 
    i2 <- 0 
    for (i in (notsource)) { 
      i2 <- i2 + 1 
      f <- which( tree[ ,2] == i|tree[ ,3] == i ) 
      if (tree[f,2] == i)  { 
        tree[f,2] <- 2*n-1 + i2 
      } else { 
        tree[f,3] <- 2*n-1 + i2 
      } 
      tree[2*n-1 + i2,2] <- i 
      tree[2*n-1 + i2,1] <- (tree[f,1] + tree[i,1])/2 
    } 
    
    #Reorder nodes chronologically 
    MySort <- sort(tree[seq(n + 1,nrow(tree),1),1],decreasing=TRUE,index.return = TRUE); ind <- MySort$ix 
    for (i in (n+1):nrow(tree)) for (j in (2:3)) if (tree[i,j] > n) tree[i,j] <- n + which( ind == tree[i,j]-n ) 
    tree <- tree[c(1:n,n + ind), ] 
    tree <- cbind(tree,.computeHost(tree)) 
    return(list(ctree=tree,nam=nam))
    
  } else {

    #Optimisation method
    n <- ceiling( nrow(tree)/2 ) 
    ft=makeCtreeFromPTree(ptree)
    pTTree <- probTTree(extractTTree(ft),off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,T,allowTransPostSamp) 
    pPTree <- probPTreeGivenTTree(ft,neg) 
    try=0
    while (try<100) {
      try=try+1
      ctree2 <- .move1(ft$ctree)$tree
      ctree2=list(ctree=ctree2,nam=nam)
      pTTree2 <- probTTree(extractTTree(ctree2),off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,T,allowTransPostSamp) 
      pPTree2 <- probPTreeGivenTTree(ctree2,neg) 
      if (pTTree2 + pPTree2>pTTree+pPTree)  { 
        ft <- ctree2 
        pTTree <- pTTree2 
        pPTree <- pPTree2 
        try=0
      }
    }
    return(ft)
  }
}

