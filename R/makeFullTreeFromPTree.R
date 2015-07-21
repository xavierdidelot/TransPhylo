#' Create a transmission tree compatible with the provided phylogenetic tree
#' @param tree Phylogenetic tree
#' @param off.r First parameter of the negative binomial distribution for offspring number
#' @param off.p Second parameter of the negative binomial distribution for offspring number
#' @param pi probability of sampling an infected individual
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation length w
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation length w 
#' @param T Date when process stops (this can be Inf for fully simulated outbreaks)
#' @return A minimal non-zero probability phylogenetic+transmission tree, or an optimised version if parameters are provided
makeFullTreeFromPTree = function(tree,off.r=NA,off.p=NA,neg=NA,pi=NA,w.shape=NA,w.scale=NA,T=NA)  {

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
    #  tree[ ,1] <- tree[ ,1]-min(tree[ ,1])
    #  tree[ ,1] <- tree[ ,1]+dateLastSample-max(tree[ ,1])
    
    #Reorder nodes chronologically 
    MySort <- sort(tree[seq(n + 1,nrow(tree),1),1],decreasing=TRUE,index.return = TRUE); ind <- MySort$ix 
    for (i in (n+1):nrow(tree)) for (j in (2:3)) if (tree[i,j] > n) tree[i,j] <- n + which( ind == tree[i,j]-n ) 
    tree <- tree[c(1:n,n + ind), ] 
    tree <- cbind(tree,.hostFromFulltree(tree)) 
    return(tree)
    
  } else {

    #Optimisation method
    n <- ceiling( nrow(tree)/2 ) 
    ft=makeFullTreeFromPTree(tree)
    pTTree <- probTTree(ttreeFromFullTree(ft),off.r,off.p,pi,w.shape,w.scale,T) 
    pPTree <- probPTreeGivenTTree(ft,neg) 
    try=0
    while (try<100) {
      try=try+1
      fulltree2 <- .move1(ft)$tree
      pTTree2 <- probTTree(ttreeFromFullTree(fulltree2),off.r,off.p,pi,w.shape,w.scale,T) 
      pPTree2 <- probPTreeGivenTTree(fulltree2,neg) 
      if (pTTree2 + pPTree2>pTTree+pPTree)  { 
        ft <- fulltree2 
        pTTree <- pTTree2 
        pPTree <- pPTree2 
        try=0
      }
    }
    return(ft)
  }
}

