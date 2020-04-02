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
#' @param optiStart Method used to optimised colored tree (0=none, 1=slow, 2=fast)
#' @return A minimal non-zero probability phylogenetic+transmission tree, or an optimised version if parameters are provided
#' @export
makeCTreeFromPTree = function(ptree,off.r=NA,off.p=NA,neg=NA,pi=NA,w.shape=NA,w.scale=NA,ws.shape=NA,ws.scale=NA,T=NA,optiStart=0)  {
  nam=ptree$nam
  tree=ptree$ptree
  if (optiStart==0) {
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
    ind=order(tree[seq(n + 1,nrow(tree),1),1],decreasing=TRUE)
    for (i in (n+1):nrow(tree)) for (j in (2:3)) if (tree[i,j] > n) tree[i,j] <- n + which( ind == tree[i,j]-n ) 
    tree <- tree[c(1:n,n + ind), ] 
    tree <- cbind(tree,.computeHost(tree)) 
    l=list(ctree=tree,nam=nam)
    class(l)<-'ctree'
    return(l)
  }
    

  if (optiStart==1) {
    #Slow optimisation method
    n <- ceiling( nrow(tree)/2 ) 
    ft=makeCTreeFromPTree(ptree)
    pTTree <- probTTree((extractTTree(ft))$ttree,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,T) 
    pPTree <- probPTreeGivenTTree(ft$ctree,neg) 
    try=0
    sumbralen=0
    for (i in (n+1):nrow(tree)) for (j in 2:3) sumbralen=sumbralen+tree[tree[i,j],1]-tree[i,1]
    added=0
    maxi=sumbralen/w.shape/w.scale
    while (try<100) {
      try=try+1
      ctree2 <- move1(ft$ctree)$tree
      ctree2=list(ctree=ctree2,nam=nam)
      class(ctree2)<-'ctree'
      pTTree2 <- probTTree((extractTTree(ctree2))$ttree,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,T) 
      pPTree2 <- probPTreeGivenTTree(ctree2$ctree,neg) 
      if (pTTree2 + pPTree2>pTTree+pPTree)  { 
        ft <- ctree2 
        pTTree <- pTTree2 
        pPTree <- pPTree2 
        try=0
        added=added+1
      }
      if (added>maxi) break
    }
    return(ft)
  }
  
  if (optiStart==2) {
    #Fast optimisation method
    n <- ceiling( nrow(tree)/2) 
    ft=makeCTreeFromPTree(ptree)
    tree=ft$ctree[,1:3]
    
    #Modify tree
    for (i in 1:(nrow(tree)-2)) {
      curi=i
      cura=tree[i,1]
      father=which(tree[,2]==i | tree[,3]==i)
      col=which(tree[father,2:3]==i)+1
      datefather=tree[father,1]
      a=rgamma(1,shape = w.shape,scale=w.scale)
      while (cura-a > datefather) {
        tree=rbind(tree,c(cura-a,curi,0))
        tree[father,col]=nrow(tree)
        cura=cura-a
        curi=nrow(tree)
        a=rgamma(1,shape = w.shape,scale=w.scale)
      }
    }
    
    #Reorder nodes chronologically 
    ind=order(tree[seq(n + 1,nrow(tree),1),1],decreasing=TRUE)
    for (i in (n+1):nrow(tree)) for (j in (2:3)) if (tree[i,j] > n) tree[i,j] <- n + which( ind == tree[i,j]-n ) 
    tree <- tree[c(1:n,n + ind), ] 
    tree <- cbind(tree,.computeHost(tree)) 
    l=list(ctree=tree,nam=nam)
    class(l)<-'ctree'
    return(l)
  }

}

