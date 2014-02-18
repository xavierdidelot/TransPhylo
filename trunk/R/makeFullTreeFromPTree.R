makeFullTreeFromPTree = function(tree)  {
  #Return a non-zero probability phylogenetic+transmission tree for a given phylogenetic tree 
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
  tree[ ,1] <- tree[ ,1]-min(tree[ ,1])
  
  #Reorder nodes chronologically 
  MySort <- sort(tree[seq(n + 1,nrow(tree),1),1],decreasing=TRUE,index.return = TRUE); ind <- MySort$ix 
  for (i in seq(n+1,nrow(tree),1)) { 
    for (j in (2:3)) { 
      if (tree[i,j] > n)  { 
        tree[i,j] <- n + which( ind == tree[i,j]-n ) 
      } 
    } 
  } 
  tree <- tree[c(1:n,n + ind), ] 
  fulltree <- cbind(tree,hostFromFulltree(tree)) 
  return(fulltree)
} 
