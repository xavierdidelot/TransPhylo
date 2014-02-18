#' Extracts phylogenetic tree from a combined phylogenetic/transmission tree 
#' @param tree Combined tree
#' @return phylogenetic tree
ptreeFromFullTree = function(tree)  {
  n <- sum(tree[ ,2] + tree[ ,3] == 0) 
  tra <- n + 1 
  while (tra < nrow(tree))  { 
    if (tree[tra,3] != 0)  { 
      tra <- tra + 1 
      next 
    } 
    t <- tree[ ,2:3] 
    f <- which( t == tra ) 
    t[f] <- tree[tra,2] 
    f <- which( t > tra ) 
    t[f] <- t[f]-1 
    tree[ ,2:3] <- t 
    tree <- tree[setdiff(1:nrow(tree),tra), ] 
    tra <- n + 1 
  } 
  ptree <- tree[1:(nrow(tree)-1),1:(ncol(tree)-1)] 
  return(ptree)
} 
