#' Extracts phylogenetic tree from a combined phylogenetic/transmission tree 
#' @param ctree Combined tree
#' @return phylogenetic tree
#' @examples
#' extractPTree(simulateOutbreak())
#' @export
extractPTree = function(ctree)  {
  tree=ctree$ctree
  nam=ctree$nam
  n <- sum(tree[ ,2] + tree[ ,3] == 0) 
  tra <- n + 1 
  while (tra <= nrow(tree))  { #Removing transmission events from the tree one by one
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
    tree <- tree[-tra, , drop=FALSE] 
    tra <- n + 1 
  } 
  ptree <- tree[,1:(ncol(tree)-1),drop=FALSE] 
  l=list(ptree=ptree,nam=nam)
  class(l)<-'ptree'
  return(l)
} 
