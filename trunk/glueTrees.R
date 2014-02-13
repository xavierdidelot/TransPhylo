glueTrees = function(ttree,wtree)  {
  #Glue together the within-host trees using the transmission tree 
  n <- nrow(ttree)
  no <- n + 1 
  ma<-0
  for (i in (1:n)) {ma<-max(ma,nrow(wtree[[i]])/2)} 
  lab<-matrix(0,n,ma)
  labh<-rep(0,n)
  for (i in (n:1)) {#Consider all hosts 
    ni <- nrow(wtree[[i]])/2 
    for (j in (1:ni)) { 
      lab[i,j] <- no 
      no <- no + 1 
    } 
    labh[i] <- no-1 
  } 
  leaves <- cbind() 
  intnodes <- cbind() 
  for (i in (1:n)) { 
    tree <- rbind(wtree[[i]])
    ni <- nrow(wtree[[i]])/2 
    tree[ ,1] <- tree[ ,1] + ttree[i,1];#Add infection time to all nodes 
    leaves <- rbind(leaves,tree[1, ]) 
    tree <- rbind(tree[seq(ni + 1,nrow(rbind(tree)),1), ]) 
    f <- which( ttree[ ,3] == i );#Infected by current nodes 
    for (j in (1:nrow(rbind(tree)))) { 
      for (k in (2:3)) { 
        if (tree[j,k] == 0)  { 
          next 
        } 
        if (tree[j,k] == 1)  { 
          tree[j,k] <- i }
        else if (tree[j,k] <= ni)  { 
          tree[j,k] <- labh[f[tree[j,k]-1]] 
        } else { tree[j,k] <- lab[i,tree[j,k]-ni] 
        } 
      } 
    } 
    intnodes <- rbind(tree,intnodes) 
  } 
  fulltree <- rbind(leaves,intnodes) 
  #  fulltree <- cbind(fulltree,hostFromFulltree(fulltree)) 
} 
