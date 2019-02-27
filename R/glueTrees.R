#Glue together some within-host trees using a transmission tree 
.glueTrees = function(ttree,wtree)  {
  nsam <- length(which(!is.na(ttree[,2])))
  n <- nrow(ttree)
  nh <- n-nsam
  no <- nsam + 1 
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
    tree <- wtree[[i]]
    ni <- nrow(wtree[[i]])/2 
    tree[ ,1] <- tree[ ,1] + ttree[i,1];#Add infection time to all nodes 
    if (!is.na(ttree[i,2])) leaves <- rbind(leaves,tree[1, ]) 
    tree <- tree[(ni+1):nrow(tree), ,drop=FALSE] #keep only internal nodes of current within-host tree
    f <- which( ttree[ ,3] == i )#Infected by current nodes 
    for (j in (1:nrow(tree))) for (k in (2:3)) { 
        if (tree[j,k] == 0)  {next} 
        if (!is.na(ttree[i,2])&&tree[j,k] == 1)  {tree[j,k] <- i}
        else if (tree[j,k] <= ni)  { 
          tree[j,k] <- labh[f[tree[j,k]-!is.na(ttree[i,2])]] 
        } else { tree[j,k] <- lab[i,tree[j,k]-ni] 
        } 
    } 
    intnodes <- rbind(tree,intnodes) 
  } 
  tree <- rbind(leaves,intnodes) 
  ind = order(tree[(nsam+1):nrow(tree),1],decreasing=TRUE)
  for (i in (nsam+1):nrow(tree)) for (j in (2:3)) if (tree[i,j] > nsam) tree[i,j] <- nsam + which( ind == tree[i,j]-nsam ) 
  tree <- tree[c(1:nsam,nsam + ind), ] 
  tree <- cbind(tree,.computeHost(tree)) 
  return(tree)
} 
