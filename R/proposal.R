#Proposal distribution used in MCMC to update the transmission tree
.proposal = function(tree)  {
  n <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
  host <- tree[ ,4]
  fathers <- seq(0, nrow(tree)+1) 
  fathers[tree[ ,2] + 1] <- 1:nrow(tree) 
  fathers[tree[ ,3] + 1] <- 1:nrow(tree) 
  fathers <- fathers[-1] 
  
  #Choose a transmission event 
  f <- which( tree[ ,2] > 0&tree[ ,3] == 0 ) 
  w <- 1 + floor(runif(1) * length(f)) 
  w <- f[w] 
  
  if (w == nrow(tree))  { 
    #Move for the transmission to the index case 
    r <- (runif(1)-0.5)/1 
    mini <- -tree[nrow(tree)-1,1] 
    if (r < mini)  { 
      r <- mini + (mini-r) 
    } 
    tree[1:(nrow(tree)-1),1] <- tree[1:(nrow(tree)-1),1] + r 
  } else { 
    #Move for the transmission to someone else 
    infector <- host[w] 
    infected <- host[tree[w,2]] 
    #First remove the transmission node 
    f <- fathers[w] 
    if (tree[f,2] == w)  { 
      tree[f,2] <- tree[w,2] 
      fathers[tree[w,2]] <- f 
    } else { 
      tree[f,3] <- tree[w,2] 
      fathers[tree[w,2]] <- f 
    } 
    tree[w, ] <- 0 
    #Second add it back on the path from infector to infected 
    islocs <- rep(0, nrow(tree)) 
    path1 <- infector 
    islocs[path1] <- 1 
    while (path1 < nrow(tree))  { 
      path1 <- fathers[path1] 
      islocs[path1] <- 1 
    } 
    path2 <- infected 
    islocs[path2] <- 1 
    while (path2 < nrow(tree))  { 
      path2 <- fathers[path2] 
      islocs[path2] <- 1-islocs[path2] 
    } 
    locs <- which(islocs==1) 
    lens <- rep(0,length(locs)) 
    for (i in 1:length(locs)) { 
      lens[i] <- tree[locs[i],1]-tree[fathers[locs[i]],1] 
    } 
    r <- runif(1) * sum(lens) 
    for (i in 1:length(locs)) { 
      if (r > lens[i])  { 
        r <- r-lens[i] 
      } else { 
        tree[w,1] <- tree[locs[i],1]-r 
        j <- fathers[locs[i]] 
        if (tree[j,2] == locs[i])  { 
          tree[j,2] <- w 
        } else { tree[j,3] <- w 
        } 
        tree[w,2] <- locs[i] 
        break 
      } 
    } 
  } 
  
  #Reorder nodes chronologically 
  MySort <- sort(tree[seq(n + 1,nrow(tree),1),1],decreasing=TRUE,index.return = TRUE); ind <- MySort$ix 
  for (i in (n+1):nrow(tree)) for (j in (2:3)) if (tree[i,j] > n) tree[i,j] <- n + which( ind == tree[i,j]-n ) 
  tree <- tree[c(1:n,n + ind), ] 
  tree <- cbind(tree[ ,1:3],.hostFromFulltree(tree)) 
  return(tree)
} 
