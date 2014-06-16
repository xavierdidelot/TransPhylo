#Proposal distribution used in MCMC to update the transmission tree
.proposal = function(fulltree)  {
  n <- sum(fulltree[ ,2] == 0&fulltree[ ,3] == 0) 
  host <- fulltree[ ,4]
  fathers <- seq(0, nrow(fulltree)+1) 
  fathers[fulltree[ ,2] + 1] <- 1:nrow(fulltree) 
  fathers[fulltree[ ,3] + 1] <- 1:nrow(fulltree) 
  fathers <- fathers[-1] 
  
  #Choose a transmission event 
  f <- which( fulltree[ ,2] > 0&fulltree[ ,3] == 0 ) 
  w <- 1 + floor(runif(1) * length(f)) 
  w <- f[w] 
  
  if (w == nrow(fulltree))  { 
    #Move for the transmission to the index case 
    r <- (runif(1)-0.5)/1 
    mini <- -fulltree[nrow(fulltree)-1,1] 
    if (r < mini)  { 
      r <- mini + (mini-r) 
    } 
    fulltree[1:(nrow(fulltree)-1),1] <- fulltree[1:(nrow(fulltree)-1),1] + r 
  } else { 
    #Move for the transmission to someone else 
    infector <- host[w] 
    infected <- host[fulltree[w,2]] 
    #First remove the transmission node 
    f <- fathers[w] 
    if (fulltree[f,2] == w)  { 
      fulltree[f,2] <- fulltree[w,2] 
      fathers[fulltree[w,2]] <- f 
    } else { 
      fulltree[f,3] <- fulltree[w,2] 
      fathers[fulltree[w,2]] <- f 
    } 
    fulltree[w, ] <- 0 
    #Second add it back on the path from infector to infected 
    islocs <- rep(0, nrow(fulltree)) 
    path1 <- infector 
    islocs[path1] <- 1 
    while (path1 < nrow(fulltree))  { 
      path1 <- fathers[path1] 
      islocs[path1] <- 1 
    } 
    path2 <- infected 
    islocs[path2] <- 1 
    while (path2 < nrow(fulltree))  { 
      path2 <- fathers[path2] 
      islocs[path2] <- 1-islocs[path2] 
    } 
    locs <- which(islocs==1) 
    lens <- rep(0,length(locs)) 
    for (i in 1:length(locs)) { 
      lens[i] <- fulltree[locs[i],1]-fulltree[fathers[locs[i]],1] 
    } 
    r <- runif(1) * sum(lens) 
    for (i in 1:length(locs)) { 
      if (r > lens[i])  { 
        r <- r-lens[i] 
      } else { 
        fulltree[w,1] <- fulltree[locs[i],1]-r 
        j <- fathers[locs[i]] 
        if (fulltree[j,2] == locs[i])  { 
          fulltree[j,2] <- w 
        } else { fulltree[j,3] <- w 
        } 
        fulltree[w,2] <- locs[i] 
        break 
      } 
    } 
  } 
  
  #Reorder nodes chronologically 
  #MySort <- sort(fulltree[seq(from = n + 1,to = nrow(fulltree),by = 1),1],decreasing=TRUE,index.return = TRUE); ind <- MySort$ix 
  #for (i in seq(n+1,nrow(fulltree),1)) { 
  #  for (j in 2:3) { 
  #    if (fulltree[i,j] > n)  { 
  #      fulltree[i,j] <- n + which( ind == fulltree[i,j]-n ) 
  #    } 
  #  } 
  #} 
  #fulltree <- fulltree[c(1:n,n+ind), ] 
  res <- cbind(fulltree[ ,1:3],.hostFromFulltree(fulltree)) 
  return(res)
} 
