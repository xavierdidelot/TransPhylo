hostFromFulltree = function(fulltree)  {
  # Here is the return arguments list return(host) 
  #Build vector 'host' indicating in which host each node of the fulltree is found 
  fathers <- rep(0, nrow(fulltree) + 1) 
  fathers[fulltree[ ,2] + 1] <- 1:nrow(fulltree) 
  fathers[fulltree[ ,3] + 1] <- 1:nrow(fulltree) 
  fathers <- fathers[seq(2,length(fathers),1)] 
  host <- rep(0, nrow(fulltree)) 
  n <- sum(cbind(fulltree[ ,2]) == 0&cbind(fulltree[ ,3]) == 0) 
  for (i in (1:n)) { 
    j <- i 
    while (1)  { 
      host[j] <- i 
      j <- fathers[j] 
      if (fulltree[j,3] == 0)  { 
        break 
      } 
    } 
  } 
  f <- n + which( fulltree[seq(n + 1,nrow(fulltree)-1,1),3] == 0 ); 
  for (i in (seq(1,length(f),1))) { 
    j <- f[i] 
    tocol <- c() 
    while (host[j] == 0)  { 
      tocol <- c(tocol,j) 
      j <- fathers[j] 
    } 
    host[tocol] <- host[j] 
  } 
  return(host)
} 
