#Build vector 'host' indicating in which host each node of the fulltree is found 
.hostFromFulltree = function(fulltree)  {
  fathers <- rep(NA, nrow(fulltree)) 
  fathers[fulltree[ ,2] + 1] <- 1:nrow(fulltree) 
  fathers[fulltree[ ,3] + 1] <- 1:nrow(fulltree) 
  fathers <- fathers[-1] 
  host <- rep(0, nrow(fulltree)) 
  nsam <- sum(fulltree[ ,2] == 0&fulltree[ ,3] == 0) 
  for (i in 1:nsam) { 
    j <- i 
    while (1)  { 
      if (host[j]>0) print('Error: two leaves in same host')
      host[j] <- i 
      j <- fathers[j] 
      if (fulltree[j,3] == 0) break 
    } 
  } 
  if (nsam==1) return(host)

  dispo=nsam+1
  f <- which( fulltree[,3] == 0 & fulltree[,2]>0 & (1:nrow(fulltree))<nrow(fulltree) ) #transmission events other than root
  for (i in 1:length(f)) { 
    j <- f[i] 
    tocol <- c() 
    while (fulltree[fathers[j],3]>0&&fathers[j]<nrow(fulltree))  { 
      tocol <- c(tocol,j) 
      j <- fathers[j] 
    } 
    if (host[j]==0) {host[j]=dispo;dispo=dispo+1}
    host[tocol] <- host[j] 
  } 
  return(host)
} 
