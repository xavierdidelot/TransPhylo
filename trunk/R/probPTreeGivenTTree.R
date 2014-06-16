#' Calculate the probability of a phylogenetic tree given a transmission tree
#' @param fulltree Combined phylogenetic/transmission tree
#' @param neg Within-host coalescent rate
#' @return Probability of phylogeny given transmission tree
probPTreeGivenTTree = function(fulltree,neg)  {
  n <- max(fulltree[,4])
  prob <- 0 
  fathers <- rep(0, nrow(fulltree) + 1) 
  fathers[fulltree[ ,2] + 1] <- 1:nrow(fulltree) 
  fathers[fulltree[ ,3] + 1] <- 1:nrow(fulltree) 
  fathers <- fathers[-1] 
  for (i in (1:n)) { 
    subtree <- .extractSubtree(fulltree,i,fathers) 
    prob <- prob + .probSubtree(subtree,neg) 
  } 
  return(prob)
} 

.extractSubtree = function(fulltree,which,fathers)  {
  #Take all nodes in host 
  host <- cbind(fulltree[ ,4]) 
  ind <- 1:nrow(fulltree) 
  ind <- ind[which(host == which)] 
  #Add father of oldest node 
  j <- which.min(fulltree[ind])
  ind <- c(ind,fathers[ind[j]]) 
  #Create subtree 
  subtree <- matrix(0, length(ind), 2) 
  invind <- rep(0, nrow(fulltree)) 
  invind[ind] <- 1:length(ind) 
  for (i in 1:length(ind)) { 
    subtree[i,1] <- fulltree[ind[i],1] 
    if (i < length(ind))  { 
      subtree[i,2] <- invind[fathers[ind[i]]] 
    } 
  } 
  return(subtree)
} 

.probSubtree = function(tab,rate)  {
  #tab(:,1)=age at bottom;tab(:,2)=father;rate=coalescence rate 
  #Return the log-prior probability of a subtree 
  #This is an extension to Eq1 of Drummond et al(2002)Genetics 
  #161:1307-1320 that accounts for condition TMRCA<INCUBATION_PERIOD 
  p <- 0 
  tab[ ,1] <- max(tab[ ,1])-tab[ ,1];#convert times to ages
  isiso <- rep(1, nrow(tab)) 
  isiso[tab[1:(nrow(tab)-1),2]] <- 0 
  iso <- which(isiso==1) 
  ex <- rep(0, nrow(tab));#Keep track of which nodes are active 
  MySort <- sort(tab[ ,1],index.return = TRUE); s <- MySort$x; ind <- MySort$ix 
  cur <- iso[1];while (tab[cur,2] > 0) {ex[cur] <- 1;cur <- tab[cur,2]};ex[length(ex)] <- 1; 
  for (l in (.seqML(2,length(iso),1))) {#For all others leaves in increasing order of age 
    isanc <- rep(0, nrow(tab));#Ancestors of the current node 
    anc <- iso[l] 
    while (tab[anc,2] > 0)  { 
      isanc[anc] <- 1 
      anc <- tab[anc,2] 
    } 
    bra1 <- 0 
    bra2 <- 0 
    start <- FALSE 
    found <- FALSE 
    curage <- 0 
    k <- 0 
    for (i in (1:(length(s)))) { 
      if (ind[i] == iso[l])  { 
        start <- TRUE 
        curage <- tab[iso[l],1] 
      } 
      if (ex[ind[i]]==0) {next;} 
      if (start)  { 
        if (found == FALSE)  { 
          bra1 <- bra1 + k*(tab[ind[i],1]-curage) 
          if (isanc[ind[i]])  { 
            found <- TRUE 
          } 
        } 
        bra2 <- bra2 + k*(tab[ind[i],1]-curage) 
      } 
      curage <- tab[ind[i],1] 
      if (isiso[ind[i]])  { 
        k <- k + 1 
      } else { k <- k-ex[ind[i]] + 1 
      } 
    } 
    p <- p-log(rate)-bra1/rate-log(1-exp(-bra2/rate))
    if (k != 1)  { 
      print('errorHere')
    } 
    if (bra1 >= bra2)  { 
      print('errorThere')
    } 
    cur <- iso[l];while (tab[cur,2] > 0) {if (ex[cur] == 1) {ex[cur] <- 2;break;};ex[cur] <- 1;cur <- tab[cur,2];} 
  }
  return(p)
}
