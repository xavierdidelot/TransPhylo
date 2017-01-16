#' Calculate the probability of a phylogenetic tree given a transmission tree
#' @param ctree Combined phylogenetic/transmission tree
#' @param neg Within-host coalescent rate
#' @return Probability of phylogeny given transmission tree
#' @export
probPTreeGivenTTree = function(ctree,neg)  {
  if (is.list(ctree)) ctree=ctree$ctree
  n <- max(ctree[,4])
  prob <- 0 
  for (i in (1:n)) { 
    if (length(which(ctree[,4]==i))==1) next
    subtree <- .extractSubtree(ctree,i) 
    prob <- prob + .probSubtree(subtree,neg) 
  } 
  return(prob)
} 

.extractSubtree = function(ctree,w)  {
  #Take all nodes in host 
  ind <- which(ctree[,4] == w)
  #Add father of oldest node 
  ind <- c(ind,which(ctree[,2]==ind[length(ind)]|ctree[,3]==ind[length(ind)]))
  #Create subtree 
  subtree <- matrix(0, length(ind), 2) 
  for (i in 1:length(ind)) { 
    subtree[i,1] <- ctree[ind[i],1]
    subtree[which(ind==ctree[ind[i],2]),2]=i
    subtree[which(ind==ctree[ind[i],3]),2]=i
  } 
  return(subtree)
} 

.probSubtree=function(tab,rate)  {
  #tab(:,1)=times at bottom;tab(:,2)=father;rate=coalescence rate 
  #Return the log-prior probability of a subtree 
  #This is an extension to Eq1 of Drummond et al(2002) Genetics 161:1307-1320 that accounts for condition TMRCA<INCUBATION_PERIOD 
  p <- 0 
  tab[ ,1] <- max(tab[ ,1])-tab[ ,1];#convert times to ages
  isiso <- rep(1, nrow(tab)) 
  isiso[tab[1:(nrow(tab)-1),2]] <- 0 
  ex <- rep(0, nrow(tab))#Keep track of which nodes are active 
  #MySort <- sort(tab[ ,1],index.return = TRUE); ind <- MySort$ix 
  ind=c(1+which(tab[2:nrow(tab),1]<tab[1,1]),1,1+which(tab[2:nrow(tab),1]>tab[1,1]))#The tab matrix should be already ordered, except maybe the first leaf corresponding to a sample
  iso=c(which(tab[,1]<tab[1,1]&isiso),1,which(tab[,1]>tab[1,1]&isiso))#List of leaves in increasing order of age
  cur <- iso[1]#Start with youngest leaf
  while (tab[cur,2] > 0) {ex[cur] <- 1;cur <- tab[cur,2]}#Activate path to root
  ex[length(ex)] <- 1#Activate root
  for (l in iso[-1]) {#For all leaves in increasing order of age
    isanc <- rep(0, nrow(tab))#Ancestors of the current node 
    anc <- l
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
    for (i in ind) {
      if (i == l)  { 
        start <- TRUE 
        curage <- tab[l,1] 
        next
      } 
      if (ex[i]==0) next #Ignore non-existent nodes
      if (start)  { 
        if (!found)  { 
          bra1 <- bra1 + k*(tab[i,1]-curage) 
          if (isanc[i])  found <- TRUE 
        } 
        bra2 <- bra2 + k*(tab[i,1]-curage) 
      } 
      curage <- tab[i,1] 
      if (isiso[i])  k <- k + 1 else k <- k-ex[i] + 1
    } 
    p <- p-log(rate)-bra1/rate-log(1-exp(-bra2/rate))
    if (k != 1) print('errorHere')
    if (bra1 >= bra2)  print('errorThere')
    #Make all ancestors of current node active
    cur <- l
    while (TRUE) {
      if (ex[cur] == 1) {ex[cur] <- 2;break}
      ex[cur] <- 1
      cur <- tab[cur,2]
    } 
  }
  if (min(ex)==0) {print('error: some nodes have not been activated');print(tab)}
  if (length(which(ex==2))!=(length(iso)-1)) print('error: some internal nodes have not been activated twice')
  return(p)
}
#.probSubtree=memoise(.probSubtree0)

