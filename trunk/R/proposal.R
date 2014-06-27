#Proposal distribution used in rjMCMC to update the transmission tree
#Returns the proposed tree as well as qr=proposal ratio part of the Metropolis-Hastings ratio
.proposal = function(tree)  {
  which=sample.int(3,1)
  if (which==1) return(.move1(tree))#Add
  if (which==2) return(.move2(tree))#Remove
  if (which==3) return(.move3(tree))#Update
}

.move1 = function(tree) {
  #MOVE1: Add a new transmission event
  fathers<-.makeFathers(tree)
  totbralen=sum(head(tree[,1],-2)-head(tree[fathers,1],-1))
  loc=runif(1)*totbralen
  bra=1
  while (loc>(tree[bra,1]-tree[fathers[bra],1])) {
    loc=loc- (tree[bra,1]-tree[fathers[bra],1])
    bra=bra+1}
  tree=rbind(tree,c(tree[fathers[bra],1]+loc,bra,0,0))
  tree[fathers[bra],1+which(tree[fathers[bra],2:3]==bra)]=nrow(tree)
  tree=.reordernodes(tree)
  return(list(tree=tree,qr=1))
}

.move2 = function(tree) {
  #MOVE2: remove a transmission event if possible
  nsam <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
  ntraeve=sum( tree[ ,2] > 0&tree[ ,3] == 0)
  if (nsam==ntraeve) return(list(tree=tree,qr=1))

  #Choose a transmission event that can be removed
  host <- tree[ ,4]
  fathers<-.makeFathers(tree)
  w=0
  while (w==0||w==nrow(tree)||(infector<=nsam && infected<=nsam)){
  w <- sample(which( tree[ ,2] > 0&tree[ ,3] == 0 ) ,1)
  infector <- host[w] 
  infected <- host[tree[w,2]]}
  
  tree[fathers[w],1+which(tree[fathers[w],2:3]==w)]=tree[w,2]
  tree[which(tree[,2]>w),2]=tree[which(tree[,2]>w),2]-1
  tree[which(tree[,3]>w),3]=tree[which(tree[,3]>w),3]-1
  tree=tree[-w,]
  tree=.reordernodes(tree)
  return(list(tree=tree,qr=1))
}

.move3 = function(tree) {
  #MOVE3: update a transmission event
    nsam <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
    host <- tree[ ,4]
    fathers <- .makeFathers(tree)
    
    #Choose a transmission event
    w <- sample(which( tree[ ,2] > 0&tree[ ,3] == 0 ) ,1)
    infector <- host[w] 
    infected <- host[tree[w,2]]   
    
    if (w == nrow(tree))  { 
     #Update age of transmission to the index case 
      r <- (runif(1)-0.5)/1 
      mini <- -tree[nrow(tree)-1,1] 
      if (r < mini)  { 
        r <- mini + (mini-r) 
      } 
      tree[1:(nrow(tree)-1),1] <- tree[1:(nrow(tree)-1),1] + r 
      
    } else if (infector<=nsam && infected<=nsam) {
       #Update transmission event from sampled case to sampled case
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
      
    } else {
      #Update of transmission event involving an unsampled individual
      #TODO: add a move proposing to move a star a bit up or down
      }
  
  tree=.reordernodes(tree)
  return(list(tree=tree,qr=1))  
}

.reordernodes = function(tree) {
  #Reorder nodes chronologically 
  nsam <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
  MySort <- sort(tree[(nsam+1):nrow(tree),1],decreasing=TRUE,index.return = TRUE); ind <- MySort$ix 
  for (i in (nsam+1):nrow(tree)) for (j in (2:3)) if (tree[i,j] > nsam) tree[i,j] <- nsam + which( ind == tree[i,j]-nsam ) 
  tree <- tree[c(1:nsam,nsam + ind), ] 
  tree <- cbind(tree[ ,1:3],.hostFromFulltree(tree)) 
}

.makeFathers = function(tree) {
fathers <- rep(NA, nrow(tree)) 
fathers[tree[ ,2] + 1] <- 1:nrow(tree) 
fathers[tree[ ,3] + 1] <- 1:nrow(tree) 
fathers <- fathers[-1] 
return(fathers)
}




