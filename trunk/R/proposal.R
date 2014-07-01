#Proposal distribution used in rjMCMC to update the transmission tree
#Returns the proposed tree as well as qr=proposal ratio part of the Metropolis-Hastings ratio
.proposal = function(tree)  {
  which=sample.int(3,1)
  if (which==1) return(.move1(tree))#Add
  if (which==2) return(.move2(tree))#Remove
  if (which==3) return(.move3(tree))#Update
}

#MOVE1: Add a new transmission event
.move1 = function(tree) {
  fathers <- rep(NA, nrow(tree));fathers[tree[ ,2:3] + 1] <- 1:nrow(tree);fathers <- fathers[-1] 
  totbralen=sum(head(tree[,1],-2)-head(tree[fathers,1],-1))
  loc=runif(1)*totbralen
  bra=1
  while (loc>(tree[bra,1]-tree[fathers[bra],1])) {
    loc=loc- (tree[bra,1]-tree[fathers[bra],1])
    bra=bra+1}
  nsam <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
  tree=.treeAdd(tree,tree[fathers[bra],1]+loc,bra,fathers[bra])
  tree <- cbind(tree[ ,1:3],.hostFromFulltree(tree)) 
  ntraeve=sum( tree[ ,2] > 0&tree[ ,3] == 0)
  qr=totbralen/(ntraeve-nsam)
  return(list(tree=tree,qr=qr))
}

#MOVE2: remove a transmission event if possible
.move2 = function(tree) {
  nsam <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
  ntraeve=sum(tree[ ,2]  > 0&tree[ ,3] == 0)
  if (nsam==ntraeve) return(list(tree=tree,qr=1))#Nothing to remove
  
  #Choose a transmission event that can be removed
  host <- tree[ ,4]
  fathers <- rep(NA, nrow(tree));fathers[tree[ ,2:3] + 1] <- 1:nrow(tree);fathers <- fathers[-1] 
  totbralen=sum(head(tree[,1],-2)-head(tree[fathers,1],-1))
  w=0
  while (w==0||w==nrow(tree)||(infector<=nsam && infected<=nsam)){
    w <- sample(which( tree[ ,2] > 0&tree[ ,3] == 0 ) ,1)
    infector <- host[w] 
    infected <- host[tree[w,2]]}
  
  #Remove it
  tree=.treeRem(tree,w,fathers[w])
  tree <- cbind(tree[ ,1:3],.hostFromFulltree(tree)) 
  qr=(ntraeve-nsam)/totbralen
  return(list(tree=tree,qr=qr))
}

#MOVE3: update a transmission event
.move3 = function(tree) {
  nsam <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
  host <- tree[ ,4]
  fathers <- rep(NA, nrow(tree));fathers[tree[ ,2:3] + 1] <- 1:nrow(tree);fathers <- fathers[-1] 
  
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
    
    #Remove the transmission node 
    tree=.treeRem(tree,w,fathers[w])
    fathers <- rep(NA, nrow(tree));fathers[tree[ ,2:3] + 1] <- 1:nrow(tree);fathers <- fathers[-1] 

    #Choose where to add it back on the path from infector to infected 
    islocs <- rep(0, nrow(tree)) 
    path <- infector 
    islocs[path] <- 1 
    while (path < nrow(tree))  { 
      path <- fathers[path] 
      islocs[path] <- 1 
    } 
    path <- infected 
    islocs[path] <- 1 
    while (path < nrow(tree))  { 
      path <- fathers[path] 
      islocs[path] <- 1-islocs[path] 
    } 
    locs <- which(islocs==1) 
    lens <- rep(0,length(locs)) 
    for (i in 1:length(locs)) { 
      lens[i] <- tree[locs[i],1]-tree[fathers[locs[i]],1] 
    } 
    r <- runif(1) * sum(lens) 
    i=1
    while (r>lens[i]) {r<-r-lens[i];i=i+1}
    
    #Add it back
    tree=.treeAdd(tree,tree[locs[i],1]-r,locs[i],fathers[locs[i]])
    tree <- cbind(tree[ ,1:3],.hostFromFulltree(tree)) 
    
  } else {
    #Update of transmission event involving an unsampled individual
    #TODO: add a move proposing to move a star a bit up or down
  }
  
  return(list(tree=tree,qr=1))  
}

#.reordernodes = function(tree) {
#  #Reorder nodes chronologically 
#  nsam <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
#  MySort <- sort(tree[(nsam+1):nrow(tree),1],decreasing=TRUE,index.return = TRUE); ind <- MySort$ix 
#  for (i in (nsam+1):nrow(tree)) for (j in (2:3)) if (tree[i,j] > nsam) tree[i,j] <- nsam + which( ind == tree[i,j]-nsam ) 
#  tree <- tree[c(1:nsam,nsam + ind), ] 
#  tree <- cbind(tree[ ,1:3],.hostFromFulltree(tree)) 
#}

#Add a transmission node onto a tree
.treeAdd = function(tree,age,child,father) {
  nsam <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
  w=nsam+min(which(tree[(nsam+1):nrow(tree),1]<age))
  tree[which(tree[,2]>=w),2]=tree[which(tree[,2]>=w),2]+1
  tree[which(tree[,3]>=w),3]=tree[which(tree[,3]>=w),3]+1
  tree[father,1+which(tree[father,2:3]==child)]=w
  tree=rbind(tree[1:(w-1),],c(age,child,0,0),tree[w:nrow(tree),])
}

#Remove a transmission node from a tree
.treeRem = function(tree,w,father) {
  tree[father,1+which(tree[father,2:3]==w)]=tree[w,2]
  tree[which(tree[,2]>w),2]=tree[which(tree[,2]>w),2]-1
  tree[which(tree[,3]>w),3]=tree[which(tree[,3]>w),3]-1
  tree=tree[-w,]
}

